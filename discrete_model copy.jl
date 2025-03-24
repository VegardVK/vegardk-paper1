using DataFrames
using XLSX
using JuMP
using CPLEX


include("C:/Users/vegardvk/vscodeProjects/bernstein/get_hydro_data.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/helper_functions.jl")

function read_input_data(include_first_stage)
    # Read input data
    input_parameters = (
        wind_df = DataFrame(XLSX.readtable("output/wind_ts_data.xlsx", "Sheet1", infer_eltypes=true)),
        hydro_df = DataFrame(XLSX.readtable("output/hydro_data.xlsx", "Sheet1", infer_eltypes=true)),
        inflow_df = DataFrame(XLSX.readtable("output/inflow_data.xlsx", "Sheet1", infer_eltypes=true)),
        load_df = DataFrame(XLSX.readtable("output/load_data.xlsx", "Sheet1", infer_eltypes=true)),
        plant_df = DataFrame(XLSX.readtable("output/plant_data.xlsx", "Sheet1", infer_eltypes=true)),
        line_df = DataFrame(XLSX.readtable("output/line_data.xlsx", "Sheet1", infer_eltypes=true)),
    )
    @unpack plant_df, load_df, line_df = input_parameters
    

    S_input = unique(load_df.scenario)
    if !(include_first_stage)
        S_input = S_input[2:end]
    end
    area_grouped_plants = groupby(plant_df, :area)
    input_sets = (
        A = unique(plant_df.area),
        P = unique(plant_df.plant_id),
        T = unique(load_df.timestep),
        L = unique(line_df.line_id),
        S = S_input,
        P_w = unique(plant_df[plant_df.fuel_type .== "Wind", :plant_id]),
        P_t = unique(plant_df[plant_df.fuel_type .== "Thermal", :plant_id]),
        P_h = unique(plant_df[plant_df.fuel_type .== "Hydro", :plant_id]),
        P_a = Dict(g.area[1] => unique(g.plant_id) for g in area_grouped_plants),
        # Dicts with the lines that go in/out of area a
        L_in = Dict(a => unique(line_df[line_df.area_to .== a, :line_id]) for a in A), 
        L_out =  Dict(a => unique(line_df[line_df.area_from .== a, :line_id]) for a in A),
        line_pairs = Dict(pair => filter(row -> row.line_pair == pair, eachrow(line_df)).line_id for pair in unique(line_df.line_pair)),
        L_P = unique(line_df.line_pair),

    )
    return input_parameters, input_sets
    # return wind_df, hydro_df, inflow_df, load_df, plant_df, line_df, A, P, T, L, P_w, P_t, P_h, P_a, L_in, L_out, I_disch, I_spill, I_bypass, C_shedding, C_dumping, C_startup, C_spill, C_bypass, S, L_P, line_pairs
end


function define_and_solve_model(include_first_stage=true)
    input_parameters, input_sets = read_input_data(include_first_stage)
    @unpack wind_df, hydro_df, inflow_df, load_df, plant_df, line_df = input_parameters
    @unpack A, P, T, L, P_w, P_t, P_h, P_a, L_in, L_out, S, L_P, line_pairs = input_sets
    C_shedding, C_dumping, C_startup, C_spill, C_bypass = get_cost_parameters()
    I_disch, I_spill, I_bypass = find_connected_plants(hydro_df)

    Δt = 24/T[end]
    π = 1/S[end]

    # display(load_df)
    model = Model()
    set_optimizer(model, CPLEX.Optimizer)
    @variable(model, production[s in S, p in P, t in T] >= 0)
    @variable(model, up_activation[s in S, p in P, t in T] >= 0)
    @variable(model, down_activation[s in S, p in P, t in T] >= 0)

    @variable(model, 0 ≤ load_shedding[s in S, a in A, t in T] ≤ load_df[(load_df.area .== a) .& (load_df.timestep .== t) .& (load_df.scenario .== s), :load][1])
    @variable(model, power_dumping[S in S, a in A, t in T] ≥ 0) 
    @variable(model, transmission[s in S, l in L, t in T] ≥ 0)


    # Unit commitment without binary variables
    # @variable(model, 0 ≤ x1[p in P, t in T] ≤ 1)
    # @variable(model, 0 ≤ x2[p in P, t in T] ≤ 1)
    # @variable(model, startup[p in P, t in T] ≥ 0)
    # @constraint(model, production_linking[p in P, t in T], production[p, t] == x1[p, t] * p_dict["gen_lb"][p] + x2[p, t] * (p_dict["gen_ub"][p] - p_dict["gen_lb"][p]))
    # @constraint(model, start_before_prod[p in P, t in T], x1[p, t] ≥ x2[p, t])
    # @constraint(model, startup_count[p in P, t in T[2:end]], x1[p, t] - x1[p, t-1]  == startup[p, t])


    #second stage production constraints
    @constraint(model, second_stage_gen_ub[s in S[2:end], p in union(P_t, P_h), t in T], production[s, p, t] == production[0, p, t] + up_activation[s, p, t] - down_activation[s, p, t])

    @constraint(model, wind_prod[s in S, p in P_w, t in T], production[s, p, t] ≤ wind_df[(wind_df.plant_id .== p) .& (wind_df.timestep .== t) .& (wind_df.scenario .== 0), :wind_power][1])
    # @constraint(model, wind_prod[p in P_w, t in T], production[0, p, t] <= wind_df[(wind_df.plant_id .== p) .& (wind_df.timestep .== t) .& (wind_df.scenario .== 0), :wind_power][1])  
    # @constraint(model, second_stage_wind[s in S[2:end], p in P_w, t in T], production[s, p, t] <= wind_df[(wind_df.plant_id .== p) .& (wind_df.timestep .== t) .& (wind_df.scenario .== s), :wind_power][1])

    #Hydropower variables
    @variable(model, flow_disch[s in S, p in P_h, t in T] ≥ 0) # Antar at alle moduler bare har ett discharge-segment
    @variable(model, flow_bypass[s in S, p in P_h, t in T] ≥ 0)
    @variable(model, flow_spill[s in S, p in P_h, t in T] ≥ 0)
    @variable(model, total_flow_in[s in S, p in P_h, t in T] ≥ 0)
    @variable(model, total_flow_out[s in S, p in P_h, t in T] ≥ 0)
    @variable(model, volume[s in S, p in P_h, t in 0:T[end]] ≥ 0)

    #Hydropower constraints
    @constraint(model, controlled_inflow[s in S, p in P_h, t in T], total_flow_in[s, p, t] == sum(flow_disch[s, i, t] for i in I_disch[p]) 
                                                                                + sum(flow_bypass[s, i, t] for i in I_bypass[p]) 
                                                                                + sum(flow_spill[s, i, t] for i in I_spill[p])
                                                                                + inflow_df[(inflow_df.plant_id .== p) .& (inflow_df.timestep .== t), :inflow][1])

    @constraint(model, controlled_outflow[s in S, p in P_h, t in T], total_flow_out[s, p, t] == flow_disch[s, p, t] + flow_bypass[s, p, t] + flow_spill[s, p, t])
    @constraint(model, starting_reservoir_constraint[s in S, p in P_h], volume[s, p, 0] == hydro_df[hydro_df.plant_id .== p, :starting_reservoir][1])
    @constraint(model, reservoir_balance[s in S, p in P_h, t in 1:(T[end])], volume[s, p, t] - volume[s, p, t-1] == Δt * (total_flow_in[s, p, t] - total_flow_out[s, p, t]))
    @constraint(model, hydro_production[s in S, p in P_h, t in T], production[s, p, t] == hydro_df[hydro_df.plant_id .== p, :enekv][1] * flow_disch[s, p, t] * 3.6) # kwh/M3 * M3/s * 3600 s/h * 1/1000 MW/kW
    @constraint(model, vol_ub[s in S, p in P_h, t in T], volume[s, p, t] ≤ hydro_df[hydro_df.plant_id .== p, :kap_mag][1])
    @constraint(model, bypass_ub[s in S, p in P_h, t in T], flow_bypass[s, p, t] ≤ hydro_df[hydro_df.plant_id .== p, :kap_forb][1])
    @constraint(model, prod_ub[s in S, p in P_h, t in T], flow_disch[s, p, t] ≤ hydro_df[hydro_df.plant_id .== p, :kap_gen_m3s][1])
    @constraint(model, spill_ub[s in S, p in P_h, t in T], flow_spill[s, p, t] ≤ hydro_df[hydro_df.plant_id .== p, :kap_spill][1])

    # Bare transmisjon på én av linjene mellom to områder
    @variable(model, line_active[s in S, l in L, t in T], Bin)
    @constraint(model, one_directional_flow[s in S, pair in L_P, t in T], line_active[s, line_pairs[pair][1], t] + line_active[s, line_pairs[pair][2], t] .== 1)
    @constraint(model, transmission_ub[s in S, l in L, t in T], transmission[s, l, t] ≤ line_active[s, l, t] * line_df[line_df.line_id .== l, :capacity][1])

    if include_first_stage
        @variable(model, up_reserve[p in P, t in T] >= 0)
        @variable(model, down_reserve[p in P, t in T] >= 0)

        @constraint(model, min_up_reserve[s in S, p in P, t in T], up_reserve[p, t] ≥ up_activation[s, p, t])
        @constraint(model, min_down_reserve[s in S, p in P, t in T], down_reserve[p, t] ≥ down_activation[s, p, t])
        buffer_balancing_capacity = 0
        @constraint(model, area_up_reserve[s in S, t in T, a in A], sum(up_reserve[p, t] for p in P_a[a]) ≥ buffer_balancing_capacity + sum(up_activation[s, p, t] for p in P_a[a]))
        @constraint(model, area_down_reserve[s in S, t in T, a in A], sum(down_reserve[p, t] for p in P_a[a]) ≥ buffer_balancing_capacity + sum(down_activation[s, p, t] for p in P_a[a]))
        
        @variable(model, status[p in P, t in T], Bin)
        @variable(model, startup[p in P, t in T], Bin)

        # Production constraints thermal and wind
        @constraint(model, gen_ub[p in union(P_t, P_h), t in T], production[0, p, t] + up_reserve[p, t] ≤ plant_df[plant_df.plant_id .== p, :gen_ub][1] * status[p, t])
        @constraint(model, gen_lb[p in union(P_t, P_h), t in T], production[0, p, t] - down_reserve[p, t] ≥ plant_df[plant_df.plant_id .== p, :gen_lb][1] * status[p, t])
        # @constraint(model, wind_reserve_ub[p in P_w, t in T], up_reserve[p, t] <= 0)
        # @constraint(model, wind_reserve_lb[p in P_w, t in T], down_reserve[p, t] <= 0)
        @constraint(model, gen_on_off[p in union(P_t, P_h), t in T[2:end]], startup[p, t] ≥ status[p, t] - status[p, t-1])

        @constraint(model, energy_balance[a in A, t in T], 
            sum(production[0, p, t] for p in P_a[a]) 
            + load_shedding[0, a, t] - power_dumping[0, a, t] 
            + sum(transmission[0, l, t] * 0.99 for l in L_in[a])
            - sum(transmission[0, l, t] for l in L_out[a]) 
            == load_df[(load_df.area .== a) .& (load_df.timestep .== t) .& (load_df.scenario .== 0), :load][1])
        @constraint(model, energy_balance_rt[s in S[2:end], a in A, t in T],
            sum(production[s, p, t] for p in P_a[a])
            + load_shedding[s, a, t] - power_dumping[s, a, t]
            + sum(transmission[s, l, t] * 0.99 for l in L_in[a])
            - sum(transmission[s, l, t] for l in L_out[a])
            == load_df[(load_df.area .== a) .& (load_df.timestep .== t) .& (load_df.scenario .== s), :load][1])
        
        @objective(model, Min, 
            π * Δt * sum(production[s, p, t] * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in P_t for t in T for s in S[2:end])
            # + π * Δt * sum(up_activation[s, p, t] * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in P_t for t in T for s in S[2:end])
            # - π * Δt * sum(down_activation[s, p, t] * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in P_t for t in T for s in S[2:end])
            + 0.1 * Δt * sum(up_reserve[p, t] * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in union(P_t, P_h) for t in T)
            + 0.1 * Δt * sum(down_reserve[p, t] * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in union(P_t, P_h) for t in T)
            + π * Δt * sum(load_shedding[s, a, t] * C_shedding for a in A for t in T for s in S[2:end]) 
            + π * Δt * sum(power_dumping[s, a, t] * C_dumping for a in A for t in T for s in S[2:end])
            + π * Δt * sum(C_bypass * flow_bypass[s, p, t] + C_spill * flow_spill[s, p, t] for s in S[2:end] for t in T for p in P_h)
            + π * sum((-volume[s, p, T[end]] + volume[0, p, 0]) * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in P_h for s in S[2:end])
            + sum(startup[p, t] * C_startup for p in P_t for t in T))
    else
        @constraint(model, energy_balance_rt[s in S, a in A, t in T],
            sum(production[s, p, t] for p in P_a[a])
            + load_shedding[s, a, t] - power_dumping[s, a, t]
            + sum(transmission[s, l, t] * 0.99 for l in L_in[a])
            - sum(transmission[s, l, t] for l in L_out[a])
            == load_df[(load_df.area .== a) .& (load_df.timestep .== t) .& (load_df.scenario .== s), :load][1])
    
        @objective(model, Min, 
            π * Δt * sum(production[s, p, t] * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in P_t for t in T for s in S)
            # + π * Δt * sum(up_activation[s, p, t] * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in P_t for t in T for s in S[2:end])
            # - π * Δt * sum(down_activation[s, p, t] * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in P_t for t in T for s in S[2:end])
            # + 0.1 * Δt * sum(up_reserve[p, t] * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in union(P_t, P_h) for t in T)
            # + 0.1 * Δt * sum(down_reserve[p, t] * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in union(P_t, P_h) for t in T)
            + π * Δt * sum(load_shedding[s, a, t] * C_shedding for a in A for t in T for s in S) 
            + π * Δt * sum(power_dumping[s, a, t] * C_dumping for a in A for t in T for s in S)
            + π * Δt * sum(C_bypass * flow_bypass[s, p, t] + C_spill * flow_spill[s, p, t] for s in S for t in T for p in P_h)
            + π * sum((-volume[s, p, T[end]] + volume[0, p, 0]) * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in P_h for s in S)
            + sum(startup[p, t] * C_startup for p in P_t for t in T))
    end
    # print_model_info(model)
    # println(objective_function(model))
    optimize!(model)
    obj_value = objective_value(model)
    println("Objective value: $obj_value")
    return model
end 


function write_results(model)

    wind_df, hydro_df, inflow_df, load_df, plant_df, line_df, A, P, T, L, P_w, P_t, P_h, P_a, L_in, L_out, I_disch, I_spill, I_bypass, C_shedding, C_dumping, C_startup, C_spill, C_bypass, S = read_input_data()

    production = value.(model[:production])
    up_activation = value.(model[:up_activation])
    down_activation = value.(model[:down_activation])
    up_reserve = value.(model[:up_reserve])
    down_reserve = value.(model[:down_reserve])
    
    production_df = DataFrame(id=Int[], scenario=Int[], plant_id=Int[], timestep=Int[], area=Int[],
                                         fuel_type=String[], production=Float64[], up_activation=Float64[], down_activation=Float64[], ub=Float64[], lb=Float64[])
    id_counter = 1
    for p in P
        for t in T
            for s in S
                row = (id_counter, s, p, t, plant_df[plant_df.plant_id .== p, :area][1], 
                        plant_df[plant_df.plant_id .== p, :fuel_type][1], round(production[s, p, t], digits=2),
                        round(up_activation[s, p, t], digits=2), round(down_activation[s, p, t], digits=2), 
                        production[0, p, t] + up_reserve[p, t], production[0, p, t] - down_reserve[p, t])
                push!(production_df, row)
                id_counter += 1
            end
        end
    end
    # XLSX.writetable("discrete_results/results.xlsx", df, overwrite=true, sheetname="production", anchor_cell="A1")


    uc = value.(model[:status])
    startup = value.(model[:startup])
    first_stage_prod_df =  DataFrame(id=Int[], plant_id=Int[], timestep=Int[], up_reserve=Float64[], down_reserve=Float64[], status=Int[], startup=Int[])
    id_counter = 1
    for p in union(P_t, P_h)
        for t in T
            row = (id_counter, p, t, round(up_reserve[p, t], digits=2), round(down_reserve[p, t], digits=2), round(uc[p, t]), round(startup[p, t]))
            push!(first_stage_prod_df, row)
            id_counter += 1
        end
    end

    transmission = value.(model[:transmission])
    transmission_df = DataFrame(id=Int[], scenario=Int[], line_id=Int[], timestep=Int[], area_from=Int[], area_to=Int[], transmission=Float64[])
    id_counter = 1
    for l in L
        for t in T
            for s in S
                row = (id_counter, s, l, t, line_df[line_df.line_id .== l, :area_from][1], line_df[line_df.line_id .==l, :area_to][1], transmission[s, l, t])
                push!(transmission_df, row)
                id_counter += 1
            end
        end
    end
    shedding = value.(model[:load_shedding])
    dumping = value.(model[:power_dumping])
    area_df = DataFrame(id=Int[], scenario=Int[], area=Int[], timestep=Int[], load=Float64[], load_shedding=Float64[], power_dumping=Float64[])
    id_counter = 1
    for a in A
        for t in T
            for s in S
                row = (id_counter, s, a, t, load_df[(load_df.area .== a) .& (load_df.timestep .== t) .& (load_df.scenario .== s), :load][1], shedding[s, a, t], dumping[s, a, t])
                push!(area_df, row)
                id_counter +=1
            end
        end
    end
    # XLSX.writetable("discrete_results/results.xlsx", transmission_df, overwrite=true, sheetname="transmission", anchor_cell="A1")

    discharge = value.(model[:flow_disch])
    spill = value.(model[:flow_spill])
    bypass = value.(model[:flow_bypass])
    inflow = value.(model[:total_flow_in])
    outflow = value.(model[:total_flow_out])
    volume_res = value.(model[:volume])

    hydro_results_df = DataFrame(id=Int[], scenario = Int[], plant_id = Int[], timestep=Int[], discharge=Float64[], spill=Float64[], bypass=Float64[], controlled_inflow=Float64[], controlled_outflow=Float64[], volume_res=Float64[])
    id_counter = 0
    for p in P_h
        for t in T
            for s in S
                row = (id_counter, s, p, t, discharge[s, p, t], spill[s, p, t], bypass[s, p, t], inflow[s, p, t], outflow[s, p, t], volume_res[s, p, t])
                push!(hydro_results_df, row)
                id_counter += 1
            end
        end
    end
    CSV.write("discrete_results/transmission.csv", transmission_df)
    CSV.write("discrete_results/production.csv", production_df)
    CSV.write("discrete_results/area_results.csv", area_df)
    CSV.write("discrete_results/hydro_results.csv", hydro_results_df)

    XLSX.writetable("discrete_results/results.xlsx", "transmission" => transmission_df, "first_stage" => first_stage_prod_df, "production" => production_df, "area_results" => area_df, "hydro_results" => hydro_results_df, overwrite=true)

end


function write_expanded_results(s)
    sheetnames = []
    filename = "discrete_results/results.xlsx"
    XLSX.openxlsx(filename) do xf
        sheetnames = XLSX.sheetnames(xf)
    end
    output = Dict{String, DataFrame}()
    fractional_timesteps = 0:(1/s):24
    for sheet in sheetnames
        println("Writing expanded discrete results: $sheet")
        data = DataFrame(XLSX.readtable(filename, sheet, infer_eltypes=true))
        if sheet == "transmission"
            sort!(data, [:line_id, :scenario, :timestep])
            group_symbols = [:line_id, :scenario]
        elseif sheet == "first_stage"
            sort!(data, [:plant_id, :timestep])
            group_symbols = [:plant_id]
        elseif sheet == "production" || sheet == "hydro_results"
            sort!(data, [:plant_id, :scenario, :timestep])
            group_symbols = [:plant_id, :scenario]
        elseif sheet == "area_results"
            sort!(data, [:area, :scenario, :timestep])
            group_symbols = [:area, :scenario]
        end

        
        T = unique(data.timestep)
        timesteps_per_hour = div(T[end], 24)
        repeats_per_datapoint = div(s, timesteps_per_hour)
        expanded_data = DataFrame()
        for col in names(data)
            expanded_data[!, col] = repeat(data[!, col], inner=repeats_per_datapoint)
        end
        grouped = groupby(expanded_data, group_symbols)
        results = DataFrame()
        for g in grouped
            first_row = g[1, :] |> DataFrame
            g = vcat(first_row, g)
            if length(fractional_timesteps) != nrow(g)
                error("Mismatch between expected and actual timestep lengths for group.")
            end
            g.timestep_fractional = fractional_timesteps
            append!(results, g)
        end
        # expanded_data.timestep_fractional = [div(timestep - 1, timesteps_per_hour) + (1 / s) * (mod(i-1, s)) for (i, timestep) in enumerate(expanded_data.timestep)]
        # display(expanded_data)
        results.timestep_fractional .= round.(results.timestep_fractional, digits=4)
        output[sheet] = results
    end
    for (sheetname, df) in output
        CSV.write("discrete_results/$sheetname" * "_expanded.csv", df)
    end
    
    # XLSX.writetable("discrete_results/results_expanded.xlsx",   "transmission" => output["transmission"], 
    #                                                             "first_stage" => output["first_stage"], 
    #                                                             "production" => output["production"], 
    #                                                             "area_results" => output["area_results"], 
    #                                                             "hydro_results" => output["hydro_results"], overwrite=true)

    # XLSX.writetable("discrete_results/results_expanded.xlsx", modified_sheets, overwrite=true)
    # XLSX.openxlsx("discrete_results/results_expanded.xlsx", mode="w") do new_xf
    #     for (sheetname, df) in modified_sheets
    #         println("Writing sheet: $sheetname")
    #         XLSX.writedata!(new_xf, sheetname, Tables.matrix(df))
end

