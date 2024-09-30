using DataFrames
using XLSX
using JuMP
using CPLEX


include("C:/Users/vegardvk/vscodeProjects/bernstein/get_hydro_data.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/helper_functions.jl")

function read_input_data()
    # Read input data
    wind_df = DataFrame(XLSX.readtable("output/wind_ts_data.xlsx", "Sheet1", infer_eltypes=true))
    hydro_df = DataFrame(XLSX.readtable("output/hydro_data.xlsx", "Sheet1", infer_eltypes=true))
    inflow_df = DataFrame(XLSX.readtable("output/inflow_data.xlsx", "Sheet1", infer_eltypes=true))
    load_df = DataFrame(XLSX.readtable("output/load_data.xlsx", "Sheet1", infer_eltypes=true))
    plant_df = DataFrame(XLSX.readtable("output/plant_data.xlsx", "Sheet1", infer_eltypes=true))
    line_df = DataFrame(XLSX.readtable("output/line_data.xlsx", "Sheet1", infer_eltypes=true))

    A = unique(plant_df.area)
    P = unique(plant_df.plant_id)
    T = unique(load_df.timestep)
    L = unique(line_df.line_id)

    P_w = unique(plant_df[plant_df.fuel_type .== "Wind", :plant_id])
    P_t = unique(plant_df[plant_df.fuel_type .== "Thermal", :plant_id])
    P_h = unique(plant_df[plant_df.fuel_type .== "Hydro", :plant_id])

    area_grouped_plants = groupby(plant_df, :area)
    P_a = Dict(g.area[1] => unique(g.plant_id) for g in area_grouped_plants)
    
    # Dicts with the lines that go in/out of area a
    L_in = Dict(a => unique(line_df[line_df.area_to .== a, :line_id]) for a in A) 
    L_out =  Dict(a => unique(line_df[line_df.area_from .== a, :line_id]) for a in A)

    I_disch, I_spill, I_bypass = find_connected_plants(hydro_df)

    C_shedding, C_dumping, C_startup = get_cost_parameters()
    return wind_df, hydro_df, inflow_df, load_df, plant_df, line_df, A, P, T, L, P_w, P_t, P_h, P_a, L_in, L_out, I_disch, I_spill, I_bypass, C_shedding, C_dumping, C_startup
end


function define_and_solve_model()
    wind_df, hydro_df, inflow_df, load_df, plant_df, line_df, A, P, T, L, P_w, P_t, P_h, P_a, L_in, L_out, I_disch, I_spill, I_bypass, C_shedding, C_dumping, C_startup = read_input_data()

    model = Model()
    set_optimizer(model, CPLEX.Optimizer)

    @variable(model, production[p in P, t in T] >= 0)
    @variable(model, load_shedding[a in A, t in T] >= 0)
    @variable(model, power_dumping[a in A, t in T] >= 0)
    @variable(model, transmission[l in L, t in T] >= 0)

    # @variable(model, 0 ≤ x1[p in P, t in T] ≤ 1)
    # @variable(model, 0 ≤ x2[p in P, t in T] ≤ 1)
    # @variable(model, startup[p in P, t in T] ≥ 0)
    # @constraint(model, production_linking[p in P, t in T], production[p, t] == x1[p, t] * p_dict["gen_lb"][p] + x2[p, t] * (p_dict["gen_ub"][p] - p_dict["gen_lb"][p]))
    # @constraint(model, start_before_prod[p in P, t in T], x1[p, t] ≥ x2[p, t])
    # @constraint(model, startup_count[p in P, t in T[2:end]], x1[p, t] - x1[p, t-1]  == startup[p, t])

    @variable(model, status[p in P, t in T], Bin)
    @variable(model, startup[p in P, t in T], Bin)
    @constraint(model, gen_ub[p in P, t in T], production[p, t] ≤ plant_df[plant_df.plant_id .== p, :gen_ub][1] * status[p, t])
    @constraint(model, wind_prod[p in P_w, t in T], production[p, t] == wind_df[(wind_df.plant_id .== p) .& (wind_df.timestep .== t), :wind_power][1])
    @constraint(model, gen_lb[p in P, t in T], production[p, t] ≥ plant_df[plant_df.plant_id .== p, :gen_lb][1] * status[p, t])
    @constraint(model, gen_on_off[p in P, t in T[2:end]], startup[p, t] ≥ status[p, t] - status[p, t-1])

    @variable(model, flow_disch[p in P_h, t in T] ≥ 0) # Antar at alle moduler bare har ett discharge-segment
    @variable(model, flow_bypass[p in P_h, t in T] ≥ 0)
    @variable(model, flow_spill[p in P_h, t in T] ≥ 0)
    @variable(model, total_flow_in[p in P_h, t in T] ≥ 0)
    @variable(model, total_flow_out[p in P_h, t in T] ≥ 0)
    @variable(model, volume[p in P_h, t in 0:T[end]] ≥ 0)

    @constraint(model, controlled_inflow[p in P_h, t in T], total_flow_in[p, t] == sum(flow_disch[i, t] for i in I_disch[p]) 
                                                                                + sum(flow_bypass[i, t] for i in I_bypass[p]) 
                                                                                + sum(flow_spill[i, t] for i in I_spill[p])
                                                                                + inflow_df[(inflow_df.plant_id .== p) .& (inflow_df.timestep .== t), :inflow][1])

    @constraint(model, controlled_outflow[p in P_h, t in T], total_flow_out[p, t] == flow_disch[p, t] + flow_bypass[p, t] + flow_spill[p, t])
    @constraint(model, starting_reservoir_constraint[p in P_h], volume[p, 0] == hydro_df[hydro_df.plant_id .== p, :starting_reservoir][1])

    @constraint(model, reservoir_balance[p in P_h, t in 1:(T[end])], volume[p, t] - volume[p, t-1] == total_flow_in[p, t] - total_flow_out[p, t])

    @constraint(model, hydro_production[p in P_h, t in T], production[p, t] == hydro_df[hydro_df.plant_id .== p, :enekv][1] * flow_disch[p, t] * 3.6) # kwh/M3 * M3/s * 3600 s/h * 1/1000 MW/kW
    @constraint(model, vol_ub[p in P_h, t in T], volume[p, t] ≤ hydro_df[hydro_df.plant_id .== p, :kap_mag][1])
    @constraint(model, bypass_ub[p in P_h, t in T], flow_bypass[p, t] ≤ hydro_df[hydro_df.plant_id .== p, :kap_forb][1])
    @constraint(model, prod_ub[p in P_h, t in T], flow_disch[p, t] ≤ hydro_df[hydro_df.plant_id .== p, :kap_gen_m3s][1])
    @constraint(model, spill_ub[p in P_h, t in T], flow_spill[p, t] ≤ hydro_df[hydro_df.plant_id .== p, :kap_spill][1])
    @constraint(model, transmission_ub[l in L, t in T], transmission[l, t] ≤ line_df[line_df.line_id .== l, :capacity][1])
    @constraint(model, transmission_lb[l in L, t in T], transmission[l, t] ≥ -line_df[line_df.line_id .== l, :capacity][1])

    @constraint(model, energy_balance[a in A, t in T], sum(production[p, t] for p in P_a[a]) 
                + load_shedding[a, t] - power_dumping[a, t] 
                + sum(transmission[l, t] * 0.99 for l in L_in[a])
                - sum(transmission[l, t] for l in L_out[a]) 
                == load_df[(load_df.area .== a) .& (load_df.timestep .== t), :Forbruk][1])

    @objective(model, Min, sum(production[p, t] * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in P_t for t in T)
                    + sum((-volume[p, T[end]] + volume[p, 0]) * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in P_h) 
                    + sum(load_shedding[a, t] * C_shedding for a in A for t in T) 
                    + sum(power_dumping[a, t] * C_dumping for a in A for t in T)
                    + sum(startup[p, t] * C_startup for p in P_t for t in T))

    # print_model_info(model)
    # println(objective_function(model))
    optimize!(model)
    obj_value = objective_value(model)
    println("Objective value: $obj_value")
    return model
end


function write_results(model)

    wind_df, hydro_df, inflow_df, load_df, plant_df, line_df, A, P, T, L, P_w, P_t, P_h, P_a, L_in, L_out, I_disch, I_spill, I_bypass, C_shedding, C_dumping, C_startup = read_input_data()

    production = value.(model[:production])
    uc = value.(model[:status])
    startup = value.(model[:startup])
    production_df = DataFrame(id=Int[], plant_id=Int[], timestep=Int[], area=Int[], fuel_type=String[], production=Float64[], status=Int[], startup=Int[])
    id_counter = 1
    for p in P
        for t in T
            row = (id_counter, p, t, plant_df[plant_df.plant_id .== p, :area][1], plant_df[plant_df.plant_id .== p, :fuel_type][1], round(production[p, t], digits=2), round(uc[p, t]), round(startup[p, t]))
            push!(production_df, row)
            id_counter += 1
        end
    end
    # XLSX.writetable("discrete_results/results.xlsx", df, overwrite=true, sheetname="production", anchor_cell="A1")

    transmission = value.(model[:transmission])
    transmission_df = DataFrame(id=Int[], line_id=Int[], timestep=Int[], area_from=Int[], area_to=Int[], transmission=Float64[])
    id_counter = 1
    for l in L
        for t in T
            row = (id_counter, l, t, line_df[line_df.line_id .== l, :area_from][1], line_df[line_df.line_id .==l, :area_to][1], transmission[l, t])
            push!(transmission_df, row)
            id_counter += 1
        end
    end

    shedding = value.(model[:load_shedding])
    dumping = value.(model[:power_dumping])
    area_df = DataFrame(id=Int[], area=Int[], timestep=Int[], load=Float64[], load_shedding=Float64[], power_dumping=Float64[])
    id_counter = 1
    for a in A
        for t in T
            row = (id_counter, a, t, load_df[(load_df.area .== a) .& (load_df.timestep .== t), :Forbruk][1], shedding[a, t], dumping[a, t])
            push!(area_df, row)
            id_counter +=1
        end
    end
    # XLSX.writetable("discrete_results/results.xlsx", transmission_df, overwrite=true, sheetname="transmission", anchor_cell="A1")

    discharge = value.(model[:flow_disch])
    spill = value.(model[:flow_spill])
    bypass = value.(model[:flow_bypass])
    inflow = value.(model[:total_flow_in])
    outflow = value.(model[:total_flow_out])
    volume_res = value.(model[:volume])

    hydro_results_df = DataFrame(id=Int[], plant_id = Int[], timestep=Int[], discharge=Float64[], spill=Float64[], bypass=Float64[], controlled_inflow=Float64[], controlled_outflow=Float64[], volume_res=Float64[])
    id_counter = 0
    for p in P_h
        for t in T
            row = (id_counter, p, t, discharge[p, t], spill[p, t], bypass[p, t], inflow[p, t], outflow[p, t], volume_res[p, t])
            push!(hydro_results_df, row)
            id_counter += 1
        end
    end
    CSV.write("discrete_results/transmission.csv", transmission_df)
    CSV.write("discrete_results/production.csv", production_df)
    CSV.write("discrete_results/area_results.csv", area_df)
    CSV.write("discrete_results/hydro_results.csv", hydro_results_df)

    XLSX.writetable("discrete_results/results.xlsx", "transmission" => transmission_df, "production" => production_df, "area_results" => area_df, "hydro_results" => hydro_results_df, overwrite=true)

end



