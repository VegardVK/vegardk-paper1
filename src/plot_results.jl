using DataFrames
using XLSX
using CSV
using UnPack

include(joinpath(@__DIR__, "helper_functions.jl"))
include(joinpath(@__DIR__, "find_bernstein_weights.jl"))


function import_data_discrete()
    model_results =(
        prod_df = DataFrame(XLSX.readtable("discrete_results/results.xlsx", "production", infer_eltypes=true)),
        first_stage_df = DataFrame(XLSX.readtable("discrete_results/results.xlsx", "first_stage", infer_eltypes=true)),
        flow_df = DataFrame(XLSX.readtable("discrete_results/results.xlsx", "transmission", infer_eltypes=true)),
        area_df = DataFrame(XLSX.readtable("discrete_results/results.xlsx", "area_results", infer_eltypes=true)),
        hydro_df = DataFrame(XLSX.readtable("discrete_results/results.xlsx", "hydro_results", infer_eltypes=true)),
        line_df = DataFrame(XLSX.readtable("output/line_data.xlsx", "Sheet1", infer_eltypes=true)),
    )
    @unpack prod_df, flow_df, line_df = model_results

    A = unique(prod_df.area)
    area_grouped_plants = groupby(prod_df, :area)
    input_sets = (
        A = unique(prod_df.area),
        P = unique(prod_df.plant_id),
        T = unique(prod_df.timestep),
        L = unique(flow_df.line_id),
        S = unique(prod_df.scenario),
        L_in = Dict(a => unique(line_df[line_df.area_to .== a, :line_id]) for a in A) ,
        L_out =  Dict(a => unique(line_df[line_df.area_from .== a, :line_id]) for a in A),
        P_a = Dict(g.area[1] => unique(g.plant_id) for g in area_grouped_plants)
    )
    return model_results, input_sets
    # return prod_df, flow_df, area_df, hydro_df, line_df, first_stage_df, A, P, T, L, S, L_in, L_out, P_a
end

function import_data_continuous()
    prod_weights_df = DataFrame(XLSX.readtable("continuous_results/result_weights.xlsx", "production", infer_eltypes=true))
    first_stage_prod = DataFrame(XLSX.readtable("output/production_weights.xlsx", "Sheet1", infer_eltypes=true))
    prod_aggregated_df = DataFrame(XLSX.readtable("continuous_results/results_aggregated.xlsx", "production", infer_eltypes=true))
    area_df = DataFrame(XLSX.readtable("continuous_results/result_weights.xlsx", "area_results", infer_eltypes=true))
    hydro_df = DataFrame(XLSX.readtable("continuous_results/results_aggregated.xlsx", "hydro_results", infer_eltypes=true))
    hydro_weights_df = DataFrame(XLSX.readtable("continuous_results/result_weights.xlsx", "hydro_results", infer_eltypes=true))
    discrete_results_df = CSV.read("discrete_results/discrete_objective.csv", DataFrame)

    A = unique(prod_weights_df.area)
    P = unique(prod_weights_df.plant_id)
    T = unique(prod_weights_df.timestep)
    B = unique(prod_weights_df.b)
    S = unique(prod_weights_df.scenario)

    return prod_weights_df, first_stage_prod, prod_aggregated_df, area_df, hydro_df, hydro_weights_df, discrete_results_df, A, P, T, B, S
end


function plot_energy_balance()
    model_results, input_sets  = import_data_discrete()
    @unpack prod_df, flow_df, area_df, hydro_df, line_df = model_results
    @unpack A, P, T, L, L_in, L_out, P_a = input_sets

    prod_df_wide = unstack(select(prod_df, :timestep, :production, :plant_id), :plant_id, :production)
    area_load_df_wide = unstack(select(area_df, :timestep, :load, :area), :area, :load)
    export_df = unstack(select(flow_df, :area_from, :timestep, :transmission), :area_from, :transmission, combine=sum)
    import_df = unstack(select(flow_df, :area_to, :timestep, :transmission), :area_to, :transmission, combine=sum)
    dumping_df = unstack(select(area_df, :area, :timestep, :power_dumping), :area, :power_dumping)
    shedding_df = unstack(select(area_df, :area, :timestep, :load_shedding), :area, :load_shedding)

    area_balance_df = DataFrame()
    for a in A
        area_plants = P_a[a]
        area_plants_production_df = select(prod_df_wide, ["$a" for a in area_plants])
        # area_balance_df = hcat(area_balance_df, area_plants_production_df)
        area_balance_df.production = sum(eachcol(area_plants_production_df))

        area_balance_df.load = area_load_df_wide[!, "$a"]
        trade_df = DataFrame()
        if "$a" in names(export_df)
            trade_df.exported_power = export_df[!, "$a"]
        else
            trade_df.exported_power = zeros(nrow(area_balance_df))
        end
        if "$a" in names(import_df)
            trade_df.imported_power = import_df[!, "$a"]
        else
            trade_df.imported_power = zeros(nrow(area_balance_df))
        end
        area_balance_df.net_position = trade_df[!, :exported_power] .- trade_df[!, :imported_power]
        area_balance_df.load_shedding = shedding_df[!, "$a"]
        area_balance_df.power_dumping = dumping_df[!, "$a"]

        area_balance_df.timestep = unique(prod_df.timestep)
        # plot_all_columns_df(area_plants_production_df, "Power plants in area $a")
        plot_all_columns_df(area_balance_df, "Energy balance area $a")
    end
end

function print_simple_results()
    prod_df, flow_df, area_df, hydro_df, line_df, A, P, T, L, L_in, L_out, P_a = import_data_discrete()
    for a in A
        area_load_shedding = sum([area_df[area_df.area .== a, :load_shedding]])

    end
end

function plot_hydro_balance()
    prod_df, flow_df, area_df, hydro_df, line_df, A, P, T, L, L_in, L_out, P_a = import_data_discrete()
    P_h = unique(hydro_df.plant_id)
    for p in P_h
        # display(hydro_df[hydro_df.plant_id .== p, :])
        pp_df = hydro_df[hydro_df.plant_id .== p, :]
        pp_df = select(pp_df, :timestep, :controlled_inflow, :controlled_outflow, :volume_res)
        pp_df.controlled_outflow = -pp_df.controlled_outflow
        plot_all_columns_df(pp_df, "Hydro balance power plant $p")
        # display(pp_df)
    end
end

function calculate_objective_components_discrete()
    model_results, input_sets  = import_data_discrete()
    @unpack prod_df, flow_df, area_df, hydro_df, line_df, first_stage_df = model_results
    @unpack A, P, T, L, S, L_in, L_out, P_a = input_sets

    Δt = 24/T[end]
    π = 1/S[end]

    plant_df = DataFrame(XLSX.readtable("output/plant_data.xlsx", "Sheet1", infer_eltypes=true))
    hydro_data_df = DataFrame(XLSX.readtable("output/hydro_data.xlsx", "Sheet1", infer_eltypes=true))

    area_grouped_plants = groupby(plant_df, :area)
    P_a = Dict(g.area[1] => unique(g.plant_id) for g in area_grouped_plants)

    C_shedding, C_dumping, C_startup, C_spill, C_bypass = get_cost_parameters()
    P_t = unique(plant_df[plant_df.fuel_type .== "Thermal", :plant_id])
    P_h = unique(plant_df[plant_df.fuel_type .== "Hydro", :plant_id])

    # prod_df_wide = unstack(select(prod_df, :timestep, :production, :plant_id), :plant_id, :production)

    day_ahead_fuel_costs = round(Δt * sum(first_stage_df[(first_stage_df.timestep .== t) .& (first_stage_df.plant_id .== p), :first_stage_prod][1] 
            * plant_df[plant_df.plant_id .== p, :fuel_price][1]  for p in P_t for t in T), digits=2)

    up_activation_costs = round(π * Δt * sum(prod_df[(prod_df.scenario .== s) .& (prod_df.timestep .== t).& (prod_df.plant_id .== p), :up_activation][1] 
            * plant_df[plant_df.plant_id .== p, :fuel_price][1]  for p in P_t for t in T for s in S), digits=2)

    down_activation_costs = - round(π * Δt * sum(prod_df[(prod_df.scenario .== s) .& (prod_df.timestep .== t).& (prod_df.plant_id .== p), :down_activation][1] 
            * plant_df[plant_df.plant_id .== p, :fuel_price][1]  for p in P_t for t in T for s in S), digits=2)

    alternative_fuel_cost = day_ahead_fuel_costs + up_activation_costs + down_activation_costs

    total_fuel_costs = round(Δt * π * sum(prod_df[(prod_df.timestep .== t) .& (prod_df.plant_id .== p) .& (prod_df.scenario .== s), :production][1] * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in P_t for t in T for s in S), digits=2)

    startup_costs = round(sum(first_stage_df[(first_stage_df.timestep .== t) .& (first_stage_df.plant_id .== p), :startup][1] * C_startup for p in P_t for t in T), digits=2)
    # startup_costs = sum(prod_df[(prod_df.timestep .== t) .& (prod_df.plant_id .== p), :startup][1] * C_startup for p in P_t for t in T)
    
    dumping_costs = round(Δt * π * sum(area_df[(area_df.area .== a) .& (area_df.timestep .== t) .& (area_df.scenario .== s), :power_dumping][1] * C_dumping for a in A for t in T for s in S), digits=2)
    shed_costs = round(Δt * π * sum(area_df[(area_df.area .== a) .& (area_df.timestep .== t) .& (area_df.scenario .== s), :load_shedding][1] * C_shedding for a in A for t in T for s in S), digits=2)
    volume_costs = round(π * sum((-hydro_df[(hydro_df.plant_id .== p) .& (hydro_df.timestep .== T[end]) .& (hydro_df.scenario .== s), :volume_res][1]
                                       + hydro_data_df[hydro_data_df.plant_id .== p, :starting_reservoir][1]) 
                                       * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in P_h for s in S), digits=2)

    bypass_costs = round(π * Δt * sum(C_bypass * hydro_df[(hydro_df.plant_id .== p) .& (hydro_df.scenario .== s) .& (hydro_df.timestep .== t), :bypass][1] for t in T for p in P_h for s in S), digits=2)
    spill_costs = round(π * Δt * sum(C_spill * hydro_df[(hydro_df.plant_id .== p) .& (hydro_df.scenario .== s) .& (hydro_df.timestep .== t), :spill][1] for t in T for p in P_h for s in S), digits=2)

    up_reserve_costs = round(0.1 * Δt * sum(first_stage_df[(first_stage_df.timestep .== t) .& (first_stage_df.plant_id .== p), :up_reserve][1] 
                                * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in union(P_t, P_h) for t in T), digits=2)
    down_reserve_costs = round(0.1 * Δt * sum(first_stage_df[(first_stage_df.timestep .== t) .& (first_stage_df.plant_id .== p), :down_reserve][1] 
                                * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in union(P_t, P_h) for t in T), digits=2)                                                     
    objective = total_fuel_costs + startup_costs + dumping_costs + shed_costs + volume_costs + up_reserve_costs + down_reserve_costs + bypass_costs + spill_costs

    println("Objective: $objective")
    println("\t Production costs: $total_fuel_costs")
    println("\t Water costs: $volume_costs")
    println("\t Startup costs: $startup_costs")
    println("\t shed_costs: $shed_costs")
    println("\t dump costs: $dumping_costs")
    println("\t Up-reserve costs: $up_reserve_costs")
    println("\t Down-reserve costs: $down_reserve_costs")
    println("\t Spill costs: $spill_costs")
    println("\t Bypass costs: $bypass_costs")

    println("\nAlternative fuel cost calculation:")
    println("\t Day-ahead fuel costs: $day_ahead_fuel_costs")
    println("\t Up-activation costs: $up_activation_costs")
    println("\t Down-activation costs: $down_activation_costs")
    println("\t Total: $alternative_fuel_cost")

    df = DataFrame()
    for s in S
        for a in A
            shed_costs = Δt * π * sum(area_df[(area_df.area .== a) .& (area_df.timestep .== t) .& (area_df.scenario .== s), :load_shedding][1] * C_shedding for t in T)
            dump_costs = Δt * π * sum(area_df[(area_df.area .== a) .& (area_df.timestep .== t) .& (area_df.scenario .== s), :power_dumping][1] * C_dumping for t in T)

            fuel_costs = 0
            startup_costs = 0
            water_costs = 0
            bypass_costs = 0
            spill_costs = 0
            up_reserve_costs = 0
            down_reserve_costs = 0
            day_ahead_costs = 0
            up_activation_costs = 0
            down_activation_costs = 0
            for p in P_a[a]

                if p in P_t
                    # fuel_costs += Δt * π * sum(prod_df[(prod_df.timestep .== t) .& (prod_df.plant_id .== p) .& (prod_df.scenario .== s), :production][1] * plant_df[plant_df.plant_id .== p, :fuel_price][1] for t in T)
                    startup_costs += sum(first_stage_df[(first_stage_df.timestep .== t) .& (first_stage_df.plant_id .== p), :startup][1] * C_startup for t in T)
                    fuel_price = plant_df[plant_df.plant_id .== p, :fuel_price][1]
                elseif p in P_h
                    water_costs += π * (-hydro_df[(hydro_df.plant_id .== p) .& (hydro_df.timestep .== T[end]) .& (hydro_df.scenario .== s), :volume_res][1] 
                        + hydro_data_df[hydro_data_df.plant_id .== p, :starting_reservoir][1]) * plant_df[plant_df.plant_id .== p, :fuel_price][1]

                        bypass_costs += round(π * Δt * sum(C_bypass * hydro_df[(hydro_df.plant_id .== p) .& (hydro_df.scenario .== s) .& (hydro_df.timestep .== t), :bypass][1] for t in T), digits=2)
                        spill_costs += round(π * Δt * sum(C_spill * hydro_df[(hydro_df.plant_id .== p) .& (hydro_df.scenario .== s) .& (hydro_df.timestep .== t), :spill][1] for t in T), digits=2)
                    fuel_price = get_water_value()    
                end
                if p in union(P_t, P_h)
                    up_reserve_costs += Δt * 0.1 * sum(first_stage_df[(first_stage_df.timestep .== t) .& (first_stage_df.plant_id .== p), :up_reserve][1] * plant_df[plant_df.plant_id .== p, :fuel_price][1] for t in T)
                    down_reserve_costs += Δt * 0.1 * sum(first_stage_df[(first_stage_df.timestep .== t) .& (first_stage_df.plant_id .== p), :down_reserve][1] * plant_df[plant_df.plant_id .== p, :fuel_price][1] for t in T)
                    day_ahead_costs += round(π * Δt *  sum(first_stage_df[(first_stage_df.timestep .== t) .& (first_stage_df.plant_id .== p), :first_stage_prod][1] * fuel_price for t in T), digits=2)
                    up_activation_costs += round(π * Δt * sum(prod_df[(prod_df.scenario .== s) .& (prod_df.timestep .== t).& (prod_df.plant_id .== p), :up_activation][1] * fuel_price for t in T), digits=2)
                    down_activation_costs += - round(π * Δt * sum(prod_df[(prod_df.scenario .== s) .& (prod_df.timestep .== t).& (prod_df.plant_id .== p), :down_activation][1] *fuel_price for t in T), digits=2)
                end
            end

            area_cost_df = DataFrame(
                component = ["Day-ahead costs", "up-activation cost", "down-activation cost", "Shedding", "Dumping", "Water", "Start", "Up-reserve", "Down-reserve", "Bypass cost", "Spill cost"],
                values = [day_ahead_costs, up_activation_costs, down_activation_costs, shed_costs, dump_costs, water_costs, startup_costs, up_reserve_costs, down_reserve_costs, bypass_costs, spill_costs],
                area = a, scenario = s)
            append!(df, area_cost_df)
        end
    end
    CSV.write("discrete_results/discrete_objective.csv", df)
end

function calculate_objective_components_continuous()
    plant_df = DataFrame(XLSX.readtable("output/plant_data.xlsx", "Sheet1", infer_eltypes=true))
    hydro_data_df = DataFrame(XLSX.readtable("output/hydro_data.xlsx", "Sheet1", infer_eltypes=true))
    prod_weights_df, first_stage_prod, prod_aggregated_df, area_df, hydro_df, hydro_weights_df, discrete_results_df, A, P, T, B, S = import_data_continuous()
    Δt = 24/T[end]
    π = 1/length(S)

    area_grouped_plants = groupby(plant_df, :area)
    P_a = Dict(g.area[1] => unique(g.plant_id) for g in area_grouped_plants)

    C_shedding, C_dumping, C_startup, C_spill, C_bypass = get_cost_parameters()
    P_t = unique(plant_df[plant_df.fuel_type .== "Thermal", :plant_id])
    P_h = unique(plant_df[plant_df.fuel_type .== "Hydro", :plant_id])

    day_ahead_fuel_costs = round(1/(B[end]+1) * Δt * sum(first_stage_prod[(first_stage_prod.scenario .== 0) .& (first_stage_prod.timestep .== t) .& (first_stage_prod.plant_id .== p) .& (first_stage_prod.b .==b), :production][1] 
                                                        * plant_df[plant_df.plant_id .== p, :fuel_price][1]  for p in P_t for b in B for t in T), digits=2)

    up_activation_costs = round(1/(B[end]+1) * π * Δt * sum(prod_weights_df[(prod_weights_df.scenario .== s) .& (prod_weights_df.timestep .== t).& (prod_weights_df.plant_id .== p) .& (prod_weights_df.b .==b), :up_activation][1] 
                                                            * plant_df[plant_df.plant_id .== p, :fuel_price][1]  for p in P_t for b in B for t in T for s in S), digits=2)

    down_activation_costs = - round(1/(B[end]+1) * π * Δt * sum(prod_weights_df[(prod_weights_df.scenario .== s) .& (prod_weights_df.timestep .== t).& (prod_weights_df.plant_id .== p) .& (prod_weights_df.b .==b), :down_activation][1] 
                                                            * plant_df[plant_df.plant_id .== p, :fuel_price][1]  for p in P_t for b in B for t in T for s in S), digits=2)

    alternative_fuel_cost = day_ahead_fuel_costs + up_activation_costs + down_activation_costs
    total_fuel_costs = round(1/(B[end]+1) * π * Δt * sum(prod_weights_df[(prod_weights_df.scenario .== s) .& (prod_weights_df.timestep .== t) .& (prod_weights_df.plant_id .== p) .& (prod_weights_df.b .==b), :production][1] * plant_df[plant_df.plant_id .== p, :fuel_price][1]  for p in P_t for b in B for t in T for s in S), digits=2)
    shed_costs =  round(1/(B[end]+1) * π * Δt * sum(area_df[(area_df.scenario .== s) .& (area_df.timestep .== t) .& (area_df.area .== a) .& (area_df.b .== b), :shedding][1] * C_shedding for b in B for a in A for t in T for s in S), digits=2)
    dump_costs = round(1/(B[end]+1) * π * Δt * sum(area_df[(area_df.scenario .== s) .& (area_df.timestep .== t) .& (area_df.area .== a) .& (area_df.b .== b), :dumping][1] * C_dumping for b in B for a in A for t in T for s in S), digits=2)
    water_costs = round(π * sum((-hydro_df[(hydro_df.scenario .== s) .& (hydro_df.plant_id .== p) .& (hydro_df.timestep .== T[end]), :volume_res][1] + hydro_data_df[hydro_data_df.plant_id .== p, :starting_reservoir][1]) * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in P_h for s in S), digits=2)
    bypass_costs = round(1/(B[end]+1) * π * Δt * sum(C_bypass * hydro_weights_df[(hydro_weights_df.plant_id .== p) .& (hydro_weights_df.scenario .== s) .& (hydro_weights_df.timestep .== t) .& (hydro_weights_df.b .== b), :bypass][1] for t in T for b in B for p in P_h for s in S), digits=2)
    spill_costs = round(1/(B[end]+1) * π * Δt * sum(C_spill * hydro_weights_df[(hydro_weights_df.plant_id .== p) .& (hydro_weights_df.scenario .== s) .& (hydro_weights_df.timestep .== t) .& (hydro_weights_df.b .== b), :spill][1] for t in T for b in B for p in P_h for s in S), digits=2)
    #startup_costs = sum(prod_aggregated_df[(prod_aggregated_df.timestep .== t) .& (prod_aggregated_df.plant_id .== p), :startup][1] * C_startup for p in P_t for t in T)
    objective = total_fuel_costs + shed_costs + dump_costs + water_costs + bypass_costs + spill_costs # + startup_costs

    println("Objective: $objective")
    println("\t Production costs: $total_fuel_costs")
    println("\t Water costs: $water_costs")
    println("\t shed_costs: $shed_costs")
    println("\t dump costs: $dump_costs")
    println("\t Spill costs: $spill_costs")
    println("\t Bypass costs: $bypass_costs")
    println("\nAlternative fuel cost calculation:")
    println("\t Day-ahead fuel costs: $day_ahead_fuel_costs")
    println("\t Up-activation costs: $up_activation_costs")
    println("\t Down-activation costs: $down_activation_costs")
    println("\t Total: $alternative_fuel_cost")
    # println("\t Startup costs: $startup_costs")

    df = DataFrame()
    for a in A
        for s in S
            shed_costs = round(1/(B[end]+1) * π * Δt * sum(area_df[(area_df.scenario .== s) .& (area_df.timestep .== t) .& (area_df.area .== a) .& (area_df.b .== b), :shedding][1] * C_shedding for b in B for t in T), digits=2)
            dump_costs = round(1/(B[end]+1) * π * Δt * sum(area_df[(area_df.scenario .== s) .& (area_df.timestep .== t) .& (area_df.area .== a) .& (area_df.b .== b), :dumping][1] * C_dumping for b in B for t in T), digits=2)
            
            fuel_costs = 0
            day_ahead_costs = 0
            up_activation_costs = 0
            down_activation_costs = 0

            startup_costs = 0
            water_costs = 0
            spill_costs = 0
            bypass_costs = 0
            for p in P_a[a]
                if p in P_t
                    day_ahead_costs += round(1/(B[end]+1)* π * Δt * sum(first_stage_prod[(first_stage_prod.scenario .== 0) .& (first_stage_prod.timestep .== t) .& (first_stage_prod.plant_id .== p) .& (first_stage_prod.b .==b), :production][1] 
                                                * plant_df[plant_df.plant_id .== p, :fuel_price][1] for b in B for t in T), digits=2)

                    up_activation_costs += round(1/(B[end]+1) * π * Δt * sum(prod_weights_df[(prod_weights_df.scenario .== s) .& (prod_weights_df.timestep .== t).& (prod_weights_df.plant_id .== p) .& (prod_weights_df.b .==b), :up_activation][1] 
                                                * plant_df[plant_df.plant_id .== p, :fuel_price][1] for b in B for t in T), digits=2)

                    down_activation_costs += - round(1/(B[end]+1) * π * Δt * sum(prod_weights_df[(prod_weights_df.scenario .== s) .& (prod_weights_df.timestep .== t).& (prod_weights_df.plant_id .== p) .& (prod_weights_df.b .==b), :down_activation][1] 
                                                * plant_df[plant_df.plant_id .== p, :fuel_price][1] for b in B for t in T), digits=2)

                    # fuel_costs += round(1/(B[end]+1) * π * Δt * sum(prod_weights_df[(prod_weights_df.scenario .== s) .& (prod_weights_df.timestep .== t) .& (prod_weights_df.plant_id .== p) .& (prod_weights_df.b .==b), :production][1] * plant_df[plant_df.plant_id .== p, :fuel_price][1] for b in B for t in T), digits=2)
                    # startup_costs += sum(prod_aggregated_df[(prod_aggregated_df.timestep .== t) .& (prod_aggregated_df.plant_id .== p), :startup][1] * C_startup for t in T)
                elseif p in P_h
                    water_value = get_water_value()
                    day_ahead_costs += round(1/(B[end]+1) * π * Δt * sum(first_stage_prod[(first_stage_prod.scenario .== s) .& (first_stage_prod.timestep .== t) .& (first_stage_prod.plant_id .== p) .& (first_stage_prod.b .==b), :production][1] 
                                         * water_value for b in B for t in T), digits=2) #10 plant_df[plant_df.plant_id .== p, :fuel_price][1] for b in B for t in T), digits=2)

                    up_activation_costs += round(1/(B[end]+1) * π * Δt * sum(prod_weights_df[(prod_weights_df.scenario .== s) .& (prod_weights_df.timestep .== t).& (prod_weights_df.plant_id .== p) .& (prod_weights_df.b .==b), :up_activation][1] 
                                         * water_value for b in B for t in T), digits=2) #*  plant_df[plant_df.plant_id .== p, :fuel_price][1] for b in B for t in T), digits=2)

                    down_activation_costs += - round(1/(B[end]+1) * π * Δt * sum(prod_weights_df[(prod_weights_df.scenario .== s) .& (prod_weights_df.timestep .== t).& (prod_weights_df.plant_id .== p) .& (prod_weights_df.b .==b), :down_activation][1] 
                                         * water_value for b in B for t in T), digits=2)#* plant_df[plant_df.plant_id .== p, :fuel_price][1] for b in B for t in T), digits=2)
                    
                    bypass_costs += round(1/(B[end]+1) * π * Δt * sum(C_bypass * hydro_weights_df[(hydro_weights_df.plant_id .== p) .& (hydro_weights_df.scenario .== s) .& (hydro_weights_df.timestep .== t) .& (hydro_weights_df.b .== b), :bypass][1] for t in T for b in B), digits=2)
                    spill_costs += round(1/(B[end]+1) * π * Δt * sum(C_spill * hydro_weights_df[(hydro_weights_df.plant_id .== p) .& (hydro_weights_df.scenario .== s) .& (hydro_weights_df.timestep .== t) .& (hydro_weights_df.b .== b), :spill][1] for t in T for b in B), digits=2)
                    
                    water_costs += round(π *
                        (-hydro_df[(hydro_df.plant_id .== p) .& (hydro_df.timestep .== T[end]) .& (hydro_df.scenario .== s), :volume_res][1] + hydro_data_df[hydro_data_df.plant_id .== p, :starting_reservoir][1]) 
                        * plant_df[plant_df.plant_id .== p, :fuel_price][1], digits=2)
                end
            end
            area_cost_df = DataFrame(
                component = ["Day-ahead costs", "up-activation cost", "down-activation cost", "Shedding", "Dumping", "Water", "Bypass cost", "Spill cost"],
                values = [day_ahead_costs, up_activation_costs, down_activation_costs, shed_costs, dump_costs, water_costs, bypass_costs, spill_costs],
                area = a,
                scenario = s)
            append!(df, area_cost_df)
        end
    end
    append!(df, discrete_results_df[discrete_results_df.component .== "Start", :])
    append!(df, discrete_results_df[discrete_results_df.component .== "Up-reserve", :])
    append!(df, discrete_results_df[discrete_results_df.component .== "Down-reserve", :])
    # display(df)
    # df = round.(df, digits=0)
    CSV.write("continuous_results/continuous_objective.csv", df)
end

function plot_continuous_and_discrete_bounds(bernstein_degree, P)
    continuous_values_df = DataFrame(XLSX.readtable("output/ub_and_lb.xlsx", "Sheet1", infer_eltypes=true))
    # measuring_points_per_hour = continuous_values[(continuous_values.plant_id .== 1) .& (continuous_values.timestep) .& ()]
    discrete_production_df = DataFrame(XLSX.readtable("discrete_results/results.xlsx", "production", infer_eltypes=true))
    s = 0
    # P = unique(continuous_values_df.plant_id)
    if length(P) == 0
        P = unique(continuous_values_df.plant_id)
    end
    steps_per_hour = div(unique(continuous_values_df.timestep)[end],24)
    plotting_timesteps = [0, 10*steps_per_hour]
    column_values = [:lb, :ub]
    for p in P
        df = DataFrame()
        df[!, "timestep"] = continuous_values_df[(continuous_values_df.plant_id .== p) .& (continuous_values_df.scenario .== s), :timestep_fractional]
        sampling_points = div((nrow(df)-1),24)
        for column in column_values
            continuous_ts = continuous_values_df[(continuous_values_df.plant_id .== p) .& (continuous_values_df.scenario .== s), column]
            discrete_ts_low_resolution = discrete_production_df[(discrete_production_df.plant_id .== p) .& (discrete_production_df.scenario .== s), column]
            discrete_ts_high_resolution = change_ts_resolution_to_second(discrete_ts_low_resolution, sampling_points, false)
            # display(discrete_ts_high_resolution)
            pushfirst!(discrete_ts_high_resolution, discrete_ts_high_resolution[1])
            df[!, "$column - discrete"] = discrete_ts_high_resolution
            df[!, "$column - continuous"] = continuous_ts
        end
        df = df[(df.timestep .>= plotting_timesteps[1]) .& (df.timestep .<= plotting_timesteps[2]), :]
        plot_all_columns_df(df, "Power plant $p, bernstein degree $bernstein_degree", "p$p-b$bernstein_degree.png")
    end
end
# plot_continuous_and_discrete_bounds()

# calculate_objective_components_discrete()
# calculate_objective_components_continuous()
# plot_hydro_balance()
# plot_energy_balance()

