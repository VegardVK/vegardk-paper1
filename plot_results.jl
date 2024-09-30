using DataFrames
using XLSX

include("C:/Users/vegardvk/vscodeProjects/bernstein/helper_functions.jl")


function import_data_discrete()
    prod_df = DataFrame(XLSX.readtable("discrete_results/results.xlsx", "production", infer_eltypes=true))
    flow_df = DataFrame(XLSX.readtable("discrete_results/results.xlsx", "transmission", infer_eltypes=true))
    area_df = DataFrame(XLSX.readtable("discrete_results/results.xlsx", "area_results", infer_eltypes=true))
    hydro_df = DataFrame(XLSX.readtable("discrete_results/results.xlsx", "hydro_results", infer_eltypes=true))
    line_df = DataFrame(XLSX.readtable("output/line_data.xlsx", "Sheet1", infer_eltypes=true))

    A = unique(prod_df.area)
    P = unique(prod_df.plant_id)
    T = unique(prod_df.timestep)
    L = unique(flow_df.line_id)
    L_in = Dict(a => unique(line_df[line_df.area_to .== a, :line_id]) for a in A) 
    L_out =  Dict(a => unique(line_df[line_df.area_from .== a, :line_id]) for a in A)

    area_grouped_plants = groupby(prod_df, :area)
    P_a = Dict(g.area[1] => unique(g.plant_id) for g in area_grouped_plants)
    return prod_df, flow_df, area_df, hydro_df, line_df, A, P, T, L, L_in, L_out, P_a
end

function import_data_continuous()
    prod_weights_df = DataFrame(XLSX.readtable("continuous_results/result_weights.xlsx", "production", infer_eltypes=true))
    prod_aggregated_df = DataFrame(XLSX.readtable("continuous_results/results_aggregated.xlsx", "production", infer_eltypes=true))
    area_df = DataFrame(XLSX.readtable("continuous_results/result_weights.xlsx", "area_results", infer_eltypes=true))
    hydro_df = DataFrame(XLSX.readtable("continuous_results/results_aggregated.xlsx", "hydro_results", infer_eltypes=true))

    A = unique(prod_weights_df.area)
    P = unique(prod_weights_df.plant_id)
    T = unique(prod_weights_df.timestep)
    B = unique(prod_weights_df.b)

    return prod_weights_df, prod_aggregated_df, area_df, hydro_df, A, P, T, B
end


function plot_energy_balance()
    prod_df, flow_df, area_df, hydro_df, line_df, A, P, T, L, L_in, L_out, P_a = import_data_discrete()

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
    prod_df, flow_df, area_df, hydro_df, line_df, A, P, T, L, L_in, L_out, P_a = import_data_discrete()
    plant_df = DataFrame(XLSX.readtable("output/plant_data.xlsx", "Sheet1", infer_eltypes=true))
    hydro_data_df = DataFrame(XLSX.readtable("output/hydro_data.xlsx", "Sheet1", infer_eltypes=true))

    C_shedding, C_dumping, C_startup = get_cost_parameters()
    P_t = unique(plant_df[plant_df.fuel_type .== "Thermal", :plant_id])
    P_h = unique(plant_df[plant_df.fuel_type .== "Hydro", :plant_id])

    prod_df_wide = unstack(select(prod_df, :timestep, :production, :plant_id), :plant_id, :production)

    fuel_costs = sum(prod_df[(prod_df.timestep .== t) .& (prod_df.plant_id .== p), :production][1] * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in P_t for t in T)
    startup_costs = sum(prod_df[(prod_df.timestep .== t) .& (prod_df.plant_id .== p), :startup][1] * C_startup for p in P_t for t in T)
    dumping_costs = sum(area_df[(area_df.area .== a) .& (area_df.timestep .== t), :power_dumping][1] * C_dumping for a in A for t in T)
    shed_costs = sum(area_df[(area_df.area .== a) .& (area_df.timestep .== t), :load_shedding][1] * C_shedding for a in A for t in T)
    volume_costs = sum((-hydro_df[(hydro_df.plant_id .== p) .& (hydro_df.timestep .== T[end]), :volume_res][1] + hydro_data_df[hydro_data_df.plant_id .== p, :starting_reservoir][1]) * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in P_h) 
    objective = fuel_costs + startup_costs + dumping_costs + shed_costs + volume_costs

    println("Objective: $objective")
    println("\t Production costs: $fuel_costs")
    println("\t Water costs: $volume_costs")
    println("\t Startup costs: $startup_costs")
    println("\t shed_costs: $shed_costs")
    println("\t dump costs: $dumping_costs")
end

function calculate_objective_components_continuous()
    plant_df = DataFrame(XLSX.readtable("output/plant_data.xlsx", "Sheet1", infer_eltypes=true))
    hydro_data_df = DataFrame(XLSX.readtable("output/hydro_data.xlsx", "Sheet1", infer_eltypes=true))

    prod_weights_df, prod_aggregated_df, area_df, hydro_df, A, P, T, B = import_data_continuous()

    C_shedding, C_dumping, C_startup = get_cost_parameters()
    P_t = unique(plant_df[plant_df.fuel_type .== "Thermal", :plant_id])
    P_h = unique(plant_df[plant_df.fuel_type .== "Hydro", :plant_id])

    fuel_costs = round(1/(B[end]+1) * sum(prod_weights_df[(prod_weights_df.timestep .== t) .& (prod_weights_df.plant_id .== p) .& (prod_weights_df.b .==b), :production][1] * plant_df[plant_df.plant_id .== p, :fuel_price][1]  for p in P_t for b in B for t in T), digits=2)
    shed_costs =  round(1/(B[end]+1) * sum(area_df[(area_df.timestep .== t) .& (area_df.area .== a) .& (area_df.b .== b), :shedding][1] * C_shedding for b in B for a in A for t in T), digits=2)
    dump_costs = round(1/(B[end]+1) * sum(area_df[(area_df.timestep .== t) .& (area_df.area .== a) .& (area_df.b .== b), :dumping][1] * C_dumping for b in B for a in A for t in T), digits=2)
    water_costs = round(sum((-hydro_df[(hydro_df.plant_id .== p) .& (hydro_df.timestep .== T[end]), :volume_res][1] + hydro_data_df[hydro_data_df.plant_id .== p, :starting_reservoir][1]) * plant_df[plant_df.plant_id .== p, :fuel_price][1] for p in P_h), digits=2)
    startup_costs = sum(prod_aggregated_df[(prod_aggregated_df.timestep .== t) .& (prod_aggregated_df.plant_id .== p), :startup][1] * C_startup for p in P_t for t in T)
    objective = fuel_costs + shed_costs + dump_costs + water_costs + startup_costs
    
    println("Objective: $objective")
    println("\t Production costs: $fuel_costs")
    println("\t Water costs: $water_costs")
    println("\t shed_costs: $shed_costs")
    println("\t dump costs: $dump_costs")
    println("\t Startup costs: $startup_costs")
    df = DataFrame(
        component = ["Fuel", "Shedding", "Dumping", "Water", "Start"],
        values = [fuel_costs, shed_costs, dump_costs, water_costs, startup_costs],
        # Objective = objective,
        # Fuel = fuel_costs,
        # Shedding = shed_costs,
        # Dumping = dump_costs,
        # Water = water_costs,
        # Start = startup_costs 
    )
    display(df)
    # df = round.(df, digits=0)
    CSV.write("continuous_results/objective.csv", df)
end

# calculate_objective_components_discrete()
# calculate_objective_components_continuous()
# plot_hydro_balance()
# plot_energy_balance()

