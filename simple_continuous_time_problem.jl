using DataFrames
using JuMP
using CSV
# using Ipopt
# using GLPK
using CPLEX
using AxisArrays


include("C:/Users/vegardvk/vscodeProjects/bernstein/find_bernstein_weights.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/helper_functions.jl")


function define_parameters(bernstein_degree, time_steps, power_plants, sampling_points, areas)
    demand_weights_df = CSV.read("input/load_weights.csv", DataFrame)
    demand_array = zeros(bernstein_degree + 1, areas, time_steps)
    for a in 1:areas
        demand_array[:, a, :] .= demand_weights_df[(1+(a - 1)*(bernstein_degree + 1)):a*(bernstein_degree+1), :]
    end

    L_in = Dict(
        1 => [],
        2 => [1],
        3 => [2, 3]
    )
    L_out = Dict(
        1 => [1, 2],
        2 => [3],
        3 => []
    )
    cap_line = Dict(
        1 => 50,
        2 => 50,
        3 => 50
    )

    data = Dict(
        # "nP" => power_plants,
        "nB" => bernstein_degree,
        "nT" => time_steps,
        "nA" => areas,
        # "P_a" => P_a,
        # "P_t" => P_t,
        "cap_line" => cap_line,
        "L" => 3,
        "L_in" => L_in,
        "L_out" => L_out,
        "sampling_points" => sampling_points,
        # "C" => Dict(1 => 1, 2 => 2, 3=> 999, 4=> 1000),
        # "C_shedding" => C_shedding,
        # "C_dumping" => C_dumping,
        # "C_startup" => C_startup,

        # "P_max" => Dict(1 => 1, 2 => 1, 3=> 1000),
        "demand_weights" => demand_array[:, :, :],
    )

    add_power_plant_data!(data)
    get_module_data!(data)
    # data["P_a"][3] = data["P_h"]

    inflow_weights_df = CSV.read("input/inflow_weights.csv", DataFrame)
    hydro_plants = data["P_h"]
    hydro_plants = sort(hydro_plants)
    inflow_dict = Dict()
    for (count, p) in enumerate(hydro_plants)
        inflow_dict[p] = inflow_weights_df[(1+(count - 1)*(bernstein_degree + 1)):count*(bernstein_degree+1), :]
    end
    data["inflow_weights"] = inflow_dict

    wind_weights_df = CSV.read("input/wind_ts_weights.csv", DataFrame)
    wind_plants = data["P_w"]
    # wind_plants = sort(wind_plants)
    wind_dict = Dict()
    for (count, p) in enumerate(wind_plants)
        wind_dict[p] = wind_weights_df[(1+(count-1)*(bernstein_degree+1)):count*(bernstein_degree+1), :]
    end
    data["wind_weights"] = wind_dict
    return data
end

function define_model(bernstein_degree, time_steps, nP, sampling_points, areas)
    p_dict = define_parameters(bernstein_degree, time_steps, nP, sampling_points, areas)
    model = Model()
    # set_optimizer(model, Ipopt.Optimizer)
    # set_optimizer(model, GLPK.Optimizer)
    set_optimizer(model, CPLEX.Optimizer)

    P = p_dict["P"]
    P_h = p_dict["P_h"]
    P_t = p_dict["P_t"]
    P_w = p_dict["P_w"]

    P_a = p_dict["P_a"]
    B = 0:p_dict["nB"]
    B_I = 0:(p_dict["nB"]+1)
    T = 1:p_dict["nT"]
    T_extended = 1:(p_dict["nT"] + 1)
    A = 1:p_dict["nA"]
    L = 1:p_dict["L"]
    L_in = p_dict["L_in"]
    L_out = p_dict["L_out"]
    I_bypass = p_dict["I_bypass"] # Dictionary where I_bypass[p] gives set of powerplants that bypass into p
    I_disch = p_dict["I_disch"]
    I_spill = p_dict["I_spill"]

    @variable(model, prod[p in P, b in B, t in T] ≥ 0)
    @variable(model, transmission[l in L, b in B, t in T])

    @variable(model, shedding[b in B, a in A, t in T] ≥ 0)
    @variable(model, dumping[b in B, a in A, t in T] ≥ 0)

    @variable(model, flow_disch[p in P_h, b in B, t in T] ≥ 0) # Antar at alle moduler bare har ett discharge-segment
    @variable(model, flow_bypass[p in P_h, b in B, t in T] ≥ 0)
    @variable(model, flow_spill[p in P_h, b in B, t in T] ≥ 0)
    @variable(model, total_flow_in[p in P_h, b in B, t in T] ≥ 0)
    @variable(model, total_flow_out[p in P_h, b in B, t in T] ≥ 0)
    @variable(model, net_inflow[p in P_h, b in B, t in T])
    @variable(model, volume_end[p in P_h, t in 0:T[end]] ≥ 0)
    @variable(model, volume[p in P_h, b in B_I, t in T])

    @constraint(model, wind_production[p in P_w, b in B, t in T], prod[p, b, t] == p_dict["wind_weights"][p][b+1, t])

    @constraint(model, controlled_inflow[p in P_h, b in B, t in T], total_flow_in[p, b, t] == sum(flow_disch[i, b, t] for i in I_disch[p]) 
                                                                                            + sum(flow_bypass[i, b, t] for i in I_bypass[p]) 
                                                                                            + sum(flow_spill[i, b, t] for i in I_spill[p]))

    @constraint(model, controlled_outflow[p in P_h, b in B, t in T], total_flow_out[p, b, t] == flow_disch[p, b, t] + flow_bypass[p, b, t] + flow_spill[p, b, t])

    @constraint(model, hydro_production[p in P_h, b in B, t in T], prod[p, b, t] == p_dict["enekv"][p] * flow_disch[p, b, t])

    @constraint(model, bypass_ub[p in P_h, b in B, t in T], flow_bypass[p, b, t] ≤ p_dict["kap_forb"][p])
    @constraint(model, prod_ub[p in P_h, b in B, t in T], flow_disch[p, b, t] ≤ p_dict["kap_gen"][p])
    @constraint(model, spill_ub[p in P_h, b in B, t in T], flow_spill[p, b, t] ≤ p_dict["kap_spill"][p])

    @constraint(model, net_inflow_calculation[p in P_h, b in B, t in T], net_inflow[p, b, t] == total_flow_in[p, b, t] 
                                                                                                + p_dict["inflow_weights"][p][b+1, t] 
                                                                                                - total_flow_out[p, b, t])

    @constraint(model, starting_reservoir[p in P_h], volume_end[p, 0] == p_dict["starting_reservoir"][p])
    @constraint(model, volume_calculation[p in P_h, t in T], volume_end[p, t] - volume_end[p, t-1] == (1/(B[end]+1)) * sum(net_inflow[p, b, t] for b in B))
    @constraint(model, simple_vol_ub[p in P_h, t in T], volume_end[p, t] ≤ p_dict["kap_mag"][p])

    @constraint(model, volume_weights[p in P_h, b in B_I, t in T], volume[p, b, t] == (1/(B[end]+1)) * sum(net_inflow[p, b2, t] for b2 in B[1:b]))
    @constraint(model, volume_ub[p in P_h, b in B_I, t in T], volume_end[p, t-1] + volume[p, b, t] ≤ p_dict["kap_mag"][p])
    # @constraint(model, volume_lb[p in P_h, b in B, t in T], volume[p, ])
    # @constraint(model, vol_ub[p in P_h, b in B, t in T], volume[p, b, t] ≤ p_dict["kap_mag"][p])


    @variable(model, status[p in P_t, t in T_extended], Bin)
    @variable(model, startup[p in P_t, t in T_extended], Bin)
    @constraint(model, startup_count[p in P_t, t in T_extended[2:end]], startup[p, t] ≥ status[p, t] - status[p, t-1])

    @constraint(model, unit_commitment1a[p in P_t, t in T], prod[p, 0, t] ≤ p_dict["gen_ub"][p] * status[p, t])
    @constraint(model, unit_commitment1b[p in P_t, t in T], prod[p, 0, t] ≥ p_dict["gen_lb"][p] * status[p, t])

    @constraint(model, unit_commitment2a[p in P_t, t in T], prod[p, B[end], t] ≤ p_dict["gen_ub"][p] * status[p, t+1])
    @constraint(model, unit_commitment2b[p in P_t, t in T], prod[p, B[end], t] ≥ p_dict["gen_lb"][p] * status[p, t+1])

    @constraint(model, unit_commitment3a[p in P_t, t in T], prod[p, 1, t] ≤ p_dict["gen_ub"][p] * status[p, t])
    @constraint(model, unit_commitment3b[p in P_t, t in T], prod[p, 1, t] ≥ p_dict["gen_lb"][p] * status[p, t])

    @constraint(model, unit_commitment4a[p in P_t, t in T], prod[p, B[end]-1, t] ≤ p_dict["gen_ub"][p] * status[p, t+1])
    @constraint(model, unit_commitment4b[p in P_t, t in T], prod[p, B[end]-1, t] ≥ p_dict["gen_lb"][p] * status[p, t+1])

    @constraint(model, continuity_constraint1[p in P_t, t in T[1:end-1]], prod[p, bernstein_degree, t] == prod[p, 0, t+1])
    @constraint(model, continuity_constraint2[p in P_t, t in T[1:end-1]], prod[p, bernstein_degree, t] - prod[p, bernstein_degree-1, t] == prod[p, 1, t+1] - prod[p, 0, t+1])
    
    @constraint(model, energy_balance[b in B, t in T, a in A], 
                sum(prod[p, b, t] for p in P_a[a])
                + sum(transmission[l, b, t] for l in L_in[a])
                - sum(transmission[l, b, t] for l in L_out[a])
                + shedding[b, a, t]
                - dumping[b, a, t] 
                == p_dict["demand_weights"][b+1, a, t])
    # @constraint(model, min_production[p in P, b in B, t in T], prod[p, b, t] >= 0)

    # @constraint(model, max_production[p in P, t in T, b in B], prod[p, b, t]  <=  p_dict["gen_ub"][p])
    @constraint(model, max_production[p in union(P_t, P_h), t in T, b in B], prod[p, b, t]  <=  p_dict["gen_ub"][p])
    @constraint(model, max_transmission[l in L, b in B, t in T], transmission[l, b, t] ≤ p_dict["cap_line"][l])
    @constraint(model, min_transmission[l in L, b in B, t in T], transmission[l, b, t] ≥ -p_dict["cap_line"][l])

    # @objective(model, Min, sum(prod[p, b, t] * get_bernstein_val(bernstein_degree, b, s/sampling_points) * p_dict["C"][p] for p in P for b in B for t in T for s in 0:sampling_points))
    @objective(model, Min, 1/(B[end]+1) * (sum(prod[p, b, t] * p_dict["fuel_price"][p]  for p in P_t for b in B for t in T)
                + sum(shedding[b, a, t] * p_dict["C_shedding"] + dumping[b, a, t] * p_dict["C_dumping"] for b in B for a in A for t in T))
                + sum((-volume_end[p, T[end]] + volume_end[p, 0]) * p_dict["fuel_price"][p] for p in P_h)
                + sum(startup[p, t] * p_dict["C_startup"] for p in P_t for t in T))

    # println(model)
    # print_model_info(model)

    optimize!(model)
    return model, p_dict
end

function plot_transmission_flow(line_dict)
    # println(line_dict)
    df = DataFrame()
    for (key, val) in line_dict
        df[!, "$key"] = val
    end

    plot_all_columns_df(df, "Flow on transmission lines", "Continuous_flow_transmission_lines.png")
end


function plot_model_results(model, p_dict)
    # P = 1:p_dict["nP"]
    # P_a = p_dict["P_a"]
    A = 1:p_dict["nA"]
    sampling_points = p_dict["sampling_points"]

    line_dict = Dict()
    transmission_weights = value.(model[:transmission][:, :, :])
    for l in 1:p_dict["L"]
        line_df = dense_array_to_df(transmission_weights[l, :,:])
        line_list = get_converted_list(line_df, sampling_points)    
        line_dict[l] = line_list
    end
    plot_transmission_flow(line_dict)

    list_length = length(line_dict[1])
    for a in A
        df = DataFrame()
        temp_df = DataFrame()
        temp_df[!, "Export"] = zeros(list_length)
        temp_df[!, "Import"] = zeros(list_length)
        
        for l in p_dict["L_in"][a]
            temp_df[!, "Import"] .+= [min(0, val) for val in line_dict[l]]
            temp_df[!, "Export"] .-= [max(0, val) for val in line_dict[l]]
        end

        for l in p_dict["L_out"][a]
            temp_df[!, "Export"] .+= [min(0, val) for val in line_dict[l]]
            temp_df[!, "Import"] .-= [max(0, val) for val in line_dict[l]]
        end
        df[!, "Net position"] = temp_df[!, :Export] - temp_df[!, :Import]
        
        # demand_df = p_dict["demand_weights"]
        # demand_array = get_converted_list(demand_df, sampling_points)
        # df[!, "Demand"] = demand_array

        demand_array = p_dict["demand_weights"]
        demand_list = get_converted_list(demand_array[:, a, :], sampling_points)
        df[!, "Demand"] = demand_list

        shedding_weights = value.(model[:shedding][:,a,:])
        shedding_weights_df = dense_array_to_df(shedding_weights)
        shedding_array = get_converted_list(shedding_weights_df, sampling_points)
        df[!, "Shedding"] = shedding_array

        dumping_weights = value.(model[:dumping][:,a, :])
        dumping_weights_df = dense_array_to_df(dumping_weights)
        dumping_array = get_converted_list(dumping_weights_df, sampling_points)
        df[!, "Dumping"] = dumping_array

        prod_df = DataFrame()
        for p in p_dict["P_a"][a]
            power_plant_weights = value.(model[:prod][p, :, :])
            power_plant_weights_df = dense_array_to_df(power_plant_weights)
            prod_array = get_converted_list(power_plant_weights_df, sampling_points)
            prod_df[!, "Power plant $p"] = prod_array
            # df[!, "Power plant $p"] = prod_array
        end
        df[!, "Production"] = sum(eachcol(prod_df)) 

        # df = round.(df, digits=2)
        # df[!, "Diff"] = df[!, "Production"] - df[!, "Demand"]

        # display(df)
        plot_all_columns_df(df, "Continuous model, area $a", "continous_area$a.png")


    end
    println(objective_value(model))

end


function plot_hydro_balance(model, bernstein_degree, time_steps, nP, sampling_points, areas)
    p_dict = define_parameters(bernstein_degree, time_steps, nP, sampling_points, areas)
    T = 1:time_steps

    volume = value.(model[:volume_end])
    total_flow_in = value.(model[:total_flow_in])
    total_flow_out = value.(model[:total_flow_out])
    volume_change = value.(model[:volume])
    net_inflow = value.(model[:net_inflow])

    all_reservoirs_df = DataFrame()
    for p in p_dict["P_h"]
        volume_change_weights_df = dense_array_to_df(volume_change[p, :, :])
        volume_change_array = get_converted_list(volume_change_weights_df, sampling_points)

        net_inflow_weights_df = dense_array_to_df(net_inflow[p, :, :])
        net_inflow_array = get_converted_list(net_inflow_weights_df, sampling_points)

        single_reservoir_df = DataFrame()
        volume_array = zeros(1 + sampling_points*time_steps)
        volume_array2 = zeros(1 + sampling_points*time_steps)
        volume_array[1] = volume[p, 0]
        volume_array2[1] = volume[p, 0]
        index = 2

        for t in 1:time_steps
            for s in 1: sampling_points
                volume_array[index] = volume[p, t-1] + volume_change_array[index]
                volume_array2[index] = volume[p, t]
                index += 1
            end
        end

        # df[!, "volume"] = collect(volume[p, :])
        single_reservoir_df[!, "volume_continuous"] = volume_array
        single_reservoir_df[!, "volume_discrete"] = volume_array2
        all_reservoirs_df[!, "Reservoir $p"] = volume_array

        flow_in_weights_df = dense_array_to_df(total_flow_in[p, :, :])
        flow_in_array = get_converted_list(flow_in_weights_df, sampling_points)
        # single_reservoir_df[!, "total_flow_in"] = flow_in_array

        flow_out_weights_df = dense_array_to_df(total_flow_out[p, :, :])
        flow_out_array = get_converted_list(flow_out_weights_df, sampling_points)
        # single_reservoir_df[!, "total_flow_out"] = flow_out_array
        if p == 49927
            println("p: $p")
            println("\tVolume array: ", volume_array2)
            println("\tVolume change array: ", round.(volume_change_array, digits=0))
            println("\tNet inflow array: ", round.(net_inflow_array, digits=0))
            plot_all_columns_df(single_reservoir_df, "Continuous model, hydro balance, plant $p", "continuous_hydro_balance$p.png")
        end
    end
    plot_all_columns_df(all_reservoirs_df, "Continuous model, reservoirs", "continuous_all_reservoirs.png")
end

function print_simple_results(model, bernstein_degree, time_steps, nP, sampling_points, areas)
    p_dict = define_parameters(bernstein_degree, time_steps, nP, sampling_points, areas)

    A = 1:p_dict["nA"]
    P_t = p_dict["P_t"]
    P_h = p_dict["P_h"]
    B = 0:p_dict["nB"]
    T = 1:p_dict["nT"]
    L = 1:p_dict["L"]

    prod = value.(model[:prod])
    shedding = value.(model[:shedding])
    dumping = value.(model[:dumping])
    volume_end = value.(model[:volume_end])
    startup = value.(model[:startup])

    line_dict = Dict()
    transmission_weights = value.(model[:transmission][:, :, :])
    for l in 1:p_dict["L"]
        line_df = dense_array_to_df(transmission_weights[l, :,:])
        line_list = get_converted_list(line_df, sampling_points)    
        line_dict[l] = line_list
    end

    println("Transmission:")
    for l in L
        println("\tLine $l: ")
        println("\t\tPositive direction: ", sum([t for t in line_dict[l] if t>0]))
        println("\t\tNegative direction: ", sum([t for t in line_dict[l] if t<0]))
    end

    fuel_costs = sum((1/(B[end]+1)) * prod[p, b, t] * p_dict["fuel_price"][p] for p in P_t for b in B for t in T)
    shedding_costs = sum(1/(B[end]+1) * shedding[b, a, t] * p_dict["C_shedding"] for b in B for a in A for t in T)
    dumping_costs = sum(1/(B[end]+1) * dumping[b, a, t] * p_dict["C_dumping"] for b in B for a in A for t in T)
    volume_change_costs = sum((-volume_end[p, T[end]] + volume_end[p, 0]) * p_dict["fuel_price"][p] for p in P_h)
    startup_costs = sum(startup[p, t] * p_dict["C_startup"] for p in P_t for t in T)
    objective = objective_value(model)


    println("Objective: $objective")
    println("\tTotal fuel costs: $fuel_costs" )
    println("\tStartup costs: $startup_costs")
    println("\tShedding costs: $shedding_costs")
    println("\tDumping costs: $dumping_costs")
    println("\tVolume change costs: $volume_change_costs")


end

function plot_unit_commitment(model)
    p_dict = define_parameters()
    T = p_dict["T"]
    thermal_uc = value.(model[:status])

    for a in p_dict["A"]
        df = DataFrame()
        peak = zeros(size(thermal_uc, 2))
        for p in p_dict["P_a"][a]
            if !(p in p_dict["P_t"]) 
                continue
            end
            uc = collect(thermal_uc[p, :])
            if sum(uc) != 0
                df[!, "$p"] = uc + peak
                peak .+= uc
            end
        end
        if ncol(df) > 0
            plot_all_columns_df(df, "Continuous model, unit commitment, area $a", "continous_uc_area$a.png", true)
        end
    end
end

function main()
    nT = 24
    nB = 3
    nP = 15
    nA = 3
    nW = 4
    sampling_points = 10

    # find_and_write_demand_weights(nB, nT, nA)
    # find_and_write_inflow_weights(nB, nT)
    # find_and_write_capacity_weights(1, nP)
    # find_and_write_wind_weights(nB, nT)

    model, p_dict = define_model(nB, nT, nP, sampling_points, nA)
    plot_model_results(model, p_dict)

    print_simple_results(model, nB, nT, nP, sampling_points, nA)
    # plot_hydro_balance(model, nB, nT, nP, sampling_points, nA)
    # plot_unit_commitment(model)
end
main()