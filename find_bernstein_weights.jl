using Combinatorics
using JuMP
using Plots
using Ipopt
using CPLEX

using DataFrames
using CSV
using AxisArrays
using XLSX

include("C:/Users/vegardvk/vscodeProjects/bernstein/helper_functions.jl")

function change_ts_resolution_to_second(input_ts, measuring_points, include_startpoint=true)
    steps_per_hour = div(length(input_ts),24)
    output_ts = []
    for t in 1:24
        hourly_ts = input_ts[Int((t-1)*steps_per_hour+1):Int(t*steps_per_hour)]
        repeated_hourly_ts = repeat(hourly_ts, inner=[div(measuring_points, steps_per_hour)])
        if include_startpoint
            # pushfirst!(repeated_hourly_ts, hourly_ts[1])
            push!(repeated_hourly_ts, hourly_ts[end])
        end
        append!(output_ts, repeated_hourly_ts)
    end
    return output_ts
end

function find_bernstein_weights(target_ts, bernstein_degree, timesteps, val_name="weights", measuring_points_hour=3600, cont_constraints=true, curve_under_target=false, curve_over_target=false)
    T = 1:timesteps
    B = 0:bernstein_degree
    timesteps_per_hour = div(timesteps, 24)
    measuring_points_timestep = div(measuring_points_hour, timesteps_per_hour)
    S = 0:measuring_points_timestep

    target_ts = change_ts_resolution_to_second(target_ts, measuring_points_hour)

    bernstein_curves =Dict((b, s) => get_bernstein_val(bernstein_degree, b, s/measuring_points_timestep) for s in S for b in B)
    model = Model()
    # set_optimizer(model, Ipopt.Optimizer)
    set_optimizer(model, CPLEX.Optimizer)
    
    @variable(model, weights[b in B, t in T] ≥ 0)
    @variable(model, deviation[s in S, t in T])
    @constraint(model, deviation_summation[s in S, t in T], deviation[s, t] == sum(bernstein_curves[b, s] * weights[b, t] for b in B) - target_ts[(t-1)*measuring_points_timestep + s+1])

    if curve_over_target
        @constraint(model, curve_over_target[s in S, t in T], sum(bernstein_curves[b, s] * weights[b, t] for b in B) >= target_ts[(t-1)*measuring_points_timestep + s+1])
    elseif curve_under_target
        @constraint(model, curve_under_target[s in S, t in T], sum(bernstein_curves[b, s] * weights[b, t] for b in B) <= target_ts[(t-1)*measuring_points_timestep + s+1])
    end
    if cont_constraints
        @constraint(model, continuity_constraint1[t in T[1:end-1]], weights[bernstein_degree, t] == weights[0, t+1])
        @constraint(model, continuity_constraint2[t in T[1:end-1]], weights[bernstein_degree, t] - weights[bernstein_degree-1, t] == weights[1, t+1] - weights[0, t+1])
    end
    # @objective(model, Min, sum(deviation[s, t]^2 for s in S for t in T))
    @variable(model, deviation_pos[s in S, t in T] >= 0)
    @variable(model, deviation_neg[s in S, t in T] <= 0)
    @constraint(model, deviation_pos_cons[s in S, t in T], deviation_pos[s, t] >= deviation[s, t])
    @constraint(model, deviation_neg_cons[s in S, t in T], deviation_neg[s, t] <= deviation[s, t])
    @objective(model, Min, sum(deviation_pos[s, t] - deviation_neg[s, t] for s in S, t in T))
    optimize!(model)
 
    weights = value.(weights)
    df = dense_array_to_df(weights)
    df = convert_weights_to_long_format(df, val_name)

    return df
end


function find_bounds_weights(lb_target, ub_target, bernstein_degree, timesteps, output_type, val_name="weights", measuring_points_hour=3600, cont_constraints=true)
    T = 1:timesteps
    B = 0:bernstein_degree
    timesteps_per_hour = div(timesteps, 24)
    measuring_points_timestep = div(measuring_points_hour, timesteps_per_hour)
    S = 0:measuring_points_timestep
    type = 1:2
    lb_target = change_ts_resolution_to_second(lb_target, measuring_points_hour)
    ub_target = change_ts_resolution_to_second(ub_target, measuring_points_hour)

    bernstein_curves =Dict((b, s) => get_bernstein_val(bernstein_degree, b, s/measuring_points_timestep) for s in S for b in B)
    model = Model()
    # set_optimizer(model, Ipopt.Optimizer)
    set_optimizer(model, CPLEX.Optimizer)
    @variable(model, weights[b in B, t in T, i in type] >= 0)
    @variable(model, deviation[s in S, t in T, i in type])

    @constraint(model, deviation_summation_lb[s in S, t in T], deviation[s, t, 1] == sum(bernstein_curves[b, s] * weights[b, t, 1] for b in B) - lb_target[(t-1)*measuring_points_timestep + s+1])
    @constraint(model, deviation_summation_ub[s in S, t in T], deviation[s, t, 2] == sum(bernstein_curves[b, s] * weights[b, t, 2] for b in B) - ub_target[(t-1)*measuring_points_timestep + s+1]) 
    @constraint(model, ub_over_lb[b in B, t in T],weights[b, t, 2] >= weights[b, t, 1])

    if cont_constraints
        @constraint(model, continuity_constraint1[t in T[1:end-1], i in type], weights[bernstein_degree, t, i] == weights[0, t+1, i])
        @constraint(model, continuity_constraint2[t in T[1:end-1], i in type], weights[bernstein_degree, t, i] - weights[bernstein_degree-1, t, i] == weights[1, t+1, i] - weights[0, t+1, i])
    end
    # @objective(model, Min, sum(deviation[s, t]^2 for s in S for t in T))
    @variable(model, deviation_pos[s in S, t in T, i in type] >= 0)
    @variable(model, deviation_neg[s in S, t in T, i in type] <= 0)
    @constraint(model, deviation_pos_cons[s in S, t in T, i in type], deviation_pos[s, t, i] >= deviation[s, t, i])
    @constraint(model, deviation_neg_cons[s in S, t in T, i in type], deviation_neg[s, t, i] <= deviation[s, t, i])
    @objective(model, Min, sum(deviation_pos[s, t, i] - deviation_neg[s, t, i] for s in S for t in T for i in type))
    optimize!(model)

    weights = value.(weights)
    df = dense_array_to_df(weights[:, :, output_type])
    df = convert_weights_to_long_format(df, val_name)

    return df
end


function convert_weights_to_long_format(weights, val_name)
    B = 0:length(weights[:, 1])-1
    df = stack(weights, variable_name="timestep", value_name=val_name)
    df.b = repeat(B, div(nrow(df), length(B)))
    df.timestep = parse.(Int, df.timestep)
    return df
end

function get_bernstein_val(B, b, s)
    n_choose_i = binomial(B, b)
    val = n_choose_i * (s^b)*((1-s)^(B-b))
    return val
end

function convert_df_from_weights_to_values(input_df, column_list, sampling_points)
    B = unique(input_df.b)
    T = unique(input_df.timestep)
    output_df = DataFrame()
    for column in column_list
        column_array = zeros(T[end] * sampling_points + 1)
        for b in B
            column_array[1] += get_bernstein_val(B[end], b, 0) * input_df[(b .== input_df.b) .& (1 .== input_df.timestep), column][1]
            for t in T
                for s in 1:sampling_points
                    column_array[(t-1)*sampling_points+s+1] += get_bernstein_val(B[end], b, s/sampling_points) * input_df[(b .== input_df.b) .& (t .== input_df.timestep), column][1]
                end
            end
        end
        output_df[!, column] = column_array
    end
    return output_df
end

function get_converted_list(weights_df, sampling_points)
    bernstein_degree = length(weights_df[:,1])-1
    time_steps = length(weights_df[1, :])
    converted_list = zeros(time_steps * sampling_points+1)

    for b in 0:bernstein_degree
        converted_list[1] += get_bernstein_val(bernstein_degree, b, 0) * weights_df[b+1, 1]        
        for t in 1:time_steps
            for s in 1:sampling_points
                # println((t-1)*sampling_points+s+1)
                converted_list[(t-1)*sampling_points+s+1] += get_bernstein_val(bernstein_degree, b, s/sampling_points) * weights_df[b+1, t]
            end
        end
    end
    a = 3
    return converted_list
end

function find_and_write_production_weights(bernstein_degree, column_symbols,  cont_constraints, timesteps, measuring_points, P=[], S=[])
    df = DataFrame(XLSX.readtable("discrete_results/results.xlsx", "production", infer_eltypes=true))
    output_df = DataFrame()
    if length(P) == 0
        P = unique(df.plant_id)
    end
    if length(S) == 0
        S = unique(df.scenario)
    end

    for p in P
        for s in S
            first = true
            temp_df = DataFrame(plant_id = repeat([p], timesteps*(bernstein_degree+1)), scenario=repeat([s], timesteps*(bernstein_degree+1)),
                                b = repeat(0:bernstein_degree, timesteps), timestep = repeat(1:timesteps, inner=[bernstein_degree+1]))
            for column_symbol in column_symbols
                value_array = df[(df.plant_id .== p) .& (df.scenario .== s), column_symbol]
                lb_target =  value_array = df[(df.plant_id .== p) .& (df.scenario .== s), :lb]
                ub_target =  value_array = df[(df.plant_id .== p) .& (df.scenario .== s), :ub]

                if string(column_symbol) == "ub"
                    weights_df = find_bounds_weights(lb_target, ub_target, bernstein_degree, timesteps, 2, string(column_symbol), measuring_points, cont_constraints)
                elseif string(column_symbol) == "lb"
                    weights_df = find_bounds_weights(lb_target, ub_target, bernstein_degree, timesteps, 1, string(column_symbol), measuring_points, cont_constraints)
                else
                    weights_df = find_bernstein_weights(value_array, bernstein_degree, timesteps, string(column_symbol), measuring_points, cont_constraints)
                end
                weights_df[!, column_symbol] = round.(weights_df[!, column_symbol], digits=2)
                temp_df[!, column_symbol] = weights_df[!, column_symbol]
            end
            output_df = vcat(output_df, temp_df)
        end
    end
    XLSX.writetable("output/production_weights.xlsx", output_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
end


function convert_production_weights_to_values(column_list, sampling_points=60)
    input_df = DataFrame(XLSX.readtable("output/production_weights.xlsx", "Sheet1", infer_eltypes=true))
    display(input_df)
    P = unique(input_df.plant_id)
    S = unique(input_df.scenario)
    T = unique(input_df.timestep)
    output_df = DataFrame()
    for p in P
        for s in S
            temp_df = convert_df_from_weights_to_values(input_df[(input_df.plant_id .== p) .& (input_df.scenario .== s), :], column_list, sampling_points) 
            temp_df.plant_id .= p
            temp_df.scenario .= s
            temp_df.timestep = vcat([1], [t  for t in T for s in 1:sampling_points])
            temp_df.timestep_fractional=0:(1/sampling_points):T[end]
            append!(output_df, temp_df)
        end
    end
    XLSX.writetable("output/ub_and_lb.xlsx", output_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
end

function find_and_write_demand_weights(bernstein_degree)
    load_df = DataFrame(XLSX.readtable("output/load_data.xlsx", "Sheet1", infer_eltypes=true))
    weights_df = DataFrame()
    A = unique(load_df.area)
    S = unique(load_df.scenario)
    for a in A
        for s in S
            area_load = load_df[(load_df.area .== a) .& (load_df.scenario .== s), :Forbruk]
            parameters = define_parameters(area_load, bernstein_degree, 100)
            df = define_problem(parameters)
            df.b = 0:bernstein_degree
            df_long = stack(df, Not(:b), variable_name="timestep", value_name="Forbruk")
            df_long.timestep = parse.(Int, df_long.timestep)
            df_long.area .= a
            df_long.scenario .= s
            weights_df = vcat(weights_df, df_long)
        end
    end
    weights_df.Forbruk = round.(weights_df.Forbruk, digits=2)
    XLSX.writetable("output/load_weights.xlsx", weights_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
end

function find_and_write_shedding_weights(bernstein_degree, cont_constraints)
    ts_df = DataFrame(XLSX.readtable("discrete_results/results.xlsx", "area_results", infer_eltypes=true))
    weights_df = DataFrame()
    A = unique(ts_df.area)
    for a in A
        area_ts_df = ts_df[ts_df.area .== a, :load_shedding]
        parameters = define_parameters(area_ts_df, bernstein_degree, 100)
        df = define_problem(parameters, cont_constraints)
        df.b = 0:bernstein_degree
        df_long = stack(df, Not(:b), variable_name="timestep", value_name="load_shedding")
        df_long.timestep = parse.(Int, df_long.timestep)
        df_long.area .= a
        weights_df = vcat(weights_df, df_long)
    end
    weights_df.load_shedding = round.(weights_df.load_shedding, digits=2)
    XLSX.writetable("output/shedding_weights.xlsx", weights_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
end

function find_and_write_dumping_weights(bernstein_degree, cont_constraints)
    ts_df = DataFrame(XLSX.readtable("discrete_results/results.xlsx", "area_results", infer_eltypes=true))
    weights_df = DataFrame()
    A = unique(ts_df.area)
    for a in A
        area_ts_df = ts_df[ts_df.area .== a, :power_dumping]
        parameters = define_parameters(area_ts_df, bernstein_degree, 100)
        df = define_problem(parameters, cont_constraints)
        df.b = 0:bernstein_degree
        df_long = stack(df, Not(:b), variable_name="timestep", value_name="power_dumping")
        df_long.timestep = parse.(Int, df_long.timestep)
        df_long.area .= a
        weights_df = vcat(weights_df, df_long)
    end
    weights_df.power_dumping = round.(weights_df.power_dumping, digits=2)
    XLSX.writetable("output/dumping_weights.xlsx", weights_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
end


function find_and_write_wind_weights(bernstein_degree)
    wind_df = DataFrame(XLSX.readtable("output/wind_ts_data.xlsx", "Sheet1", infer_eltypes=true))
    weights_df = DataFrame()
    P = unique(wind_df.plant_id)
    S = unique(wind_df.scenario)
    for p in P
        for s in S
            plant_wind_ts = wind_df[(wind_df.plant_id .== p) .& (wind_df.scenario .== s), :wind_power]
            parameters = define_parameters(plant_wind_ts, bernstein_degree, 100)
            df = define_problem(parameters)
            df.b = 0:bernstein_degree
            df_long = stack(df, Not(:b), variable_name="timestep", value_name="wind_power")
            df_long.timestep = parse.(Int, df_long.timestep)
            df_long.plant_id .= p
            df_long.scenario .= s
            weights_df = vcat(weights_df, df_long)
        end
    end
    weights_df.wind_power = round.(weights_df.wind_power, digits=2)
    XLSX.writetable("output/wind_ts_weights.xlsx", weights_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
end


function find_and_write_inflow_weights(bernstein_degree)
    inflow_df = DataFrame(XLSX.readtable("output/inflow_data.xlsx", "Sheet1", infer_eltypes=true))
    weights_df = DataFrame()
    P = unique(inflow_df.plant_id)
    for p in P
        plant_inflow_ts = inflow_df[inflow_df.plant_id .== p, :inflow]
        parameters = define_parameters(plant_inflow_ts, bernstein_degree, 100)
        df = define_problem(parameters)
        df.b = 0:bernstein_degree
        df_long = stack(df, Not(:b), variable_name="timestep", value_name="inflow")
        df_long.timestep = parse.(Int, df_long.timestep)
        df_long.plant_id .= p
        weights_df = vcat(weights_df, df_long)
    end
    weights_df.inflow = round.(weights_df.inflow, digits=2)
    XLSX.writetable("output/inflow_weights.xlsx", weights_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
end


nB = 3
# find_and_write_demand_weights(nB)
# find_and_write_wind_weights(nB)
# find_and_write_inflow_weights(nB)


# find_and_write_capacity_weights()
# find_bernstein_weights()


# input_ts = [1, 2, 3, 4, 5, 6]
# output_ts = convert_input_ts_for_finding_weights(input_ts)
# println(output_ts)