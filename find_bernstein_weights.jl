using Combinatorics
using JuMP
using Plots
using Ipopt
using CPLEX

using DataFrames
using CSV
using AxisArrays
using XLSX

include(joinpath(@__DIR__, "helper_functions.jl"))

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
    set_optimizer_attribute(model, "CPXPARAM_ScreenOutput", 0)

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


function find_bounds_weights(lb_target, ub_target, prod_ref, max_prod, bernstein_degree, timesteps, measuring_points_hour=3600, cont_constraints=true)
    T = 1:timesteps
    B = 0:bernstein_degree
    timesteps_per_hour = div(timesteps, 24)
    measuring_points_timestep = div(measuring_points_hour, timesteps_per_hour)
    S = 0:measuring_points_timestep
    type = 1:3 #1 = lb, 2 = ub, 3 = prod_ref
    lb_target = change_ts_resolution_to_second(lb_target, measuring_points_hour)
    ub_target = change_ts_resolution_to_second(ub_target, measuring_points_hour)
    prod_ref = change_ts_resolution_to_second(prod_ref, measuring_points_hour)
    bernstein_curves =Dict((b, s) => get_bernstein_val(bernstein_degree, b, s/measuring_points_timestep) for s in S for b in B)
    model = Model()
    # set_optimizer(model, Ipopt.Optimizer)
    set_optimizer(model, CPLEX.Optimizer)
    set_optimizer_attribute(model, "CPXPARAM_ScreenOutput", 0)

    @variable(model, weights[b in B, t in T, i in type] >= 0)
    @variable(model, deviation[s in S, t in T, i in type])

    @constraint(model, deviation_summation_lb[s in S, t in T], deviation[s, t, 1] == sum(bernstein_curves[b, s] * weights[b, t, 1] for b in B) - lb_target[(t-1)*measuring_points_timestep + s+1])
    @constraint(model, deviation_summation_ub[s in S, t in T], deviation[s, t, 2] == sum(bernstein_curves[b, s] * weights[b, t, 2] for b in B) - ub_target[(t-1)*measuring_points_timestep + s+1]) 
    @constraint(model, deviation_summation_prod_ref[s in S, t in T], deviation[s, t, 3] == sum(bernstein_curves[b, s] * weights[b, t, 3] for b in B) - prod_ref[(t-1)*measuring_points_timestep + s + 1])
    
    # @constraint(model, ub_over_prod_ref[b in B, t in T],weights[b, t, 2] ≥ weights[b, t, 3])
    # @constraint(model, prod_ref_over_lb[b in B, t in T], weights[b, t, 3] ≥ weights[b, t, 1])
    @constraint(model, ub_over_lb[b in B, t in T], weights[b, t, 2] ≥ weights[b, t, 1])
    ϵ = 0.01
    if max_prod ≥ ϵ
        @constraint(model, ub_under_max_prod[b in B, t in T], weights[b, t, 2] ≤ max_prod - ϵ)
    else
        @constraint(model, ub_under_max_prod[b in B, t in T], weights[b, t, 2] ≤ max_prod)
    end
    if cont_constraints
        @constraint(model, continuity_constraint1[t in T[1:end-1], i in type[1:2]], weights[bernstein_degree, t, i] == weights[0, t+1, i])
        @constraint(model, continuity_constraint2[t in T[1:end-1], i in type[1:2]], weights[bernstein_degree, t, i] - weights[bernstein_degree-1, t, i] == weights[1, t+1, i] - weights[0, t+1, i])
    end
    # @objective(model, Min, sum(deviation[s, t]^2 for s in S for t in T))
    @variable(model, deviation_pos[s in S, t in T, i in type] >= 0)
    @variable(model, deviation_neg[s in S, t in T, i in type] <= 0)
    @constraint(model, deviation_pos_cons[s in S, t in T, i in type], deviation_pos[s, t, i] >= deviation[s, t, i])
    @constraint(model, deviation_neg_cons[s in S, t in T, i in type], deviation_neg[s, t, i] <= deviation[s, t, i])
    @objective(model, Min, sum(deviation_pos[s, t, i] - deviation_neg[s, t, i] for s in S for t in T for i in type))
    optimize!(model)

    weights = value.(weights)
    ub_df = dense_array_to_df(weights[:, :, 2])
    ub_df = convert_weights_to_long_format(ub_df, "ub")

    lb_df = dense_array_to_df(weights[:, :, 1])
    lb_df = convert_weights_to_long_format(lb_df, "lb")
    ub_df.lb = lb_df.lb

    prod_ref_df = dense_array_to_df(weights[:, :, 3])
    prod_ref_df = convert_weights_to_long_format(prod_ref_df, "production")
    ub_df.production = prod_ref_df.production
    # println(df)

    # df[!, val_name] .= floor.(df[!, val_name], digits=2)
    # println(df)
    a = 3
    return ub_df
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

function get_converted_list(weights_df, sampling_points, include_first=true)
    bernstein_degree = length(weights_df[:,1])-1
    time_steps = length(weights_df[1, :])
    converted_list = zeros(time_steps * sampling_points+1)

    for b in 0:bernstein_degree
        if include_first
            converted_list[1] += get_bernstein_val(bernstein_degree, b, 0) * weights_df[b+1, 1]        
        else
            converted_list[end] += get_bernstein_val(bernstein_degree, b, 1) * weights_df[b+1, time_steps]   
        end
        for t in 1:time_steps
            for s in 1:sampling_points
                # println((t-1)*sampling_points+s+1)
                if include_first
                    converted_list[(t-1)*sampling_points+s+1] += get_bernstein_val(bernstein_degree, b, s/sampling_points) * weights_df[b+1, t]
                else
                    converted_list[(t-1)*sampling_points+s] += get_bernstein_val(bernstein_degree, b, (s-1)/sampling_points) * weights_df[b+1, t]
                end
            end
        end
    end
    a = 3
    return converted_list
end

function find_and_write_production_weights(bernstein_degree, column_symbols,  cont_constraints, timesteps, measuring_points, P=[], S=[])
    df = DataFrame(XLSX.readtable("discrete_results/results.xlsx", "first_stage", infer_eltypes=true))
    plant_data_df = DataFrame(XLSX.readtable("output/plant_data.xlsx", "Sheet1", infer_eltypes=true))
    output_df = DataFrame()
    if length(P) == 0
        P = unique(df.plant_id)
    end
    if length(S) == 0
        S = unique(df.scenario)
    end
    wind_df = DataFrame(XLSX.readtable("output/wind_power_weights.xlsx", "Sheet1", infer_eltypes=true))
    P_w = unique(plant_data_df[plant_data_df.fuel_type .== "Wind", :plant_id])

    for p in P

        ub = plant_data_df[plant_data_df.plant_id .== p, "gen_ub"][1]
        lb_target = df[(df.plant_id .== p), :lb]
        ub_target = df[(df.plant_id .== p), :ub]
        if p in P_w
            lb_target .= 0
            cont_constraints = false
        end
        prod_ref = df[(df.plant_id .== p), :first_stage_prod]

        weights_df = find_bounds_weights(lb_target, ub_target, prod_ref, ub, bernstein_degree, timesteps, measuring_points, cont_constraints)

        for s in S
            temp_df = DataFrame(plant_id = repeat([p], timesteps*(bernstein_degree+1)), scenario=repeat([s], timesteps*(bernstein_degree+1)),
                                b = repeat(0:bernstein_degree, timesteps), timestep = repeat(1:timesteps, inner=[bernstein_degree+1]))

            println("Finding bounds for power plant $p in scenario $s")
            temp_df[!, :ub] = weights_df[!, :ub]
            temp_df[:, :lb] = weights_df[!, :lb]
            temp_df[:, :production] = weights_df[!, :production]
            # if Symbol("production") in column_symbols
            #     prod_array = df[(df.plant_id .== p) .& (df.scenario .== s), :production]
            #     weights_df = find_bernstein_weights(prod_array, bernstein_degree, timesteps, "production", measuring_points, cont_constraints)
            #     temp_df[!, :production] = weights_df[!, :production]
            # end
            output_df = vcat(output_df, temp_df)
        end
    end
    XLSX.writetable("output/production_weights.xlsx", output_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
end


function convert_production_weights_to_values(column_list, sampling_points=60)
    input_df = DataFrame(XLSX.readtable("output/production_weights.xlsx", "Sheet1", infer_eltypes=true))
    # display(input_df)
    P = unique(input_df.plant_id)
    S = unique(input_df.scenario)
    T = unique(input_df.timestep)
    output_df = DataFrame()
    sp_ts_input = div(sampling_points, div(T[end], 24))
    for p in P
        for s in S
            temp_df = convert_df_from_weights_to_values(input_df[(input_df.plant_id .== p) .& (input_df.scenario .== s), :], column_list, sp_ts_input) 
            temp_df.plant_id .= p
            temp_df.scenario .= s
            temp_df.timestep = vcat([1], [t  for t in 1:24 for s in 1:sampling_points])
            temp_df.timestep_fractional=0:(1/sampling_points):24
            append!(output_df, temp_df)
        end
    end
    output_df.timestep_fractional .= round.(output_df.timestep_fractional, digits=4)
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


function find_and_write_inflow_weights(bernstein_degree, timesteps, sampling_points)
    inflow_df = DataFrame(XLSX.readtable("output/inflow_data.xlsx", "Sheet1", infer_eltypes=true))
    weights_df = DataFrame()
    P = unique(inflow_df.plant_id)
    for p in P
        plant_inflow_ts = inflow_df[inflow_df.plant_id .== p, :inflow]
        df = find_bernstein_weights(plant_inflow_ts, bernstein_degree, timesteps, "inflow", sampling_points)
        df.plant_id .= p
        weights_df = vcat(weights_df, df)
    end
    weights_df.inflow = round.(weights_df.inflow, digits=2)
    display(weights_df)
    XLSX.writetable("output/inflow_weights.xlsx", weights_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
end

function prepare_parameter_weights(nB, nT, nS)
    # find_and_write_parameter_weights(nB, nT, nS, "output/inflow_data.xlsx", :inflow, :plant_id)
    # find_and_write_parameter_weights(nB, nT, nS, "output/wind_ts_data.xlsx", :wind_power, :plant_id)
    find_and_write_parameter_weights(nB, nT, nS, "output/5min_load_data.xlsx", :load, :area)
end

function find_and_write_parameter_weights(nB, nT, nS, input_file, val_col, idx_col)
    input_df = DataFrame(XLSX.readtable(input_file, "Sheet1", infer_eltypes=true))
    weights_df = DataFrame()
    idx_set = unique(input_df[:, idx_col])
    # display(input_df)
    S = unique(input_df.scenario)
    for i in idx_set
        for s in S
            input_ts = input_df[(input_df[:, idx_col] .== i) .& (input_df[:, :scenario] .== s), val_col]
            # input_ts ./= div(length(input_ts), nT)
            df = find_bernstein_weights(input_ts, nB, nT, string(val_col), nS)
            df[:, idx_col] .= i
            df.scenario .= s
            weights_df = vcat(weights_df, df)
        end
    end
    weights_df[:, val_col] = round.(weights_df[:, val_col], digits=2)
    println("Writing weights for $val_col")
    # display(weights_df)
    output_filename = "output/" * string(val_col) * "_weights.xlsx"
    XLSX.writetable(output_filename, weights_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
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