using Combinatorics
using JuMP
using Plots
using Ipopt
using DataFrames
using CSV
using AxisArrays
using XLSX

include("C:/Users/vegardvk/vscodeProjects/bernstein/helper_functions.jl")

function define_problem(parameters)
    B = parameters["B"]
    T = parameters["T"]
    S = parameters["S"]
    bernstein_degree = parameters["bernstein_degree"]
    
    model = Model()
    set_optimizer(model, Ipopt.Optimizer)

    @variable(model, weights[b in B, t in T] ≥ 0)
    @variable(model, deviation[s in S, t in T])

    @constraint(model, deviation_summation[s in S, t in T], deviation[s, t] == sum(parameters["bernstein_curves"][b+1, s+1] * weights[b, t] for b in B) - parameters["target_values"][s+1, t])
    @constraint(model, continuity_constraint1[t in T[1:end-1]], weights[bernstein_degree, t] == weights[0, t+1])
    @constraint(model, continuity_constraint2[t in T[1:end-1]], weights[bernstein_degree, t] - weights[bernstein_degree-1, t] == weights[1, t+1] - weights[0, t+1])

    
    @objective(model, Min, sum(deviation[s, t]^2 for s in S for t in T))
    optimize!(model)

    weights = value.(weights)
    df = dense_array_to_df(weights)
    
    # df = get_converted_df_separated_timesteps(weights, S[end])
    converted_list = get_converted_list(df, S[end])
    # display(plot(converted_list))

    # plot_all_columns_df(df)

    # print(model)
    return df

end


function define_parameters(target_values, bernstein_degree, sampling_points)
    time_steps = length(target_values)
    T = 1:time_steps
    B = 0:bernstein_degree
    S = 0:sampling_points
    bernstein_curves = zeros(bernstein_degree+1, sampling_points+1)
    target_values_expanded = zeros(sampling_points+1, time_steps)
    for s in S
        for b in B
            bernstein_curves[b+1, s+1] = get_bernstein_val(bernstein_degree, b, s/sampling_points)
        end
        for t in T
            target_values_expanded[s+1, t] = target_values[t]
        end
    end

    data = Dict(
        "bernstein_degree" => bernstein_degree,
        "bernstein_curves" => bernstein_curves,
        "target_values" => target_values_expanded,
        "B" => B,
        "S" => S,
        "T" => T
    )
    return data
end


function get_bernstein_val(B, b, s)
    n_choose_i = binomial(B, b)
    val = n_choose_i * (s^b)*((1-s)^(B-b))
    return val
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


function get_converted_df_separated_timesteps(weights, sampling_points)
    bernstein_degree = length(weights[:,1])-1
    time_steps = length(weights[0, :])
    df = DataFrame()
    for t in 1:time_steps
        # array = zeros(time_steps * sampling_points+1)
        array = fill(NaN, time_steps * sampling_points+1)
        for s in 0:sampling_points
            array[(t-1) * sampling_points+s+1] = 0
            for b in 0:bernstein_degree
                array[(t-1) * sampling_points+s+1] += get_bernstein_val(bernstein_degree, b, s/sampling_points) * weights[b, t]
            end
        end
        df[!, "$t"] = array
        # empty_df.t = array
        # hcat!(empty_df, new_df)
    end
    return df
end


function find_and_write_demand_weights(bernstein_degree, time_steps, areas)
    # file_path = "consumption.xlsx"
    # xf = XLSX.readxlsx(file_path)
    # sh = xf["Sheet1"]
    # demand = sh[2:end, 3]
    # demand = [demand[i] for i in 1:(length(demand))]
    # println(demand)
    load_df = get_load(areas, time_steps)
    weights_df = DataFrame()
    for a in 1:areas
        area_load = load_df[:, "$a"]
        parameters = define_parameters(area_load, bernstein_degree, 100)
        df = define_problem(parameters)
        weights_df = vcat(weights_df, df)
        # df = round.(df, digits=2)
    end
    CSV.write("input/load_weights.csv", weights_df)
end

function find_and_write_inflow_weights(bernstein_degree, time_steps)
    inflow_df = get_inflow(time_steps)
    weights_df = DataFrame()
    plants  = names(inflow_df)
    plants_sorted = sort(plants)
    for p in plants_sorted
        plant_inflow = inflow_df[:, "$p"]
        parameters = define_parameters(plant_inflow, bernstein_degree, 100)
        df = define_problem(parameters)
        df = round.(df, digits = 1)
        weights_df = vcat(weights_df, df)
    end
    CSV.write("input/inflow_weights.csv", weights_df)
end



function find_and_write_capacity_weights(bernstein_degree, n_power_plants)
    file_path = "input/max_production.csv"
    # capacity = [1, 2, 3, 4, 5]
    capacity = [c for c in 1:n_power_plants]
    parameters = define_parameters(capacity, bernstein_degree, 100)
    df = define_problem(parameters)
    df = round.(df, digits=2)
    CSV.write("input/capacity_weights.csv", df)
end

function find_and_write_wind_weights(bernstein_degree, time_steps)
    wind_ts_df = get_wind_ts()
    weights_df = DataFrame()
    plants = names(wind_ts_df) # Rekkefølgen på tidsseriene blir kanskje ikke riktig. Vurder å sorter dem + sorter ved lesing og kobling senere
    for p in plants
        wind_ts = wind_ts_df[1:time_steps, "$p"]
        parameters = define_parameters(wind_ts, bernstein_degree, 100)
        df = define_problem(parameters)
        df = round.(df, digits = 2)
        weights_df = vcat(weights_df, df)
    end
    CSV.write("input/wind_ts_weights.csv", weights_df)
end


nB = 3
# find_and_write_demand_weights(nB)
# find_and_write_capacity_weights()
# find_bernstein_weights()
