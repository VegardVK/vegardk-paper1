using Plots
using DataFrames
using JuMP
using AxisArrays


function plot_all_columns_df(df, title="", filename="", y_lim=false)
     
    if "timestep" in names(df)
        index_values = df.timestep
        x_axis_name = "Timestep"
    else
        index_values = 1:nrow(df)
        x_axis_name = ""
    end
    if y_lim
        a = plot(title=title, ylim =(0, :auto))
    else
        a = plot(title=title)
    end
    for col in names(df)
        if col != "timestep"
            plot!(index_values, df[:, col], label=col, xlabel=x_axis_name, legend=:outertopright)
        end
    end
    display(a)
    if filename != ""
        savefig(a, joinpath(@__DIR__, "results", filename))
    end
end

function dense_array_to_df(weights)
    cols = length(weights[1, :])
    array = collect(weights)
    column_names = ["$c" for c in 1:cols]
    df = DataFrame(array, column_names)
    return df
end

function print_model_info(model::Model)
    # Get variables and constraints

    variables = all_variables(model)
    num_binary = 0
    num_integer = 0
    num_continuous = 0
    for var in variables
        if is_binary(var)
            num_binary += 1
        elseif is_integer(var)
            num_integer += 1
        else
            num_continuous += 1
        end
    end
    println("Number of binary variables: ", num_binary)
    println("Number of integer variables: ", num_integer)
    println("Number of continuous variables: ", num_continuous)
    # num_binary = num_constraints(model, VariableRef, MOI.ZeroOne)

    constraints = all_constraints(model, include_variable_in_set_constraints=true)
    println("Number of constraints: ", length(constraints))
    print("\n\n")
end


function get_cost_parameters()
    C_shedding = 1000
    C_dumping = 50
    C_startup = 100
    C_spill = 100
    C_bypass = 100
    return C_shedding, C_dumping, C_startup, C_spill, C_bypass
end

function get_water_value()
    return 10
end

function expand_timeseries(data::DataFrame, s::Int, group_symbols, value_symbols)
    T = unique(data.timestep)
    Δt = 24/T[end]

    timesteps_per_hour = div(T[end], 24)
    repeats_per_datapoint = div(s, timesteps_per_hour)
    expanded_data = DataFrame()
    fractional_timesteps = 0:(1/s):24
    sort!(data, vcat(group_symbols, :timestep))
    for col in names(data)
        expanded_data[!, col] = repeat(data[!, col], inner=repeats_per_datapoint)
    end
    
    for col in value_symbols
        expanded_data[!, col] ./= Δt
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
        results = vcat(results, g)
    end
    results.timestep_fractional .= round.(results.timestep_fractional, digits=4)
    return results
end