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
        savefig(a, "C:/Users/vegardvk/vscodeProjects/Bernstein/results/$filename")
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
