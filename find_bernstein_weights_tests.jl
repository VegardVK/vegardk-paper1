include("C:/Users/vegardvk/vscodeProjects/bernstein/find_bernstein_weights.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/helper_functions.jl")



function test_change_ts_resolution_to_second()
    input_ts = 1:24
    measuring_points = 2
    output_ts = change_ts_resolution_to_second(input_ts, measuring_points)
    println("input_ts, $input_ts")
    println("output_ts, $output_ts")
end
# test_change_ts_resolution_to_second()


function test_find_bernstein_weights()
    target_ts = 1:24
    bernstein_degree = 3
    timesteps = 24
    measuring_points = 2
    cont_constraints = true
    weights_df = find_bernstein_weights(target_ts, bernstein_degree, timesteps, measuring_points, cont_constraints)
    display(weights_df)
end
# test_find_bernstein_weights()

function test_convert_weights_to_long_format()
    bernstein_degree = 3
    timesteps = 24
    test_weights = rand(bernstein_degree+1, timesteps)
    test_weights_df = dense_array_to_df(test_weights)
    display(test_weights_df)
    df = convert_bernstein_to_long_format(test_weights_df)
    display(df)
end
# test_convert_weights_to_long_format()

function test_convert_df_from_weights_to_values()
    target_ts = 1:24
    bernstein_degree = 3
    timesteps = 24
    measuring_points = 3600
    cont_constraints = false
    val_name = "laughs_per_hour"
    weights_df = find_bernstein_weights(target_ts, bernstein_degree, timesteps, val_name, measuring_points, cont_constraints)
    measuring_points = 100
    display(weights_df)
    values_df = convert_df_from_weights_to_values(weights_df, [val_name], measuring_points)
    display(values_df)
    plot_all_columns_df(values_df)
end
# test_convert_df_from_weights_to_values()

