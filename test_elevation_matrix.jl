include("C:/Users/vegardvk/vscodeProjects/bernstein/helper_functions.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/find_bernstein_weights.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/elevation_matrix.jl")


nb = 3

elevation_matrix = get_elevation_matrix(nb, nb+1)

find_and_write_demand_weights(nb)
load_df = DataFrame(XLSX.readtable("output/load_weights.xlsx", "Sheet1", infer_eltypes=true))
T = unique(load_df.timestep)
area_load_weights_df = load_df[load_df.area .== 1, :]

area_load_weights_unstacked_df = select!(unstack(area_load_weights_df, :b, :timestep, :Forbruk), Not(:b))
display(permutedims(area_load_weights_unstacked_df))

b_n_transposed = permutedims(area_load_weights_unstacked_df)
b_m_transposed = DataFrame()

for i in 1:(nb+2)
    b_m_transposed[!, "$i"] = sum(b_n_transposed[!, "x$j"] * elevation_matrix[j, i] for j in 1:(nb+1))
end

# b_m_transposed[!, "1"] = b_n_transposed[!, "x1"]
# b_m_transposed[!, "2"] = 0.33 * b_n_transposed[!, "x1"] + 0.67 * b_n_transposed[!, "x2"]
# b_m_transposed[!, "3"] = 0.67 * b_n_transposed[!, "x2"] + 0.33 * b_n_transposed[!, "x3"]
# b_m_transposed[!, "4"] = b_n_transposed[!, "x3"]

b_m = permutedims(b_m_transposed)
rename!(b_m, string.(1:24))
display(b_m)

b_n_array = get_converted_list(area_load_weights_unstacked_df, sampling_points)
b_m_array = get_converted_list(b_m, sampling_points)
# area_results_df = DataFrame(load = area_load_array, derivative=derivative_array, timestep = 0:(1/sampling_points):T[end])
area_results_df = DataFrame(load_B2 = b_n_array, load_B3 = b_m_array, timestep = 0:(1/sampling_points):T[end])
plot_all_columns_df(area_results_df)

# XLSX.writetable("continuous_results/test_derivative.xlsx", "derivative" => area_results_df, overwrite=true)


