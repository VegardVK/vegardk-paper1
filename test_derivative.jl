using XLSX
using DataFrames
using JuMP
using CPLEX

include("C:/Users/vegardvk/vscodeProjects/bernstein/helper_functions.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/find_bernstein_weights.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/elevation_matrix.jl")

nb = 3

k_matrix = get_k_matrix(nb)
elevation_matrix = get_elevation_matrix(nb, nb+1)

find_and_write_demand_weights(nb)
load_df = DataFrame(XLSX.readtable("output/load_weights.xlsx", "Sheet1", infer_eltypes=true))
T = unique(load_df.timestep)
area_load_weights_df = load_df[load_df.area .== 1, :]

area_load_weights_unstacked_df = select!(unstack(area_load_weights_df, :b, :timestep, :Forbruk), Not(:b))
display(permutedims(area_load_weights_unstacked_df))

b_n_transposed = permutedims(area_load_weights_unstacked_df)
b_derived_transposed = DataFrame()

for i in 1:(nb)
    b_derived_transposed[!, "$i"] = sum(b_n_transposed[!, "x$j"] * k_matrix[j, i] for j in 1:(nb+1))
end

b_derived = permutedims(b_derived_transposed)
rename!(b_derived, string.(1:24))
display(b_derived)

b_n_array = get_converted_list(area_load_weights_unstacked_df, sampling_points)
b_derived_array = get_converted_list(b_derived, sampling_points)

area_results_df = DataFrame(load_B2 = b_n_array, load_derived = b_derived_array, timestep = 0:(1/sampling_points):T[end])
plot_all_columns_df(area_results_df)

# load_df = DataFrame(XLSX.readtable("output/load_weights.xlsx", "Sheet1", infer_eltypes=true))
# T = unique(load_df.timestep)

# sampling_points = 10

# area_load_weights_df = load_df[load_df.area .== 1, :]



# derivative = DataFrame()
# for t in T
#     elem1 = 3* (-area_load_weights_df[(area_load_weights_df.timestep .== t) .& (area_load_weights_df.b .== 0), :Forbruk][1] + area_load_weights_df[(area_load_weights_df.timestep .== t) .& (area_load_weights_df.b .== 1), :Forbruk][1])
#     elem2 = 3* (-area_load_weights_df[(area_load_weights_df.timestep .== t) .& (area_load_weights_df.b .== 1), :Forbruk][1] + area_load_weights_df[(area_load_weights_df.timestep .== t) .& (area_load_weights_df.b .== 2), :Forbruk][1])
#     elem3 = 3* (-area_load_weights_df[(area_load_weights_df.timestep .== t) .& (area_load_weights_df.b .== 2), :Forbruk][1] + area_load_weights_df[(area_load_weights_df.timestep .== t) .& (area_load_weights_df.b .== 3), :Forbruk][1])
#     col = [elem1, elem2, elem3]
#     derivative[!, "$t"] = col
# end
# derivative_array = get_converted_list(derivative, sampling_points)

# area_load_weights_df = select!(unstack(area_load_weights_df, :b, :timestep, :Forbruk), Not(:b))
# area_load_array = get_converted_list(area_load_weights_df, sampling_points)
# area_results_df = DataFrame(load = area_load_array, derivative=derivative_array, timestep = 0:(1/sampling_points):T[end])
# XLSX.writetable("continuous_results/test_derivative.xlsx", "derivative" => area_results_df, overwrite=true)

# display(area_results_df)
# plot_all_columns_df(area_results_df)

