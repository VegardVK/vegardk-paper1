include("C:/Users/vegardvk/vscodeProjects/bernstein/process_input_data.jl")
# include("C:/Users/vegardvk/vscodeProjects/bernstein/continuous_model.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/continuous_model_simulator.jl")
# include("C:/Users/vegardvk/vscodeProjects/bernstein/continuous_model_testing.jl")
# include("C:/Users/vegardvk/vscodeProjects/bernstein/continuous_model_exogenous_prod.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/plot_results.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/find_bernstein_weights.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/continuous_write_results.jl")


steps_per_hour = 1
# process_wind_ts_data(steps_per_hour)
# process_load_data(steps_per_hour)
# process_plant_data(steps_per_hour)

bernstein_degree = 3
sampling_points = 600
timesteps = 24
continuity_constraint = true
column_list = [:ub, :lb]
power_plants = []
scenarios = []
# find_and_write_inflow_weights(bernstein_degree)
# find_and_write_wind_weights(bernstein_degree)
# find_and_write_demand_weights(bernstein_degree)
# find_and_write_production_weights(bernstein_degree, column_list, continuity_constraint, timesteps, sampling_points, power_plants, scenarios)


sampling_points = 60
convert_production_weights_to_values(column_list, sampling_points)
plot_continuous_and_discrete_bounds(bernstein_degree, power_plants)

# find_and_write_shedding_weights(bernstein_degree, false)
# find_and_write_dumping_weights(bernstein_degree, false)


model = define_and_solve_model([0, 1])
# println(model)
write_results(model, sampling_points, [0, 1])
# calculate_objective_components_continuous()