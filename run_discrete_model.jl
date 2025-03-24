include("C:/Users/vegardvk/vscodeProjects/bernstein/process_input_data.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/discrete_model.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/plot_results.jl")


steps_per_hour = 4
scenarios = 10
process_wind_ts_data(steps_per_hour, scenarios)
process_load_data(steps_per_hour, scenarios)
process_plant_data(steps_per_hour, scenarios)

model = define_and_solve_model(steps_per_hour, true)
write_results(model, steps_per_hour)

model = define_and_solve_model(steps_per_hour, false)
write_results(model, steps_per_hour)
calculate_objective_components_discrete()

sampling_points_per_hour = 60
write_expanded_results(sampling_points_per_hour)
