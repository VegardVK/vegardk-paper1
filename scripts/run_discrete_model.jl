include(joinpath(@__DIR__, "process_input_data.jl"))
include(joinpath(@__DIR__, "discrete_model.jl"))
include(joinpath(@__DIR__, "plot_results.jl"))

# Create output directories if they don't exist
for dir in ["output", "continuous_results", "discrete_results", "results"]
    if !isdir(dir)
        mkdir(dir)
    end
end

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
