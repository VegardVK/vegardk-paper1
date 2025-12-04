include(joinpath(@__DIR__, "process_input_data.jl"))
# include(joinpath(@__DIR__, "continuous_model.jl"))
include(joinpath(@__DIR__, "continuous_model_simulator.jl"))
# include(joinpath(@__DIR__, "continuous_model_testing.jl"))
# include(joinpath(@__DIR__, "continuous_model_exogenous_prod.jl"))
include(joinpath(@__DIR__, "plot_results.jl"))
include(joinpath(@__DIR__, "find_bernstein_weights.jl"))
include(joinpath(@__DIR__, "continuous_write_results.jl"))

# Create output directories if they don't exist
for dir in ["output", "continuous_results", "discrete_results", "results"]
    if !isdir(dir)
        mkdir(dir)
    end
end

function process_raw_data(steps_per_hour, scenarios, sampling_points)
    # process_wind_ts_data(steps_per_hour, scenarios)
    process_load_data(12, scenarios, sampling_points)
    # process_plant_data(steps_per_hour, scenarios)
end

function get_input_weights(bernstein_degree, weights_sp, res_sp, timesteps, scenarios)

    # Find weights for wind, load and inflow
    prepare_parameter_weights(bernstein_degree, timesteps, weights_sp) # weights_sp)
    
    # Find weights for upper and lower production bounds
    continuity_constraint = true
    column_list = [:ub, :lb, :production]
    power_plants = []
    scenarios = 0:scenarios
    find_and_write_production_weights(bernstein_degree, column_list, continuity_constraint, timesteps, weights_sp, power_plants, scenarios)
    convert_production_weights_to_values(column_list, res_sp)

    # Plot discrete bounds alongside converted continuous bounds
    # plot_continuous_and_discrete_bounds(bernstein_degree, power_plants)
end

function run_continuous_model(n_scen, sampling_points)
    scenarios = 1:n_scen
    model = define_and_solve_model() #scenarios)
    write_results(model, sampling_points)# , scenarios)
    calculate_objective_components_continuous()
end


steps_per_hour = 4
input_data_scenarios = 10
simulation_scenarios = input_data_scenarios-1
bernstein_degree = 4

weights_calc_sampling_points = 1800
res_sampling_points = 60
continuous_timesteps = 24 * steps_per_hour

process_raw_data(steps_per_hour, input_data_scenarios, res_sampling_points)
# Run discrete model
get_input_weights(bernstein_degree, weights_calc_sampling_points, res_sampling_points, continuous_timesteps, input_data_scenarios)
run_continuous_model(simulation_scenarios, res_sampling_points)

