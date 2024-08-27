include("C:/Users/vegardvk/vscodeProjects/bernstein/process_input_data.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/discrete_model.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/plot_results.jl")



process_wind_ts_data()
process_load_data()
process_plant_data()

model = define_and_solve_model()
write_results(model)

calculate_objective_components()
