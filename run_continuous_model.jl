include("C:/Users/vegardvk/vscodeProjects/bernstein/process_input_data.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/continuous_model.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/plot_results.jl")

model = define_and_solve_model()
# write_results(model)