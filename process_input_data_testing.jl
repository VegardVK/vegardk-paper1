include("C:/Users/vegardvk/vscodeProjects/bernstein/process_input_data.jl")
function test_process_input_data()
    steps_per_hour = 4
    scenarios = 2
    process_load_data(steps_per_hour, scenarios)
end

test_process_input_data()