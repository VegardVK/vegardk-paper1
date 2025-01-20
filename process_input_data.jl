using DataFrames
using XLSX
using CSV
using Interpolations
using Statistics
using Random
using Distributions
using Dates

include("C:/Users/vegardvk/vscodeProjects/bernstein/get_hydro_data.jl")



function process_plant_data(steps_per_hour, scenarios)
    input_df = CSV.read("input/gen.csv", DataFrame)
    input_df[!, :gen_index] = 1:nrow(input_df)

    wind_df = filter(row -> row.Category == "Wind", input_df)
    
    P_t = 1:10
    P_w = wind_df[:, "gen_index"]
    P_a = Dict(
        1 => P_t,
        2 => P_w,
        # 3 => [i for i in 21:30],
    )
    output_df = DataFrame(plant_id=Int[], area=Int[], gen_ub=Float64[], gen_lb=Float64[], fuel_type=String[], fuel_price=Float64[])
    for a in 1:2
        for p in P_a[a]
            p_max = input_df[input_df.gen_index .== p, "PMax MW"][1]
            p_min = input_df[input_df.gen_index .== p, "PMin MW"][1]
            fuel_type = input_df[input_df.gen_index .== p, "Category"][1]
            fuel_price = input_df[input_df.gen_index .== p, "Fuel Price \$/MMBTU"][1]
            if fuel_type != "Wind"
                fuel_type = "Thermal"
            end
            row = (p, a, p_max, p_min, fuel_type, fuel_price)
            push!(output_df, row)
        end
    end

    hydro_df = process_hydro_data(steps_per_hour, scenarios)
    for hydro_row in eachrow(hydro_df)
        output_row = (hydro_row["plant_id"], 3, hydro_row["kap_gen_mw"], 0, "Hydro", hydro_row["fuel_price"])
        push!(output_df, output_row)
    end
    XLSX.writetable("output/plant_data.xlsx", output_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
end

function process_wind_ts_data(steps_per_hour, scenarios)
    input_df = CSV.read("input/gen.csv", DataFrame)
    input_df[!, :gen_index] = 1:nrow(input_df)
    wind_df = filter(row -> row.Category == "Wind", input_df)

    wind_ts_df = get_wind_ts()
    name_mapping = Dict(Symbol(wind_df[!, "GEN UID"][i]) .=> Symbol(wind_df[!, "gen_index"][i]) for i in 1:nrow(wind_df))
    rename!(wind_ts_df, name_mapping)
    wind_ts_df.timestep = 1:nrow(wind_ts_df)
    df = stack(wind_ts_df, Not(:timestep), variable_name="plant_id", value_name="wind_power")
    df.plant_id = parse.(Int, df.plant_id)
    P = unique(df.plant_id)
    extended_df = DataFrame()
    if steps_per_hour != 1
        for p in P        
            extended_plant_df = DataFrame()
            extended_plant_df.wind_power = extend_ts(df[df.plant_id .== p, :wind_power], steps_per_hour)
            extended_plant_df.timestep = 1:(24*steps_per_hour)
            extended_plant_df.plant_id .= p
            append!(extended_df, extended_plant_df)
        end
        df = extended_df

    end
    target_peak_wind = 200
    grouped_df = groupby(df, :timestep)
    group_sums = combine(grouped_df, :wind_power => sum => :total_wind)
    multiplier = (target_peak_wind/maximum(group_sums.total_wind)) 
    df.wind_power .= df.wind_power .* multiplier 
    df = generate_scenarios(df, :wind_power, scenarios)

    XLSX.writetable("output/wind_ts_data.xlsx", df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")

end

function get_wind_ts(year=2020, month=1, day=1)
    file_path = "input/DAY_AHEAD_wind.csv"
    df = CSV.read(file_path, DataFrame)
    filtered_df = subset(df, :Year => ByRow(==(year)), :Month => ByRow(==(month)), :Day => ByRow(==(day)))
    select!(filtered_df, Not([:Year, :Month, :Day, :Period]))
    return filtered_df 
end

function process_hydro_data(steps_per_hour, scenarios)
    file_path = "C:/Users/vegardvk/vscodeProjects/bernstein/Input/moduldata.xlsx"
    # xf = XLSX.readxlsx(file_path)
    # sh = xf["Sheet1"]

    df = DataFrame(XLSX.readtable(file_path, "Sheet1", infer_eltypes=true))
    df = df[df.vassdrag .== "NIDELVA_H.DETD", :]
    df = select(df, :modnr, :modnavn, :kap_mag_mm3, :kap_gen_m3s, :kap_forb_m3s, :kap_gen_mw, 
                :enekv, :topo_gen, :topo_forb, :topo_flom, :tilsig_reg_mm3, :tilsig_ureg_mm3)
    df[!, "kap_spill"] .= 100000
    rename!(df, :modnr => "plant_id", :kap_mag_mm3 => "kap_mag", :kap_forb_m3s => "kap_forb")
    df.kap_mag .*= 1000
    df.kap_forb .*= 1000

    # P = df[:, "plant_id"]
    I_disch, I_spill, I_bypass = find_connected_plants(df)

    end_reservoir = 49900
    energy_price = get_water_value()
    df[!, "fuel_price"] .= 0.0
    set_wv!(I_disch, df, end_reservoir, energy_price)
    add_start_reservoir!(df)
    write_inflow_table(df, 24 * steps_per_hour, scenarios)
    XLSX.writetable("output/hydro_data.xlsx", df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
    return df
end


function process_load_data_old(steps_per_hour, scenarios)
    file_path = "input/consumption.xlsx"
    df = DataFrame(XLSX.readtable(file_path, "Sheet1", infer_eltypes=true))
    multiplier = 200
    df.Forbruk .*= multiplier
    df = select(df, :Forbruk)
    df = df[1:24*3, :]
    df.timestep = repeat(1:24, outer=3)
    df.area = repeat(1:3, inner=24)
    df = select(df, :area, :timestep, :Forbruk)
    
    A = unique(df.area)
    if steps_per_hour != 1
        extended_df = DataFrame()
        for a in A        
            extended_area_df = DataFrame()
            extended_area_df.Forbruk = extend_ts(df[df.area .== a, :Forbruk], steps_per_hour)
            extended_area_df.timestep = 1:(24*steps_per_hour)
            extended_area_df.area .= a
            append!(extended_df, extended_area_df)
        end
    end
    df = extended_df
    display(df)
    df = generate_scenarios(df, :Forbruk, scenarios)
    display(df)

    XLSX.writetable("output/load_data.xlsx", df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
end


function process_load_data(steps_per_hour, scenarios)
    caiso_load_df = CSV.read("input/CAISO load data/CAISO-netdemand-20190101.csv", DataFrame)
    nyiso_load_df = CSV.read("input/NYISO load data/20190101pal.csv", DataFrame, delim=",")
    target_load_peak_hydro = 610
    target_load_peak_thermal = 680

    rename!(nyiso_load_df, "Time Stamp" => "timestamp", "Load" => "load")
    nyiso_load_df.timestamp = map(String, nyiso_load_df.timestamp)
    nyiso_load_df.timestamp = DateTime.(nyiso_load_df.timestamp, dateformat"dd/mm/yyyy HH:MM:SS")
    filter!(row -> second(row.timestamp) == 0, nyiso_load_df)
    nyiso_load_df = nyiso_load_df[nyiso_load_df.Name .== "HUD VL", ["timestamp", "load"]]
    nyiso_load_df.load .= nyiso_load_df.load .* (target_load_peak_hydro/maximum(nyiso_load_df.load))
    nyiso_load_df = change_resolution(nyiso_load_df, steps_per_hour)
    hydro_load_df = DataFrame(
        load = nyiso_load_df.load,
        area = fill(3, nrow(nyiso_load_df)),
        timestep = 1:nrow(nyiso_load_df)
    )
    
    # hydro_load_df.load .= hydro_load_df.load .* (target_load_peak_hydro/maximum(hydro_load_df.load)) 


    rename!(caiso_load_df, "Time" => "timestamp", "Day-ahead net forecast" => "load")
    caiso_load_df.load .= caiso_load_df.load .* (target_load_peak_thermal/maximum(caiso_load_df.load))

    caiso_load_df = change_resolution(caiso_load_df, steps_per_hour)
    thermal_load_df = DataFrame(
        load = caiso_load_df.load,
        timestep = 1:nrow(nyiso_load_df),
        area = fill(1, nrow(caiso_load_df))
    )
    # thermal_load_df.load .= thermal_load_df.load .* (target_load_peak_thermal/maximum(thermal_load_df.load)) 
    
    wind_load_df = DataFrame(
        load = fill(0, nrow(caiso_load_df)),
        timestep = 1:nrow(nyiso_load_df),
        area = fill(2, nrow(caiso_load_df))
    )
    load_df = vcat(thermal_load_df, wind_load_df, hydro_load_df)
    load_df = generate_scenarios(load_df, :load, scenarios)
    # display(load_df)
    XLSX.writetable("output/load_data.xlsx", load_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")

end

function change_resolution(df::DataFrame, steps_per_hour::Int)
    # Validate steps_per_hour
    if !(steps_per_hour in [12, 4, 1])
        throw(ArgumentError("steps_per_hour must be one of 12, 4, or 1"))
    end

    # Number of timesteps to average
    group_size = 12 ÷ steps_per_hour

    # Resample DataFrame by taking averages
    grouped_df = DataFrame(
        timestamp = [df.timestamp[1 + (i - 1) * group_size] for i in 1:(nrow(df) ÷ group_size)],
        load = [mean(df.load[(1 + (i - 1) * group_size):(i * group_size)]) for i in 1:(nrow(df) ÷ group_size)],
    )
    return grouped_df 
end

function extend_ts(ts, n_steps)
    new_ts = ts_interpolation(n_steps, ts)

    for i in 1:length(ts)
        original_avg = ts[i]
        segment = new_ts[((i-1)*n_steps + 1):i*n_steps]
        segment_avg = mean(segment)
        
        scaling_factor = original_avg/segment_avg
        new_ts[((i-1)*n_steps + 1):i*n_steps] .= segment .* scaling_factor
    end
    return new_ts

end

function ts_interpolation(n_steps, input_ts)
    steplength = 1/(n_steps * 2)
    output_vals = []

    for j in 1:(length(input_ts)-1)
        timestep = 0:1
        val = [input_ts[j], input_ts[j+1]]
        itp = LinearInterpolation(timestep, val, extrapolation_bc=Line())
        
        if j == 1
            for i in [steplength + k*2*steplength for k in (-n_steps/2):(n_steps-1)]
                push!(output_vals, max(itp(i), 0))
            end
        elseif j == (length(input_ts)-1)
            for i in [steplength + k*2*steplength for k in 0:((n_steps-1) + (n_steps/2))]
                push!(output_vals, max(itp(i), 0))
            end
        else
            for i in [steplength + k*2*steplength for k in 0:(n_steps-1)]
                push!(output_vals, max(itp(i), 0))
            end
        end
    end
    return output_vals
end


function generate_scenarios(df::DataFrame, column_to_scale::Symbol, num_scenarios::Int)
    # Ensure the column to scale is converted to Float64 for compatibility
    df[!, column_to_scale] = Float64.(df[!, column_to_scale])
    
    # Create an empty DataFrame to store results
    results = copy(df)
    results.scenario = fill(0, nrow(df))

    # Define the normal distribution for scaling
    dist = Normal(1.0, 0.1)
    hard_coded_scaling_wind = [0.99, 0.92, 1.02, 0.93, 0.98, 0.84, 1.01, 0.97, 1.04, 1.25]
    hard_coded_scaling_load = [1.18, 0.82, 1.1, 0.97, 0.92, 1.04, 1.05, 1.1, 0.97, 0.87]

    # Generate scenarios
    for scenario in 1:num_scenarios
        # Create a copy of the DataFrame to keep other columns the same
        scenario_df = copy(df)

        # Scale the specified column
        # scaling_factor = rand(dist)
        if string(column_to_scale) == "wind_power"
            scaling_factor = hard_coded_scaling_wind[scenario]
        elseif string(column_to_scale) == "load"
            scaling_factor = hard_coded_scaling_load[scenario]
        end
        scenario_df[!, column_to_scale] .= df[!, column_to_scale] .* scaling_factor

        # Add the scenario number column
        scenario_df.scenario = fill(scenario, nrow(df))

        # Append to the results DataFrame
        append!(results, scenario_df)
    end
    results[!, column_to_scale] .= round.(results[!, column_to_scale], digits=1)
    return results
end

# Example DataFrame
# df = DataFrame(time_value = 1:10, other_column = repeat(["A", "B"], 5))
# num_scenarios = 2

# # Call the function
# results = generate_scenarios(df, :time_value, num_scenarios)
# println(results)

# a = [1, 2, 5, 6, 2]
# b = generate_scenarios(a, 3)
# display(b)
# process_hydro_data(4)
# process_wind_ts_data(4)
# ts_interpolation(4)
