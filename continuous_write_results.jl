using Statistics
using UnPack

# include("C:/Users/vegardvk/vscodeProjects/bernstein/continuous_model.jl")


function write_results(model, sampling_points, S_input=[])
    input_parameters, input_sets = read_input_data(S_input)
    @unpack wind_df, hydro_df, inflow_df, load_df, plant_df, line_df, prod_df, shedding_df, dumping_df = input_parameters
    @unpack A, P, T, L, B, P_w, P_t, P_h, P_a, L_in, L_out, S = input_sets

    model_results = (
        production = value.(model[:prod]),
        up_activation = value.(model[:up_activation]),
        down_activation = value.(model[:down_activation]),
        uc = value.(model[:status]),
        startup = value.(model[:startup]),
        shedding_weights = value.(model[:shedding]),
        dumping_weights = value.(model[:dumping]),
        transmission = value.(model[:transmission]),
        discharge = value.(model[:flow_disch]),
        spill = value.(model[:flow_spill]),
        bypass = value.(model[:flow_bypass]),
        total_inflow = value.(model[:total_flow_in]),
        total_outflow = value.(model[:total_flow_out]),
        volume_end = value.(model[:volume_end]),
        volume = value.(model[:volume]),
    )

    write_aggregated_results(model_results, input_sets, input_parameters)
    write_expanded_results(model_results, input_sets, input_parameters, sampling_points)
    write_weight_results(model_results, input_sets, input_parameters)
    # frequency = value.(model[:frequency])
    # frequency_d = value.(model[:frequency_d])
    # frequency_expanded_df = get_expanded_frequency_result(sampling_points, frequency, frequency_d, T)
    # CSV.write("continuous_results/frequency.csv", frequency_expanded_df)
end

function write_expanded_results(model_results, input_sets, input_parameters, sampling_points)
    hydro_expanded_df = get_expanded_hydro_results(model_results, input_sets, sampling_points)
    CSV.write("continuous_results/hydro_results_expanded.csv", hydro_expanded_df)
    
    production_expanded_df = get_expanded_prod_results(model_results, input_sets, input_parameters, sampling_points)
    CSV.write("continuous_results/production_expanded.csv", production_expanded_df)

    area_expanded_df = get_expanded_area_results(model_results, input_sets, input_parameters, sampling_points)
    CSV.write("continuous_results/area_results_expanded.csv", area_expanded_df)

    transmission_expanded_df = get_expanded_transmission_results(model_results, input_sets, input_parameters, sampling_points)
    CSV.write("continuous_results/transmission_results_expanded.csv", transmission_expanded_df)
    XLSX.writetable("continuous_results/results_expanded.xlsx", 
    "production" => production_expanded_df,
    "area_results" => area_expanded_df,
    "hydro_results" => hydro_expanded_df,
    "transmission" => transmission_expanded_df,
     overwrite=true)
end   

function write_weight_results(model_results, input_sets, input_parameters)
    hydro_weights_df = get_weight_hydro_result(model_results, input_sets, input_parameters)
    production_weights_df = get_weight_prod_results(model_results, input_sets, input_parameters)
    area_weights_df = get_weight_area_results(model_results, input_sets, input_parameters)
    XLSX.writetable("continuous_results/result_weights.xlsx", 
                    "production" => production_weights_df, 
                    "area_results" => area_weights_df, 
                    "hydro_results" => hydro_weights_df,
                    overwrite=true)
end

function write_aggregated_results(model_results, input_sets, input_parameters)
    hydro_aggregated_df = get_aggregated_hydro_results(model_results, input_sets)
    CSV.write("continuous_results/hydro_results_aggregated.csv", hydro_aggregated_df)

    production_aggregated_df = get_aggregated_prod_results(model_results, input_sets, input_parameters)
    CSV.write("continuous_results/production_aggregated.csv", production_aggregated_df)
    area_aggregated_df = get_aggregated_area_results(model_results, input_sets, input_parameters)
    CSV.write("continuous_results/area_results_aggregated.csv", area_aggregated_df)
    transmission_aggregated_df = get_aggregated_transmission_results(model_results, input_sets, input_parameters)
    CSV.write("continuous_results/transmission_results_aggregated.csv", transmission_aggregated_df)

    XLSX.writetable("continuous_results/results_aggregated.xlsx", 
                    "production" => production_aggregated_df,
                    "area_results" => area_aggregated_df,
                    "hydro_results" => hydro_aggregated_df,
                    "transmission" => transmission_aggregated_df,
                     overwrite=true)
end

function get_expanded_frequency_result(sampling_points, frequency_weights, frequency_d_weights, T)
    frequency_df = DataFrame(id=Int[], timestep=Int[], timestep_fractional=Float64[], frequency=Float64[], frequency_d=Float64[])
    frequency_weights_df = DataFrame(Array(frequency_weights), :auto)
    display(frequency_weights_df)
    frequency_array = get_converted_list(frequency_weights_df, sampling_points)

    frequency_d_weights_df = DataFrame(Array(frequency_d_weights), :auto)
    frequency_d_array = get_converted_list(frequency_d_weights_df, sampling_points)

    df = DataFrame(
        id = 1:length(frequency_array),
        timestep = vcat([1], [t for t in T for s in 1:sampling_points]),
        timestep_fractional= 0:(1/sampling_points):T[end],
        frequency = frequency_array,
        frequncy_d = frequency_d_array,
    )
    display(df)
    return df
end

function get_weight_hydro_result(model_results, input_sets, input_parameters)
    @unpack discharge, spill, bypass, total_inflow, total_outflow = model_results
    @unpack T, P_h, B, S = input_sets
    hydro_weights_df = DataFrame(id=Int[], areas_fk=Int[], plant_id=Int[], timestep=Int[], scenario = Int[], b=Int[], discharge=Float64[], spill=Float64[], bypass=Float64[], controlled_inflow=Float64[], controlled_outflow=Float64[])
    id_counter = 1
    areas_fk = 1
    for t in T
        for p in P_h
            for b in B
                for s in S
                    row = (id_counter, areas_fk, p, t, s, b, round(discharge[s, p, b, t], digits=2), 
                            round(spill[s, p, b, t], digits=2), round(bypass[s, p, b, t], digits=2),
                            round(total_inflow[s, p, b, t], digits=2), round(total_outflow[s, p, b, t], digits=2))
                    push!(hydro_weights_df, row)
                    id_counter += 1
                end
            end
        end
        areas_fk += 1
    end
    return hydro_weights_df
end

function get_aggregated_hydro_results(model_results, input_sets)
    @unpack volume_end = model_results
    @unpack T, P_h, S = input_sets
    expanded_df = get_expanded_hydro_results(model_results, input_sets, 1000, T[end])

    df_avg = combine(groupby(expanded_df, [:timestep, :plant_id, :scenario]), 
                    :discharge => mean => :discharge,
                    :bypass => mean => :bypass,
                    :spill => mean => :spill,
                    :controlled_inflow => mean => :controlled_inflow,
                    :controlled_outflow => mean => :controlled_outflow,
                    :volume => mean => :volume)
                    # :scenario => first => :scenario)
    # new_array =  [volume_end[s, p, t] for t in T for p in P_h for s in S]
    df_avg.volume_res = [volume_end[s, p, t] for t in T for p in P_h for s in S]
    return df_avg
end


function get_expanded_hydro_results(model_results, input_sets, sampling_points, output_timesteps=24)
    @unpack discharge, spill, bypass, total_inflow, total_outflow, volume, volume_end = model_results
    @unpack T, P_h, S = input_sets
    hydro_df = DataFrame(DataFrame(id=Int[], plant_id = Int[], timestep=Int[], timestep_fractional=Float64[], scenario=Int[], discharge=Float64[], spill=Float64[], bypass=Float64[], controlled_inflow=Float64[], controlled_outflow=Float64[], volume=Float64[]))
    id_counter = 1
    sp_ts_input = div(sampling_points, div(T[end], 24))
    sp_ts_output = div(sampling_points, div(output_timesteps, 24))
    n_elements = 24*sampling_points + 1

    for p in P_h
        for s in S
            plant_discharge_weights = discharge[s, p, :, :]
            plant_discharge_weights_df = DataFrame(Array(plant_discharge_weights), :auto)
            plant_discharge_array = get_converted_list(plant_discharge_weights_df, sp_ts_input)

            plant_spill_weights = spill[s, p, :, :]
            plant_spill_weights_df = DataFrame(Array(plant_spill_weights), :auto)
            plant_spill_array = get_converted_list(plant_spill_weights_df, sp_ts_input)

            plant_bypass_weights = bypass[s, p, :, :]
            plant_bypass_weights_df = DataFrame(Array(plant_bypass_weights), :auto)
            plant_bypass_array = get_converted_list(plant_bypass_weights_df, sp_ts_input)

            plant_total_inflow_weights = total_inflow[s, p, :, :]
            plant_total_inflow_weights_df = DataFrame(Array(plant_total_inflow_weights), :auto)
            plant_total_inflow_array = get_converted_list(plant_total_inflow_weights_df, sp_ts_input)

            plant_total_outflow_weights = total_outflow[s, p, :, :]
            plant_total_outflow_weights_df = DataFrame(Array(plant_total_outflow_weights), :auto)
            plant_total_outflow_array = get_converted_list(plant_total_outflow_weights_df, sp_ts_input)

            volume_weights = volume[s, p, :, :]
            volume_weights_df = DataFrame(Array(volume_weights), :auto)
            volume_array = get_converted_list(volume_weights_df, sp_ts_input, false) 
            volume_array = [volume_end[s, p, t-1] + volume_array[sp + sp_ts_input*(t-1)] for t in T for sp in 1:sp_ts_input]
            pushfirst!(volume_array, volume_end[s, p, 0])

            area_results_df = DataFrame(
                id=id_counter:(id_counter+n_elements-1),
                plant_id = fill(p, n_elements),
                scenario = fill(s, n_elements),
                timestep = vcat([1], [t  for t in 1:output_timesteps for s in 1:sp_ts_output]),
                timestep_fractional=0:(1/sampling_points):24,
                discharge =plant_discharge_array,
                spill = plant_spill_array,
                bypass = plant_bypass_array,
                controlled_inflow = plant_total_inflow_array,
                controlled_outflow = plant_total_outflow_array,
                volume = volume_array 
                )
            append!(hydro_df, area_results_df)
            id_counter += n_elements
        end
    end
    hydro_df.timestep_fractional .= round.(hydro_df.timestep_fractional, digits=4)
    return hydro_df
end

function get_weight_prod_results(model_results, input_sets, input_parameters)
    @unpack production, up_activation, down_activation = model_results
    @unpack plant_df = input_parameters
    @unpack T, P, B, S = input_sets
    production_weights_df = DataFrame(id=Int[], areas_fk=Int[], plant_id=Int[], timestep=Int[], b=Int[], scenario=Int[], area=Int[], fuel_type=String[], production=Float64[], up_activation=Float64[], down_activation=Float64[])
    id_counter = 1
    areas_fk = 1
    for t in T
        for p in P
            for b in B
                for s in S
                    row = (id_counter, areas_fk, p, t, b, s, plant_df[plant_df.plant_id .== p, :area][1], plant_df[plant_df.plant_id .== p, :fuel_type][1], round(production[s, p, b, t], digits=2), round(up_activation[s, p, b, t], digits=2), round(down_activation[s, p, b, t], digits=2))
                    push!(production_weights_df, row)
                    id_counter += 1
                end
            end
        end
        areas_fk += 1
    end
    CSV.write("continuous_results/production_weights.csv", production_weights_df)
    return production_weights_df
end

function get_expanded_prod_results(model_results, input_sets, input_parameters, sampling_points, output_timesteps=24)
    @unpack production, up_activation, down_activation, uc, startup = model_results
    @unpack plant_df = input_parameters
    @unpack T, P, B, S = input_sets
    prod_df = DataFrame(id=Int[], plant_id=Int[], scenario=Int[], timestep=Int[], timestep_fractional=Float64[], production=Float64[], up_activation=Float64[], down_activation=Float64[], area=Int[], fuel_type=String[], status=Int[], startup=Int[])
    # prod_df = DataFrame(id=Int[], plant_id=Int[], timestep=Int[], area=Int[], fuel_type=String[], production=Float64[], status=Int[], startup=Int[])
    id_counter=1
    sp_ts_input = div(sampling_points, div(T[end], 24))
    sp_ts_output = div(sampling_points, div(output_timesteps, 24))
    n_elements = 24*sampling_points + 1

    for p in P
        for s in S
            weights = production[s, p, :, :]
            weights_df = DataFrame(Array(weights), :auto)
            # weights_df = dense_array_to_df(weights)
            prod_array = get_converted_list(weights_df, sp_ts_input)

            up_activation_weights = up_activation[s, p, :, :]
            up_activation_weights_df = DataFrame(Array(up_activation_weights), :auto)
            up_activation_array = get_converted_list(up_activation_weights_df, sp_ts_input)

            down_activation_weights = down_activation[s, p, :, :]
            down_activation_weights_df = DataFrame(Array(down_activation_weights), :auto)
            down_activation_array = get_converted_list(down_activation_weights_df, sp_ts_input)

            plant_prod_df = DataFrame(id=id_counter:(id_counter+n_elements-1),
                                plant_id=fill(p, n_elements),
                                scenario=fill(s, n_elements),
                                timestep = vcat([1], [t  for t in 1:output_timesteps for s in 1:sp_ts_output]),
                                timestep_fractional=0:(1/sampling_points):24, 
                                production = round.(prod_array, digits=2),
                                up_activation = round.(up_activation_array, digits=2),
                                down_activation = round.(down_activation_array, digits =2),
                                area=fill(plant_df[plant_df.plant_id .== p, :area][1], n_elements),
                                fuel_type=fill(plant_df[plant_df.plant_id .== p, :fuel_type][1], n_elements),
                                status=vcat([round(uc[p, 1])], [round(uc[p, t]) for t in 1:T[end] for s in 1:sp_ts_input]),
                                startup=vcat([round(startup[p, 1])], [round(startup[p, t]) for t in 1:T[end] for s in 1:sp_ts_input]))
            append!(prod_df, plant_prod_df)
            id_counter += n_elements
        end
    end
    prod_df.timestep_fractional .= round.(prod_df.timestep_fractional, digits=4)
    return prod_df
end

function get_aggregated_prod_results(model_results, input_sets, input_parameters)
    @unpack T = input_sets
    expanded_df = get_expanded_prod_results(model_results, input_sets, input_parameters, 1000, T[end])
    df_avg = combine(groupby(expanded_df, [:timestep, :plant_id, :scenario]), 
                    :production => mean => :production,
                    :up_activation => mean => :up_activation,
                    :down_activation => mean => :down_activation,
                    :area => first => :area,
                    :fuel_type => first => :fuel_type,
                    :status => first => :status,
                    :startup => first => :startup)
    # new_array =  [volume_end[p, t] for t in T for p in P_h]
    # df_avg.volume_res = [volume_end[p, t] for t in T for p in P_h]
    df_avg.production = round.(df_avg.production, digits=2)
    return df_avg
end

function get_weight_area_results(model_results, input_sets, input_parameters)
    @unpack shedding_weights, dumping_weights = model_results
    @unpack T, A, B, S = input_sets
    @unpack load_df = input_parameters

    area_weights_df = DataFrame(id=Int[], area=Int[], timestep=Int[], b=Int[], scenario=Int[], load=Float64[], shedding=Float64[], dumping=Float64[])
    id_counter = 1
    for a in A
        for t in T
            for b in B
                for s in S
                    row = (id_counter, a, t, b, s, load_df[(load_df.area .== a) .& (load_df.timestep .== t) .& (load_df.scenario .== s), :load][1], round(shedding_weights[s, b, a, t], digits=2), round(dumping_weights[s, b, a, t], digits=2))
                    push!(area_weights_df, row)
                    id_counter += 1
                end
            end
        end
    end
    CSV.write("continuous_results/area_result_weights.csv", area_weights_df)
    return area_weights_df
end

function get_expanded_area_results(model_results, input_sets, input_parameters, sampling_points, output_timesteps=24)
    @unpack shedding_weights, dumping_weights = model_results
    @unpack load_df = input_parameters
    @unpack A, T, S, = input_sets
    area_df = DataFrame(id=Int[], area=Int[], scenario=Int[], timestep=Int64[], timestep_fractional=Float64[], load=Float64[], load_shedding=Float64[], power_dumping=Float64[])
    id_counter = 1
    sp_ts_input = div(sampling_points, div(T[end], 24))
    sp_ts_output = div(sampling_points, div(output_timesteps, 24))
    n_elements = 24*sampling_points + 1

    for a in A
        for s in S
            area_shedding_weights = shedding_weights[s, :, a, :]
            area_shedding_weights_df =DataFrame(Array(area_shedding_weights), :auto)
            area_shedding_array = get_converted_list(area_shedding_weights_df, sp_ts_input)

            area_dumping_weights = dumping_weights[s, :, a, :]
            area_dumping_weights_df  = DataFrame(Array(area_dumping_weights), :auto)
            area_dumping_array = get_converted_list(area_dumping_weights_df, sp_ts_input)
            
            area_load_weights_df = load_df[(load_df.area .== a) .& (load_df.scenario .== s), :]
            area_load_weights_df = select!(unstack(area_load_weights_df, :b, :timestep, :load), Not(:b))
            area_load_array = get_converted_list(area_load_weights_df, sp_ts_input)
            
            area_results_df = DataFrame(
                id=id_counter:(id_counter+n_elements-1),
                area = fill(a, n_elements),
                scenario = fill(s, n_elements),
                timestep = vcat([1], [t for t in 1:output_timesteps for s in 1:sp_ts_output]),
                timestep_fractional=0:(1/sampling_points):24, 
                load = area_load_array,
                load_shedding = area_shedding_array,
                power_dumping = area_dumping_array, 
                )
            append!(area_df, area_results_df)
            id_counter += n_elements
        end
    end
    area_df.timestep_fractional .= round.(area_df.timestep_fractional, digits=4)
    return area_df
end

function get_aggregated_area_results(model_results, input_sets, input_parameters)
    @unpack T = input_sets
    expanded_df = get_expanded_area_results(model_results, input_sets, input_parameters, 1000, T[end])
    df_avg = combine(groupby(expanded_df, [:timestep, :area, :scenario]),
                    :load => mean => :load,
                    :load_shedding => mean => :load_shedding,
                    :power_dumping => mean => :power_dumping)
    return df_avg
end

function get_expanded_transmission_results(model_results, input_sets, input_parameters, sampling_points, output_timesteps=24)
    @unpack transmission = model_results
    @unpack line_df = input_parameters
    @unpack T, L, S = input_sets
    transmission_df = DataFrame(id=Int[], scenario=Int[], line_id=Int[], timestep=Int[], timestep_fractional=Float64[], area_from=Int[], area_to=Int[], transmission=Float64[])
    id_counter = 1
    sp_ts_input = div(sampling_points, div(T[end], 24))
    sp_ts_output = div(sampling_points, div(output_timesteps, 24))
    n_elements = 24*sampling_points + 1
    for l in L
        for s in S
            line_transmission_weights = transmission[s, l, :, :]
            line_transmission_weights_df = DataFrame(Array(line_transmission_weights), :auto)
            line_transmission_array = get_converted_list(line_transmission_weights_df, sp_ts_input)

            line_transmission_df = DataFrame(
                id=id_counter:(id_counter+n_elements-1),
                line_id = fill(l, n_elements),
                scenario = fill(s, n_elements),
                timestep = vcat([1], [t  for t in 1:output_timesteps for s in 1:sp_ts_output]),
                timestep_fractional=0:(1/sampling_points):24,
                area_from = fill(line_df[line_df.line_id .== l, :area_from][1], n_elements),
                area_to = fill(line_df[line_df.line_id.== l, :area_to][1], n_elements),
                transmission = round.(line_transmission_array, digits=2)
            )
            # display(line_transmission_df)
            append!(transmission_df, line_transmission_df)
            id_counter += n_elements
        end
    end
    transmission_df.timestep_fractional .= round.(transmission_df.timestep_fractional, digits=4)
    return transmission_df
end

function get_aggregated_transmission_results(model_results, input_sets, input_parameters)
    @unpack T = input_sets
    expanded_df = get_expanded_transmission_results(model_results, input_sets, input_parameters, 10000, T[end])
    df_avg = combine(groupby(expanded_df, [:timestep, :line_id, :scenario]),
                    :transmission => mean => :transmission,
                    :area_from => first => :area_from,
                    :area_to => first => :area_to)
    df_avg.transmission = round.(df_avg.transmission, digits=2)
    return df_avg
end