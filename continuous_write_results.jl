using Statistics

function write_results(model)
    wind_df, hydro_df, inflow_df, load_df, plant_df, line_df, A, P, T, L, B, P_w, P_t, P_h, P_a, L_in, L_out, I_disch, I_spill, I_bypass, C_shedding, C_dumping, C_startup = read_input_data()
    sampling_points = 10
    production = value.(model[:prod])
    uc = value.(model[:status])
    startup = value.(model[:startup])
    shedding = value.(model[:shedding])
    dumping = value.(model[:dumping])
    transmission = value.(model[:transmission])

    discharge = value.(model[:flow_disch])
    spill = value.(model[:flow_spill])
    bypass = value.(model[:flow_bypass])
    total_inflow = value.(model[:total_flow_in])
    total_outflow = value.(model[:total_flow_out])
    volume_end = value.(model[:volume_end])
    frequency = value.(model[:frequency])
    frequency_d = value.(model[:frequency_d])

    hydro_aggregated_df = get_aggregated_hydro_results(discharge, spill, bypass, total_inflow, total_outflow, volume_end, T, P_h)
    hydro_weights_df = get_weight_hydro_result(discharge, spill, bypass, total_inflow, total_outflow, T, P_h, B)
    hydro_expanded_df = get_expanded_hydro_results(discharge, spill, bypass, total_inflow, total_outflow, T, P_h, sampling_points)
    
    CSV.write("continuous_results/hydro_results_expanded.csv", hydro_expanded_df)
    CSV.write("continuous_results/hydro_results_aggregated.csv", hydro_aggregated_df)
    
    production_weights_df = get_weight_prod_results(plant_df, production, T, P, B)
    production_expanded_df = get_expanded_prod_results(plant_df, production, uc, startup, T, P, B, sampling_points)
    production_aggregated_df = get_aggregated_prod_results(plant_df, production, uc, startup, T, P, B,)
    CSV.write("continuous_results/production_expanded.csv", production_expanded_df)
    CSV.write("continuous_results/production_aggregated.csv", production_aggregated_df)

    area_weights_df = get_weight_area_results(load_df, shedding, dumping, T, A, B)
    area_expanded_df = get_expanded_area_results(load_df, shedding, dumping, A, T, sampling_points)
    area_aggregated_df = get_aggregated_area_results(load_df, shedding, dumping, A, T)
    CSV.write("continuous_results/area_results_expanded.csv", area_expanded_df)
    CSV.write("continuous_results/area_results_aggregated.csv", area_aggregated_df)

    transmission_expanded_df = get_expanded_transmission_results(line_df, transmission, T, L, sampling_points)
    transmission_aggregated_df = get_aggregated_transmission_results(line_df, transmission, T, L)
    CSV.write("continuous_results/transmission_results_expanded.csv", transmission_expanded_df)
    CSV.write("continuous_results/transmission_results_aggregated.csv", transmission_aggregated_df)

    frequency_expanded_df = get_expanded_frequency_result(sampling_points, frequency, frequency_d, T)
    CSV.write("continuous_results/frequency.csv", frequency_expanded_df)

    XLSX.writetable("continuous_results/results_expanded.xlsx", 
                    "production" => production_expanded_df,
                    "area_results" => area_expanded_df,
                    "hydro_results" => hydro_expanded_df,
                    "transmission" => transmission_expanded_df, overwrite=true)
    XLSX.writetable("continuous_results/results_aggregated.xlsx", 
                    "production" => production_aggregated_df,
                    "area_results" => area_aggregated_df,
                    "hydro_results" => hydro_aggregated_df,
                    "transmission" => transmission_aggregated_df, overwrite=true)
    XLSX.writetable("continuous_results/result_weights.xlsx", 
                    "production" => production_weights_df, 
                    "area_results" => area_weights_df, 
                    "hydro_results" => hydro_weights_df, overwrite=true)

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

function get_weight_hydro_result(discharge, spill, bypass, total_inflow, total_outflow, T, P_h, B)
    hydro_weights_df = DataFrame(id=Int[], areas_fk=Int[], plant_id=Int[], timestep=Int[], b=Int[], discharge=Float64[], spill=Float64[], bypass=Float64[], controlled_inflow=Float64[], controlled_outflow=Float64[])
    id_counter = 1
    areas_fk = 1
    for t in T
        for p in P_h
            for b in B
                row = (id_counter, areas_fk, p, t, b, round(discharge[p, b, t], digits=2), round(spill[p, b, t], digits=2), round(bypass[p, b, t], digits=2), round(total_inflow[p, b, t], digits=2), round(total_outflow[p, b, t], digits=2))
                push!(hydro_weights_df, row)
                id_counter += 1
            end
        end
        areas_fk += 1
    end
    return hydro_weights_df
end

function get_aggregated_hydro_results(discharge, spill, bypass, total_inflow, total_outflow, volume_end, T, P_h)
    expanded_df = get_expanded_hydro_results(discharge, spill, bypass, total_inflow, total_outflow, T, P_h, 10000)

    df_avg = combine(groupby(expanded_df, [:timestep, :plant_id]), 
                    :discharge => mean => :discharge,
                    :bypass => mean => :bypass,
                    :spill => mean => :spill,
                    :controlled_inflow => mean => :controlled_inflow,
                    :controlled_outflow => mean => :controlled_outflow)
    new_array =  [volume_end[p, t] for t in T for p in P_h]
    df_avg.volume_res = [volume_end[p, t] for t in T for p in P_h]
    return df_avg
end


function get_expanded_hydro_results(discharge, spill, bypass, total_inflow, total_outflow, T, P_h, sampling_points)
    hydro_df = DataFrame(DataFrame(id=Int[], plant_id = Int[], timestep=Int[], timestep_fractional=Float64[], discharge=Float64[], spill=Float64[], bypass=Float64[], controlled_inflow=Float64[], controlled_outflow=Float64[]))
    id_counter = 1
    n_elements = T[end]*sampling_points + 1

    for p in P_h
        plant_discharge_weights = discharge[p, :, :]
        plant_discharge_weights_df = DataFrame(Array(plant_discharge_weights), :auto)
        plant_discharge_array = get_converted_list(plant_discharge_weights_df, sampling_points)

        plant_spill_weights = spill[p, :, :]
        plant_spill_weights_df = DataFrame(Array(plant_spill_weights), :auto)
        plant_spill_array = get_converted_list(plant_spill_weights_df, sampling_points)

        plant_bypass_weights = bypass[p, :, :]
        plant_bypass_weights_df = DataFrame(Array(plant_bypass_weights), :auto)
        plant_bypass_array = get_converted_list(plant_bypass_weights_df, sampling_points)

        plant_total_inflow_weights = total_inflow[p, :, :]
        plant_total_inflow_weights_df = DataFrame(Array(plant_total_inflow_weights), :auto)
        plant_total_inflow_array = get_converted_list(plant_total_inflow_weights_df, sampling_points)

        plant_total_outflow_weights = total_outflow[p, :, :]
        plant_total_outflow_weights_df = DataFrame(Array(plant_total_outflow_weights), :auto)
        plant_total_outflow_array = get_converted_list(plant_total_outflow_weights_df, sampling_points)

        area_results_df = DataFrame(
            id=id_counter:(id_counter+n_elements-1),
            plant_id = fill(p, n_elements),
            timestep = vcat([1], [t  for t in T for s in 1:sampling_points]),
            timestep_fractional=0:(1/sampling_points):T[end],
            discharge =plant_discharge_array,
            spill = plant_spill_array,
            bypass = plant_bypass_array,
            controlled_inflow = plant_total_inflow_array,
            controlled_outflow = plant_total_outflow_array 
            )
        append!(hydro_df, area_results_df)
        id_counter += n_elements
    end
    CSV.write("continuous_results/hydro_results_expanded.csv", hydro_df)
    return hydro_df
end

function get_weight_prod_results(plant_df, production_results, T, P, B)
    production_weights_df = DataFrame(id=Int[], areas_fk=Int[], plant_id=Int[], timestep=Int[], b=Int[], area=Int[], fuel_type=String[], production=Float64[])
    id_counter = 1
    areas_fk = 1
    for t in T
        for p in P
            for b in B
                row = (id_counter, areas_fk, p, t, b, plant_df[plant_df.plant_id .== p, :area][1], plant_df[plant_df.plant_id .== p, :fuel_type][1], round(production_results[p, b, t], digits=2))
                push!(production_weights_df, row)
                id_counter += 1
            end
        end
        areas_fk += 1
    end
    CSV.write("continuous_results/production_weights.csv", production_weights_df)
    return production_weights_df
end

function get_expanded_prod_results(plant_df, production_results, uc, startup, T, P, B, sampling_points)
    prod_df = DataFrame(id=Int[], plant_id=Int[], timestep=Int[], timestep_fractional=Float64[], production=Float64[], area=Int[], fuel_type=String[], status=Int[], startup=Int[])
    # prod_df = DataFrame(id=Int[], plant_id=Int[], timestep=Int[], area=Int[], fuel_type=String[], production=Float64[], status=Int[], startup=Int[])
    id_counter=1
    n_elements = T[end]*sampling_points + 1
    for p in P
        weights = production_results[p, :, :]
        weights_df = DataFrame(Array(weights), :auto)
        # weights_df = dense_array_to_df(weights)
        prod_array = get_converted_list(weights_df, sampling_points)
        plant_prod_df = DataFrame(id=id_counter:(id_counter+n_elements-1),
                             plant_id=fill(p, n_elements),
                             timestep = vcat([1], [t  for t in T for s in 1:sampling_points]),
                             timestep_fractional=0:(1/sampling_points):T[end], 
                             production=prod_array,
                             area=fill(plant_df[plant_df.plant_id .== p, :area][1], n_elements),
                             fuel_type=fill(plant_df[plant_df.plant_id .== p, :fuel_type][1], n_elements),
                             status=vcat([round(uc[p, 1])], [round(uc[p, t]) for t in T for s in 1:sampling_points]),
                             startup=vcat([round(startup[p, 1])], [round(startup[p, t]) for t in T for s in 1:sampling_points]))
        append!(prod_df, plant_prod_df)
        id_counter += n_elements
    end
    return prod_df
end

function get_aggregated_prod_results(plant_df, production_results, uc, startup, T, P, B,)
    expanded_df = get_expanded_prod_results(plant_df, production_results, uc, startup, T, P, B, 10000)
    df_avg = combine(groupby(expanded_df, [:timestep, :plant_id]), 
                    :production => mean => :production,
                    :area => first => :area,
                    :fuel_type => first => :fuel_type,
                    :status => first => :status,
                    :startup => first => :startup)
    # new_array =  [volume_end[p, t] for t in T for p in P_h]
    # df_avg.volume_res = [volume_end[p, t] for t in T for p in P_h]
    return df_avg
end

function get_weight_area_results(load_df, shedding, dumping, T, A, B)
    area_weights_df = DataFrame(id=Int[], area=Int[], timestep=Int[], b=Int[], load=Float64[], shedding=Float64[], dumping=Float64[])
    id_counter = 1
    for a in A
        for t in T
            for b in B
                row = (id_counter, a, t, b, load_df[(load_df.area .== a) .& (load_df.timestep .== t), :Forbruk][1], round(shedding[b, a, t], digits=2), round(dumping[b, a, t], digits=2))
                push!(area_weights_df, row)
                id_counter += 1
            end
        end
    end
    CSV.write("continuous_results/area_result_weights.csv", area_weights_df)
    return area_weights_df
end

function get_expanded_area_results(load_df, shedding_weights, dumping_weights, A, T, sampling_points)
    area_df = DataFrame(id=Int[], area=Int[],timestep=Int64[], timestep_fractional=Float64[], load=Float64[], load_shedding=Float64[], power_dumping=Float64[])
    id_counter = 1
    n_elements = T[end]*sampling_points + 1

    for a in A
        area_shedding_weights = shedding_weights[:, a, :]
        area_shedding_weights_df =DataFrame(Array(area_shedding_weights), :auto)
        area_shedding_array = get_converted_list(area_shedding_weights_df, sampling_points)

        area_dumping_weights = dumping_weights[:, a, :]
        area_dumping_weights_df  = DataFrame(Array(area_dumping_weights), :auto)
        area_dumping_array = get_converted_list(area_dumping_weights_df, sampling_points)
        
        area_load_weights_df = load_df[load_df.area .== a, :]
        area_load_weights_df = select!(unstack(area_load_weights_df, :b, :timestep, :Forbruk), Not(:b))
        area_load_array = get_converted_list(area_load_weights_df, sampling_points)
        
        area_results_df = DataFrame(
            id=id_counter:(id_counter+n_elements-1),
            area = fill(a, n_elements),
            timestep = vcat([1], [t  for t in T for s in 1:sampling_points]),
            timestep_fractional=0:(1/sampling_points):T[end], 
            load = area_load_array,
            load_shedding = area_shedding_array,
            power_dumping = area_dumping_array, 
            )
        append!(area_df, area_results_df)
        id_counter += n_elements
    end
    return area_df
end

function get_aggregated_area_results(load_df, shedding_weights, dumping_weights, A, T)
    expanded_df = get_expanded_area_results(load_df, shedding_weights, dumping_weights, A, T, 10000)
    df_avg = combine(groupby(expanded_df, [:timestep, :area]),
                    :load => mean => :load,
                    :load_shedding => mean => :load_shedding,
                    :power_dumping => mean => :power_dumping)
    return df_avg
end

function get_expanded_transmission_results(line_df, transmission_weights, T, L, sampling_points)
    transmission_df = DataFrame(id=Int[], line_id=Int[], timestep=Int[], timestep_fractional=Float64[], area_from=Int[], area_to=Int[], transmission=Float64[])
    id_counter = 1
    n_elements = T[end]*sampling_points + 1
    for l in L
        line_transmission_weights = transmission_weights[l, :, :]
        line_transmission_weights_df = DataFrame(Array(line_transmission_weights), :auto)
        line_transmission_array = get_converted_list(line_transmission_weights_df, sampling_points)

        line_transmission_df = DataFrame(
            id=id_counter:(id_counter+n_elements-1),
            line_id = fill(l, n_elements),
            timestep = vcat([1], [t  for t in T for s in 1:sampling_points]),
            timestep_fractional=0:(1/sampling_points):T[end],
            area_from = fill(line_df[line_df.line_id .== l, :area_from][1], n_elements),
            area_to = fill(line_df[line_df.line_id.== l, :area_to][1], n_elements),
            transmission = line_transmission_array
        )
        # display(line_transmission_df)
        append!(transmission_df, line_transmission_df)
        id_counter += n_elements
    end
    return transmission_df
end

function get_aggregated_transmission_results(line_df, transmission_weights, T, L)
    expanded_df = get_expanded_transmission_results(line_df, transmission_weights, T, L, 10000)
    df_avg = combine(groupby(expanded_df, [:timestep, :line_id]),
                    :transmission => mean => :transmission,
                    :area_from => first => :area_from,
                    :area_to => first => :area_to)
    return df_avg
end