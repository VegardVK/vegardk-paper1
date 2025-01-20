using XLSX
using DataFrames
using JuMP
using CPLEX
using UnPack

include("C:/Users/vegardvk/vscodeProjects/bernstein/get_hydro_data.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/helper_functions.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/find_bernstein_weights.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/plot_results.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/elevation_matrix.jl")



function read_input_data(S_input)
    # Read input data
    input_parameters = (
        wind_df = DataFrame(XLSX.readtable("output/wind_power_weights.xlsx", "Sheet1", infer_eltypes=true)),
        hydro_df = DataFrame(XLSX.readtable("output/hydro_data.xlsx", "Sheet1", infer_eltypes=true)),
        inflow_df = DataFrame(XLSX.readtable("output/inflow_weights.xlsx", "Sheet1", infer_eltypes=true)),
        load_df = DataFrame(XLSX.readtable("output/load_weights.xlsx", "Sheet1", infer_eltypes=true)),
        prod_df = DataFrame(XLSX.readtable("output/production_weights.xlsx", "Sheet1", infer_eltypes=true)),
        shedding_df = DataFrame(XLSX.readtable("output/shedding_weights.xlsx", "Sheet1", infer_eltypes=true)),
        dumping_df = DataFrame(XLSX.readtable("output/dumping_weights.xlsx", "Sheet1", infer_eltypes=true)),
        plant_df = DataFrame(XLSX.readtable("output/plant_data.xlsx", "Sheet1", infer_eltypes=true)),
        line_df = DataFrame(XLSX.readtable("output/line_data.xlsx", "Sheet1", infer_eltypes=true)),
    )
    @unpack plant_df, load_df, line_df = input_parameters
    A = unique(plant_df.area)
    area_grouped_plants = groupby(plant_df, :area)

    if !(length(S_input)>=0)
        S = unique(load_df.s)
        S_input = S[2:end]
    end
    input_sets = (
        A = unique(plant_df.area),
        P = unique(plant_df.plant_id),
        T = unique(load_df.timestep),
        L = unique(line_df.line_id),
        B = unique(load_df.b),
        S = S_input,
        P_w = unique(plant_df[plant_df.fuel_type .== "Wind", :plant_id]),
        P_t = unique(plant_df[plant_df.fuel_type .== "Thermal", :plant_id]),
        P_h = unique(plant_df[plant_df.fuel_type .== "Hydro", :plant_id]),
        P_a = Dict(g.area[1] => unique(g.plant_id) for g in area_grouped_plants),
        # Dicts with the lines that go in/out of area a
        L_in = Dict(a => unique(line_df[line_df.area_to .== a, :line_id]) for a in A),
        L_out =  Dict(a => unique(line_df[line_df.area_from .== a, :line_id]) for a in A),
        line_pairs = Dict(pair => filter(row -> row.line_pair == pair, eachrow(line_df)).line_id for pair in unique(line_df.line_pair)),
        L_P = unique(line_df.line_pair),
    )

    return input_parameters, input_sets
end

function define_and_solve_model(S_input=[])
    input_parameters, input_sets = read_input_data(S_input)
    @unpack wind_df, hydro_df, inflow_df, load_df, plant_df, line_df, prod_df, shedding_df, dumping_df = input_parameters
    @unpack A, P, T, L, B, P_w, P_t, P_h, P_a, L_in, L_out, S, L_P, line_pairs = input_sets
    I_disch, I_spill, I_bypass = find_connected_plants(hydro_df)
    C_shedding, C_dumping, C_startup = get_cost_parameters()

    Δt = 24/T[end]
    π = 1/length(S)
    B_I = 0:(B[end]+1)
    B_D = 0:(B[end]-1)
    elevation_matrix = get_elevation_matrix(B[end]-1, B[end])
    k_matrix = get_k_matrix(B[end])

    T_extended = 1:(T[end] + 1)
    model = Model()
    set_optimizer(model, CPLEX.Optimizer)

    @variable(model, prod[s in S, p in P, b in B, t in T] ≥ 0)
    @variable(model, up_activation[s in S, p in P, b in B, t in T] ≥ 0)
    @variable(model, down_activation[s in S, p in P, b in B, t in T] ≥ 0)

    @variable(model, transmission[s in S, l in L, b in B, t in T] ≥ 0)

    @variable(model, 0 ≤ shedding[s in S, b in B, a in A, t in T] ≤ load_df[(load_df.b .== b) .& (load_df.area .== a) .& (load_df.timestep .== t) .& (load_df.scenario .== s), :load][1]) 
    @variable(model, 0 ≤ dumping[s in S, b in B, a in A, t in T])
    @constraint(model, ub_dumping[s in S, b in B, a in A, t in T], dumping[s, b, a, t] ≤ sum(prod[s, p, b, t] for p in P_a[a]))
    
    @variable(model, flow_disch[s in S, p in P_h, b in B, t in T] ≥ 0) # Antar at alle moduler bare har ett discharge-segment
    @variable(model, flow_bypass[s in S, p in P_h, b in B, t in T] ≥ 0)
    @variable(model, flow_spill[s in S, p in P_h, b in B, t in T] ≥ 0)
    @variable(model, total_flow_in[s in S, p in P_h, b in B, t in T] ≥ 0)
    @variable(model, total_flow_out[s in S, p in P_h, b in B, t in T] ≥ 0)
    @variable(model, net_inflow[s in S, p in P_h, b in B, t in T])
    @variable(model, volume_end[s in S, p in P_h, t in 0:T[end]] ≥ 0)
    @variable(model, volume[s in S, p in P_h, b in B_I, t in T])

    @constraint(model, wind_production[s in S, p in P_w, b in B, t in T], prod[s, p, b, t] ≤ wind_df[(wind_df.plant_id .== p) .& (wind_df.timestep .== t) .&
                                                                                                        (wind_df.scenario .== s) .& (wind_df.b .== b), :wind_power][1])

    @constraint(model, controlled_inflow[s in S, p in P_h, b in B, t in T], total_flow_in[s, p, b, t] == sum(flow_disch[s, i, b, t] for i in I_disch[p]) 
                                                                                                            + sum(flow_bypass[s, i, b, t] for i in I_bypass[p]) 
                                                                                                            + sum(flow_spill[s, i, b, t] for i in I_spill[p]))

    @constraint(model, controlled_outflow[s in S, p in P_h, b in B, t in T], total_flow_out[s, p, b, t] == flow_disch[s, p, b, t] + flow_bypass[s, p, b, t] + flow_spill[s, p, b, t])

    @constraint(model, hydro_production[s in S, p in P_h, b in B, t in T], prod[s, p, b, t] == hydro_df[hydro_df.plant_id .== p, :enekv][1] * flow_disch[s, p, b, t] * 3.6)

    @constraint(model, bypass_ub[s in S, p in P_h, b in B, t in T], flow_bypass[s, p, b, t] ≤ hydro_df[hydro_df.plant_id .== p, :kap_forb][1])
    # test_set = [49908, 49915, 49930, 49929]
    # @constraint(model, prod_ub[s in S, p in test_set, b in B, t in T], flow_disch[s, p, b, t] ≤ hydro_df[hydro_df.plant_id .== p, :kap_gen_m3s][1])
    @constraint(model, prod_ub[s in S, p in P_h, b in B, t in T], flow_disch[s, p, b, t] ≤ hydro_df[hydro_df.plant_id .== p, :kap_gen_m3s][1])

    @constraint(model, spill_ub[s in S, p in P_h, b in B, t in T], flow_spill[s, p, b, t] ≤ hydro_df[hydro_df.plant_id .== p, :kap_spill][1])

    @constraint(model, net_inflow_calculation[s in S, p in P_h, b in B, t in T], net_inflow[s, p, b, t] == total_flow_in[s, p, b, t] 
                                                                                                + inflow_df[(inflow_df.plant_id .== p) .& (inflow_df.timestep .== t) .& (inflow_df.b .== b), :inflow][1] 
                                                                                                - total_flow_out[s, p, b, t])

    @constraint(model, starting_reservoir[s in S, p in P_h], volume_end[s, p, 0] == hydro_df[hydro_df.plant_id .== p, :starting_reservoir][1])
    @constraint(model, volume_calculation[s in S, p in P_h, t in T], volume_end[s, p, t] - volume_end[s, p, t-1] == Δt * (1/(B[end]+1)) * sum(net_inflow[s, p, b, t] for b in B))
    @constraint(model, simple_vol_ub[s in S, p in P_h, t in T], volume_end[s, p, t] ≤ hydro_df[hydro_df.plant_id .== p, :kap_mag][1])

    @constraint(model, volume_weights[s in S, p in P_h, b in B_I, t in T], volume[s, p, b, t] == Δt * (1/(B[end]+1)) * sum(net_inflow[s, p, b2, t] for b2 in B[1:b]))
    # @constraint(model, volume_ub[s in S, p in P_h, b in B_I, t in T], volume_end[s, p, t-1] + volume[s, p, b, t] ≤ hydro_df[hydro_df.plant_id .== p, :kap_mag][1])

    @variable(model, status[p in P, t in T_extended], Bin)
    @variable(model, startup[p in P, t in T_extended], Bin)
    @constraint(model, startup_count[p in P_t, t in T_extended[2:end]], startup[p, t] ≥ status[p, t] - status[p, t-1])
    # @constraint(model, min_production[s in S, p in union(P_t, P_h), t in T, b in B], prod[s, p, b, t] >=  
    #                                                                                         prod_df[(prod_df.scenario .== s) .& (prod_df.plant_id .== p) .& (prod_df.timestep .== t) .& (prod_df.b .== b), :lb][1]) 
    @constraint(model, min_production[s in S, p in union(P_t, P_h), t in T, b in B], prod[s, p, b, t] ≥  
                prod_df[(prod_df.scenario .== s) .& (prod_df.plant_id .== p) .& (prod_df.timestep .== t) .& (prod_df.b .== b), :lb][1])                                                                                         
    @constraint(model, max_production[s in S, p in union(P_t, P_h), t in T, b in B], prod[s, p, b, t] ≤  
                prod_df[(prod_df.scenario .== s) .& (prod_df.plant_id .== p) .& (prod_df.timestep .== t) .& (prod_df.b .== b), :ub][1])
    # display(prod_df)
    @constraint(model, balancing_activation[s in S, p in union(P_t, P_h), t in T, b in B], prod[s, p, b, t] ==
                prod_df[(prod_df.scenario .== 0) .& (prod_df.plant_id .== p) .& (prod_df.timestep .== t) .& (prod_df.b .== b), :production][1]
                + up_activation[s, p, b, t] - down_activation[s, p, b, t])

    # @constraint(model, unit_commitment1a[s in S, p in P_t, t in T], prod[s, p, 0, t] ≤ plant_df[plant_df.plant_id .== p, :gen_ub][1] * status[p, t])
    # @constraint(model, unit_commitment1b[s in S, p in P_t, t in T], prod[s, p, 0, t] ≥ plant_df[plant_df.plant_id .== p, :gen_lb][1] * status[p, t])

    # @constraint(model, unit_commitment2a[s in S, p in P_t, t in T], prod[s, p, B[end], t] ≤ plant_df[plant_df.plant_id .== p, :gen_ub][1] * status[p, t+1])
    # @constraint(model, unit_commitment2b[s in S, p in P_t, t in T], prod[s, p, B[end], t] ≥ plant_df[plant_df.plant_id .== p, :gen_lb][1] * status[p, t+1])

    # @constraint(model, unit_commitment3a[s in S, p in P_t, t in T], prod[s, p, 1, t] ≤ plant_df[plant_df.plant_id .== p, :gen_ub][1] * status[p, t])
    # @constraint(model, unit_commitment3b[s in S, p in P_t, t in T], prod[s, p, 1, t] ≥ plant_df[plant_df.plant_id .== p, :gen_lb][1] * status[p, t])

    # @constraint(model, unit_commitment4a[s in S, p in P_t, t in T], prod[s, p, B[end]-1, t] ≤ plant_df[plant_df.plant_id .== p, :gen_ub][1] * status[p, t+1])
    # @constraint(model, unit_commitment4b[s in S, p in P_t, t in T], prod[s, p, B[end]-1, t] ≥ plant_df[plant_df.plant_id .== p, :gen_lb][1] * status[p, t+1])

    @constraint(model, continuity_constraint1[s in S, p in P_t, t in T[1:end-1]], prod[s, p, B[end], t] == prod[s, p, 0, t+1])
    @constraint(model, continuity_constraint2[s in S, p in P_t, t in T[1:end-1]], prod[s, p, B[end], t] - prod[s, p, B[end]-1, t] == prod[s, p, 1, t+1] - prod[s, p, 0, t+1])


    # @variable(model, frequency[i in B, t in T])
    # @variable(model, frequency_pos[i in B, t in T] ≥ 0)
    # @variable(model, frequency_neg[i in B, t in T] ≤ 0)
    # @variable(model, frequency_d[i in B_D, t in T])

    # @constraint(model, frequency_differentiated[i in B_D, t in T], frequency_d[i, t] .== sum(frequency[j, t] * k_matrix[j+1, i+1] for j in B))
    # @constraint(model, frequency_continuous[t in T[1:end-1]], frequency[B[end], t] == frequency[0, t+1])
    # @constraint(model, freqency_ub[b in B, t in T], frequency[b, t] ≤ 2)
    # @constraint(model, frequency_lb[b in B, t in T], frequency[b, t] ≥ -2)

    # @constraint(model, frequency_pos_lb[b in B, t in T], frequency_pos[b, t] ≥ frequency[b, t])
    # @constraint(model, frequency_neg_ub[b in B, t in T], frequency_neg[b, t] ≤ frequency[b , t])

    # D_eff = 6
    # M_eff = 6

    @constraint(model, energy_balance[s in S, b in B, t in T, a in A], 
                # - M_eff * frequency[b, t]
                # - D_eff * sum(frequency_d[j, t] * elevation_matrix[j+1, b+1] for j in B_D)
                + sum(prod[s, p, b, t] for p in P_a[a])
                + 0.99 * sum(transmission[s, l, b, t] for l in L_in[a])
                - sum(transmission[s, l, b, t] for l in L_out[a])
                + shedding[s, b, a, t]
                - dumping[s, b, a, t]
                == load_df[(load_df.scenario .== s) .& (load_df.area .== a) .& (load_df.timestep .== t) .& (load_df.b .== b), :load][1])


    @variable(model, line_active[s in S, l in L, t in T], Bin)
    @constraint(model, one_directional_flow[s in S, pair in L_P, t in T], line_active[s, line_pairs[pair][1], t] + line_active[s, line_pairs[pair][2], t] .== 1)
    @constraint(model, transmission_ub[s in S, l in L, b in B, t in T], transmission[s, l, b, t] ≤ line_active[s, l, t] * line_df[line_df.line_id .== l, :capacity][1])

    # @constraint(model, max_transmission[s in S, l in L, b in B, t in T], 0 ≤ transmission[s, l, b, t] ≤  line_df[line_df.line_id .== l, :capacity][1])

    # @objective(model, Min, sum(prod[p, b, t] * get_bernstein_val(bernstein_degree, b, s/sampling_points) * p_dict["C"][p] for p in P for b in B for t in T for s in 0:sampling_points))
    bypass_penalty = 100
    spill_penalty = 100
    @objective(model, Min, 
                1/(B[end]+1) * Δt * π * (sum(prod[s, p, b, t] * plant_df[plant_df.plant_id .== p, :fuel_price][1]  for s in S for p in P_t for b in B for t in T)
                + sum(shedding[s, b, a, t] * C_shedding + dumping[s, b, a, t] * C_dumping for s in S for b in B for a in A for t in T))
                # + M_eff * sum(frequency_pos[b, t] * maximum(plant_df[:, :fuel_price]) for b in B for t in T)
                # - M_eff * sum(frequency_neg[b, t] * maximum(plant_df[:, :fuel_price]) for b in B for t in T))
                + 1/(B[end]+1) * Δt * π * sum(bypass_penalty * flow_bypass[s, p, b, t] + spill_penalty * flow_spill[s, p, b, t] for s in S for p in P_h for b in B for t in T)
                + π * sum((-volume_end[s, p, T[end]] + volume_end[s, p, 0]) * plant_df[plant_df.plant_id .== p, :fuel_price][1] for s in S for p in P_h))
                # + sum(startup[p, t] * C_startup for p in P_t for t in T))
    optimize!(model)
    obj_value = objective_value(model)
    println("Objective value: $obj_value")
    return model

end


