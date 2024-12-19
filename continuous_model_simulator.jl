using XLSX
using DataFrames
using JuMP
using CPLEX

include("C:/Users/vegardvk/vscodeProjects/bernstein/get_hydro_data.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/helper_functions.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/find_bernstein_weights.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/plot_results.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/elevation_matrix.jl")



function read_input_data()
    # Read input data
    wind_df = DataFrame(XLSX.readtable("output/wind_ts_weights.xlsx", "Sheet1", infer_eltypes=true))
    hydro_df = DataFrame(XLSX.readtable("output/hydro_data.xlsx", "Sheet1", infer_eltypes=true))
    inflow_df = DataFrame(XLSX.readtable("output/inflow_weights.xlsx", "Sheet1", infer_eltypes=true))
    load_df = DataFrame(XLSX.readtable("output/load_weights.xlsx", "Sheet1", infer_eltypes=true))
    prod_df = DataFrame(XLSX.readtable("output/production_weights.xlsx", "Sheet1", infer_eltypes=true))
    shedding_df = DataFrame(XLSX.readtable("output/shedding_weights.xlsx", "Sheet1", infer_eltypes=true))
    dumping_df = DataFrame(XLSX.readtable("output/dumping_weights.xlsx", "Sheet1", infer_eltypes=true))
    plant_df = DataFrame(XLSX.readtable("output/plant_data.xlsx", "Sheet1", infer_eltypes=true))
    line_df = DataFrame(XLSX.readtable("output/line_data.xlsx", "Sheet1", infer_eltypes=true))

    A = unique(plant_df.area)
    P = unique(plant_df.plant_id)
    T = unique(load_df.timestep)
    L = unique(line_df.line_id)
    B = unique(load_df.b)
    S = unique(load_df.scenario)

    P_w = unique(plant_df[plant_df.fuel_type .== "Wind", :plant_id])
    P_t = unique(plant_df[plant_df.fuel_type .== "Thermal", :plant_id])
    P_h = unique(plant_df[plant_df.fuel_type .== "Hydro", :plant_id])
    # println("Thermal: ", P_t)
    # println("Wind: ", P_w)
    # println("Hydro: ", P_h)
    area_grouped_plants = groupby(plant_df, :area)
    P_a = Dict(g.area[1] => unique(g.plant_id) for g in area_grouped_plants)
    
    # Dicts with the lines that go in/out of area a
    L_in = Dict(a => unique(line_df[line_df.area_to .== a, :line_id]) for a in A) 
    L_out =  Dict(a => unique(line_df[line_df.area_from .== a, :line_id]) for a in A)

    I_disch, I_spill, I_bypass = find_connected_plants(hydro_df)

    C_shedding, C_dumping, C_startup = get_cost_parameters()
    return wind_df, hydro_df, inflow_df, load_df, plant_df, line_df, A, P, T, L, B, P_w, P_t, P_h, P_a, L_in, L_out, I_disch, I_spill, I_bypass, C_shedding, C_dumping, C_startup, prod_df, shedding_df, dumping_df, S
end

function define_and_solve_model(S_input=[])
    wind_df, hydro_df, inflow_df, load_df, plant_df, line_df, A, P, T, L, B, P_w, P_t, P_h, P_a, L_in, L_out, I_disch, I_spill, I_bypass, C_shedding, C_dumping, C_startup, prod_df, shedding_df, dumping_df, S = read_input_data()
    if length(S_input)>=0
        S = S_input
    end
    Δt = 24/T[end]
    B_I = 0:(B[end]+1)
    B_D = 0:(B[end]-1)
    elevation_matrix = get_elevation_matrix(B[end]-1, B[end])
    k_matrix = get_k_matrix(B[end])

    T_extended = 1:(T[end] + 1)
    model = Model()
    set_optimizer(model, CPLEX.Optimizer)

    @variable(model, prod[s in S[2:end], p in P, b in B, t in T] ≥ 0)
    @variable(model, transmission[l in L, b in B, t in T])

    @variable(model, shedding[s in S[2:end], b in B, a in A, t in T] ≥ 0) 
    @variable(model, dumping[s in S[2:end], b in B, a in A, t in T] ≥ 0)

    @variable(model, flow_disch[s in S[2:end], p in P_h, b in B, t in T] ≥ 0) # Antar at alle moduler bare har ett discharge-segment
    @variable(model, flow_bypass[s in S[2:end], p in P_h, b in B, t in T] ≥ 0)
    @variable(model, flow_spill[s in S[2:end], p in P_h, b in B, t in T] ≥ 0)
    @variable(model, total_flow_in[s in S[2:end], p in P_h, b in B, t in T] ≥ 0)
    @variable(model, total_flow_out[s in S[2:end], p in P_h, b in B, t in T] ≥ 0)
    @variable(model, net_inflow[s in S[2:end], p in P_h, b in B, t in T])
    @variable(model, volume_end[s in S[2:end], p in P_h, t in 0:T[end]] ≥ 0)
    @variable(model, volume[s in S[2:end], p in P_h, b in B_I, t in T])

    @constraint(model, wind_production[s in S[2:end], p in P_w, b in B, t in T], prod[s, p, b, t] == wind_df[(wind_df.plant_id .== p) .& (wind_df.timestep .== t) .&
                                                                                                        (wind_df.scenario .== s) .& (wind_df.b .== b), :wind_power][1])

    @constraint(model, controlled_inflow[s in S[2:end], p in P_h, b in B, t in T], total_flow_in[s, p, b, t] == sum(flow_disch[s, i, b, t] for i in I_disch[p]) 
                                                                                                            + sum(flow_bypass[s, i, b, t] for i in I_bypass[p]) 
                                                                                                            + sum(flow_spill[s, i, b, t] for i in I_spill[p]))

    @constraint(model, controlled_outflow[s in S[2:end], p in P_h, b in B, t in T], total_flow_out[s, p, b, t] == flow_disch[s, p, b, t] + flow_bypass[s, p, b, t] + flow_spill[s, p, b, t])

    @constraint(model, hydro_production[s in S[2:end], p in P_h, b in B, t in T], prod[s, p, b, t] == hydro_df[hydro_df.plant_id .== p, :enekv][1] * flow_disch[s, p, b, t] * 3.6)

    @constraint(model, bypass_ub[s in S[2:end], p in P_h, b in B, t in T], flow_bypass[s, p, b, t] ≤ hydro_df[hydro_df.plant_id .== p, :kap_forb][1])
    @constraint(model, prod_ub[s in S[2:end], p in P_h, b in B, t in T], flow_disch[s, p, b, t] ≤ hydro_df[hydro_df.plant_id .== p, :kap_gen_m3s][1])
    @constraint(model, spill_ub[s in S[2:end], p in P_h, b in B, t in T], flow_spill[s, p, b, t] ≤ hydro_df[hydro_df.plant_id .== p, :kap_spill][1])

    @constraint(model, net_inflow_calculation[s in S[2:end], p in P_h, b in B, t in T], net_inflow[s, p, b, t] == total_flow_in[s, p, b, t] 
                                                                                                + inflow_df[(inflow_df.plant_id .== p) .& (inflow_df.timestep .== t) .& (inflow_df.b .== b), :inflow][1] 
                                                                                                - total_flow_out[s, p, b, t])

    @constraint(model, starting_reservoir[s in S[2:end], p in P_h], volume_end[s, p, 0] == hydro_df[hydro_df.plant_id .== p, :starting_reservoir][1])
    @constraint(model, volume_calculation[s in S[2:end], p in P_h, t in T], volume_end[s, p, t] - volume_end[s, p, t-1] == Δt * (1/(B[end]+1)) * sum(net_inflow[s, p, b, t] for b in B))
    @constraint(model, simple_vol_ub[s in S[2:end], p in P_h, t in T], volume_end[s, p, t] ≤ hydro_df[hydro_df.plant_id .== p, :kap_mag][1])

    @constraint(model, volume_weights[s in S[2:end], p in P_h, b in B_I, t in T], volume[s, p, b, t] == Δt * (1/(B[end]+1)) * sum(net_inflow[s, p, b2, t] for b2 in B[1:b]))
    @constraint(model, volume_ub[s in S[2:end], p in P_h, b in B_I, t in T], volume_end[s, p, t-1] + volume[s, p, b, t] ≤ hydro_df[hydro_df.plant_id .== p, :kap_mag][1])

    @variable(model, status[p in P, t in T_extended], Bin)
    @variable(model, startup[p in P, t in T_extended], Bin)
    @constraint(model, startup_count[p in P_t, t in T_extended[2:end]], startup[p, t] ≥ status[p, t] - status[p, t-1])
    # @constraint(model, min_production[s in S[2:end], p in union(P_t, P_h), t in T, b in B], prod[s, p, b, t] >=  
    #                                                                                         prod_df[(prod_df.scenario .== s) .& (prod_df.plant_id .== p) .& (prod_df.timestep .== t) .& (prod_df.b .== b), :lb][1]) 
    @constraint(model, min_production[s in S[2:end], p in union(P_t, P_h), t in T, b in B], prod[s, p, b, t] >=  
                                                                                            prod_df[(prod_df.scenario .== s) .& (prod_df.plant_id .== p) .& (prod_df.timestep .== t) .& (prod_df.b .== b), :lb][1])                                                                                         
    @constraint(model, max_production[s in S[2:end], p in union(P_t, P_h), t in T, b in B], prod[s, p, b, t] <=  
                                                                                            prod_df[(prod_df.scenario .== s) .& (prod_df.plant_id .== p) .& (prod_df.timestep .== t) .& (prod_df.b .== b), :ub][1])

    # @constraint(model, unit_commitment1a[s in S[2:end], p in P_t, t in T], prod[s, p, 0, t] ≤ plant_df[plant_df.plant_id .== p, :gen_ub][1] * status[p, t])
    # @constraint(model, unit_commitment1b[s in S[2:end], p in P_t, t in T], prod[s, p, 0, t] ≥ plant_df[plant_df.plant_id .== p, :gen_lb][1] * status[p, t])

    # @constraint(model, unit_commitment2a[s in S[2:end], p in P_t, t in T], prod[s, p, B[end], t] ≤ plant_df[plant_df.plant_id .== p, :gen_ub][1] * status[p, t+1])
    # @constraint(model, unit_commitment2b[s in S[2:end], p in P_t, t in T], prod[s, p, B[end], t] ≥ plant_df[plant_df.plant_id .== p, :gen_lb][1] * status[p, t+1])

    # @constraint(model, unit_commitment3a[s in S[2:end], p in P_t, t in T], prod[s, p, 1, t] ≤ plant_df[plant_df.plant_id .== p, :gen_ub][1] * status[p, t])
    # @constraint(model, unit_commitment3b[s in S[2:end], p in P_t, t in T], prod[s, p, 1, t] ≥ plant_df[plant_df.plant_id .== p, :gen_lb][1] * status[p, t])

    # @constraint(model, unit_commitment4a[s in S[2:end], p in P_t, t in T], prod[s, p, B[end]-1, t] ≤ plant_df[plant_df.plant_id .== p, :gen_ub][1] * status[p, t+1])
    # @constraint(model, unit_commitment4b[s in S[2:end], p in P_t, t in T], prod[s, p, B[end]-1, t] ≥ plant_df[plant_df.plant_id .== p, :gen_lb][1] * status[p, t+1])

    # @constraint(model, continuity_constraint1[s in S[2:end], p in P_t, t in T[1:end-1]], prod[s, p, B[end], t] == prod[s, p, 0, t+1])
    # @constraint(model, continuity_constraint2[s in S[2:end], p in P_t, t in T[1:end-1]], prod[s, p, B[end], t] - prod[s, p, B[end]-1, t] == prod[s, p, 1, t+1] - prod[s, p, 0, t+1])


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

    @constraint(model, energy_balance[s in S[2:end], b in B, t in T, a in A], 
                # - M_eff * frequency[b, t]
                # - D_eff * sum(frequency_d[j, t] * elevation_matrix[j+1, b+1] for j in B_D)
                + sum(prod[s, p, b, t] for p in P_a[a])
                + 0.99 * sum(transmission[l, b, t] for l in L_in[a])
                - sum(transmission[l, b, t] for l in L_out[a])
                + shedding[s, b, a, t]
                - dumping[s, b, a, t]
                == load_df[(load_df.scenario .== s) .& (load_df.area .== a) .& (load_df.timestep .== t) .& (load_df.b .== b), :Forbruk][1])



    @constraint(model, max_transmission[l in L, b in B, t in T], 0 ≤ transmission[l, b, t] ≤  line_df[line_df.line_id .== l, :capacity][1])
    # @constraint(model, min_transmission[l in L, b in B, t in T], transmission[l, b, t] ≥ -line_df[line_df.line_id .== l, :capacity][1])

    # @objective(model, Min, sum(prod[p, b, t] * get_bernstein_val(bernstein_degree, b, s/sampling_points) * p_dict["C"][p] for p in P for b in B for t in T for s in 0:sampling_points))
    @objective(model, Min, 
                1/(B[end]+1) * Δt * (sum(prod[s, p, b, t] * plant_df[plant_df.plant_id .== p, :fuel_price][1]  for s in S[2:end] for p in P_t for b in B for t in T)
                + sum(shedding[s, b, a, t] * C_shedding + dumping[s, b, a, t] * C_dumping for s in S[2:end] for s in S[2:end] for b in B for a in A for t in T))
                # + M_eff * sum(frequency_pos[b, t] * maximum(plant_df[:, :fuel_price]) for b in B for t in T)
                # - M_eff * sum(frequency_neg[b, t] * maximum(plant_df[:, :fuel_price]) for b in B for t in T))
                + sum((-volume_end[s, p, T[end]] + volume_end[s, p, 0]) * plant_df[plant_df.plant_id .== p, :fuel_price][1] for s in S[2:end] for p in P_h))
                # + sum(startup[p, t] * C_startup for p in P_t for t in T))
    optimize!(model)
    return model

end


