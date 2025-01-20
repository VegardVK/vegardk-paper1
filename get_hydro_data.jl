using DataFrames

function add_start_reservoir!(df)
    df[!, "starting_reservoir"] .= 0.0
    for p in df[!, "plant_id"]
        if p == 49927
            df[df.plant_id .== p, "starting_reservoir"] = df[df.plant_id .== p, "kap_mag"] *0.5
        else
            df[df.plant_id .== p, "starting_reservoir"] = df[df.plant_id .== p, "kap_mag"] *0.5
        end
    end
end

function write_inflow_table(hydro_df, timesteps, scenarios)
    df = DataFrame(plant_id=Int[], timestep=Int[], inflow=Float64[])
    for row in eachrow(hydro_df)
        if row["plant_id"] == 49927
            for t in 1:timesteps
                push!(df, (row["plant_id"], t, 0))
            end
        else
            for t in 1:timesteps
                push!(df, (row["plant_id"], t, 0))
            end
        end
    end
    original_df_length = nrow(df)
    df = repeat(df, inner=(scenarios+1))
    df.scenario = repeat(0:scenarios, outer=original_df_length)
    XLSX.writetable("output/inflow_data.xlsx", df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
end

function find_connected_plants(df)
    # Makes the dicts I_disch, I_bypass and I_spill which, for a key p, 
    # gives an array of all plants which discharge, bypass or spill into plant p.
    P = df[:, "plant_id"]
    I_disch = Dict(key => [] for key in P)
    I_bypass = Dict(key => [] for key in P)
    I_spill = Dict(key => [] for key in P)

    for row in eachrow(df)
        push!(get!(I_disch, row["topo_gen"], []), row["plant_id"])
        push!(get!(I_spill, row["topo_flom"], []), row["plant_id"])
        push!(get!(I_bypass, row["topo_forb"], []), row["plant_id"])
    end

    return I_disch, I_spill, I_bypass
end

function set_wv!(I_disch, df, p, energy_price)
    for p2 in I_disch[p]
        df.fuel_price[df.plant_id .== p2] = df.fuel_price[df.plant_id .== p] +  energy_price * df.enekv[df.plant_id .== p2] * 3.6
        # println("Fuel price in $p2 = $energy_price * $(enekv[p2]) + $(fuel_price[p]) = $(fuel_price[p2])")
        if length(I_disch[p2]) > 0
            set_wv!(I_disch, df, p2, energy_price)
        end
    end
end