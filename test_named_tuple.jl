include("C:/Users/vegardvk/vscodeProjects/bernstein/get_hydro_data.jl")
using XLSX
using DataFrames
using UnPack

function get_connection_sets()
    hydro_df = DataFrame(XLSX.readtable("output/hydro_data.xlsx", "Sheet1", infer_eltypes=true))
    I_disch, I_spill, I_bypass = find_connected_plants(hydro_df)

    I = (disch = I_disch, spill = I_spill, bypass = I_bypass)
    return I
end

@unpack disch, spill = get_connection_sets()
println(disch)