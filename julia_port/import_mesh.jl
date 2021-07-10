using DelimitedFiles

function importMesh()
    PL = readdlm("points.txt", ',', Float64)'    # Points list
    CL = readdlm("tris.txt", ',', Int8)'         # Connectivity list
    EL = readdlm("edges.txt", ',', Int8)'        # Edge List
    display(PL)

    return PL, CL, EL
end

importMesh()