using Makie
using GeometryBasics
using BenchmarkTools
import DelimitedFiles: readdlm





function impShell()
    
    PL, CL, EL = importMesh()
    
    L = Array{Float64, 1}(undef,77)
    L .= sqrt.(sum(eachrow((PL[:,EL[1,:]]-PL[:,EL[2,:]]).^2)))
    
    # Convert to GeometryBasics.Mesh
    pts = [Point3f0(val[1], val[2], val[3]) for val in eachcol(PL)]
    fcs = [TriangleFace(val[1], val[2], val[3]) for val in eachcol(CL)] 
    msh = GeometryBasics.Mesh(pts, fcs)



    # Plot Mesh
    scn = Makie.mesh(msh, color = 1 : length(msh.position))
    
    Makie.wireframe!(msh)

    display(scn)
end

function importMesh()
    PL = readdlm("points.txt", ',', Float64)'    # Points list
    CL = readdlm("tris.txt", ',', Int8)'         # Connectivity list
    EL = readdlm("edges.txt", ',', Int8)'        # Edge List

    return PL, CL, EL
end


impShell()