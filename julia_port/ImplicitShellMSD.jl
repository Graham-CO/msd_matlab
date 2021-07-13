using Makie
using GeometryBasics
using BenchmarkTools
using LinearAlgebra
import DelimitedFiles: readdlm

# Solution Paramaters
global const tsim = 10              # total sim time
global const dt = 0.1               # timestep
global const p = 1                  # applied pressure
global const fixed = 1:48           # fixed nodes (edges)
global const th = 0.1               # material thickness
global const rho = 10               # material density
global const beta = 0.25            # Newmark integration beta parameter
global const gamma = 0.5            # Newmark integration gamma parameter
global const damp = 1.5             # damping coefficient

function impShell()
    
    PL, CL, EL = importMesh() # Import pointslist, connectivitylist, and edgelist from .txt DelimitedFiles

    

    
    L = Array{Float64, 1}(undef,77)  
    L .=  sqrt.(sum(eachrow((PL[:,EL[1,:]]-PL[:,EL[2,:]]).^2)))

    m = Array{Float64, 2}(undef, length(PL), 1)
    
    ts = time_ns()
    for el ∈ 1:size(CL,2)
        ea = PL[:, CL[2,el]] - PL[:, CL[1,el]]
        eb = PL[:, CL[3,el]] - PL[:, CL[1,el]]
        A = norm(×(ea,eb)/2)
        ne = [3*CL[1,el].-[2; 1; 0] 3*CL[2,el].-[2;1;0] 3*CL[3,el].-[2; 1; 0]]
        m[ne] .= th*rho*A/3
    end
    
    te = time_ns() - ts
    display(te*1e-9)
        
    
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

