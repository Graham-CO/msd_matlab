using GeometryBasics
using LinearAlgebra
# using DifferentialEquations
import DelimitedFiles: readdlm
using Plots

# Solution Parameters
global const p = 1              # applied pressure
global const fixed = 1:48       # fixed nodes (edges)
global const th = 0.1           # material thickness
global const rho = 10           # material density
global const damp = 1.5         # damping coefficient

function expShell()
    # Import pointslist, connectivitylist, and edgelist
    PL, CL, EL = importMesh()
    
    


    # Calculate nodal masses, assemble mass & damping matrices
    m = zeros(Float64, 96, 1)
    for el ∈ 1:size(CL,2)
        ea = PL[:, CL[2,el]] - PL[:, CL[1,el]]
        eb = PL[:, CL[3,el]] - PL[:, CL[1,el]]
        A = norm(cross(ea,eb)/2)
        ne = [3*CL[1,el].-[2; 1; 0] 3*CL[2,el].-[2; 1; 0] 3*CL[3,el].-[2; 1; 0]]
        m[ne] .+= th*rho*A/3    # mass is calc'd as th*density*area/3 - three is since 3 points on tri --- iterated for each triangle the mass belongs to
    end
    m[fixed] .= 1e4
    
    c = damp*m # damping
    
    # 0'd arrays for holding disp/vel/accel
    d = zeros(Float64, length(PL), 1)
    v = zeros(Float64, length(PL), 1)
    a = zeros(Float64, length(PL), 1)
    
    k = Float64[] # accum vec for sparse mat
    f = zeros(length(PL),1)


    # Calculate resting length of springs
    L = Array{Float64, 1}(undef,77)  
    L .=  sqrt.(sum(eachrow((PL[:,EL[1,:]]-PL[:,EL[2,:]]).^2))) 

    
    
    
    # Assemble force vector

    # Add pressure forces 
    for el = 1:size(CL,2)
        ea = PL[:, CL[2,el]] - PL[:, CL[1,el]]
        eb = PL[:, CL[3,el]] - PL[:, CL[1,el]]
        fe = (p/6)*cross(ea,eb)
        ne = [3*CL[1,el].-[2; 1; 0] 3*CL[2,el].-[2;1;0] 3*CL[3,el].-[2; 1; 0]]
        f[ne] = repeat(fe,3,1)    
    end
    
    # Add spring forces (and damping soon) 
    for ele = 1:size(EL,2)

        V = [ PL[:, EL[1,ele]]+d[3*EL[1,ele].-[2; 1; 0;]]   PL[:, EL[2,ele]]+d[3*EL[2,ele].-[2; 1; 0;]] ] # Coords of nodes belonging to this spring
        dL1 = sqrt(sum((V[:,2] - V[:,1]).^2)) # That minus this for node 1
        dL2 = sqrt(sum((V[:,1] - V[:,2]).^2)) # That minus this for node 2  

        ne = [3*EL[1,ele].-[2; 1; 0] 3*EL[2,ele].-[2; 1; 0]] # 
        f[ne[:,1]] .+= -L[ele]*(dL1 - L[ele])
        f[ne[:,2]] .+= -L[ele]*(dL2 - L[ele])
    
    end
    

    # params = (m, f, c, k)

    # function accel(a, v, d, params) 
    #     m, f, c, k = params

    #     @. a = 1/m  * (f - (c*v + k*d))
        
    #     return
    # end
    


    # prob = SecondOrderODEProblem(accel, v, d, params, (0.0,10.0))
    # sol = solve(prob)


end

function importMesh()

    PL = readdlm("points.txt", ',', Float64)'    # Points list
    CL = readdlm("tris.txt", ',', Int8)'         # Connectivity list
    EL = readdlm("edges.txt", ',', Int8)'        # Edge List

    return PL, CL, EL
end


expShell()
