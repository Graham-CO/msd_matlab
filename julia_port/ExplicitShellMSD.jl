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
    
    # Calculate nodal masses, assemble mass & damping vectors
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
    
    k = Float64[] # accum vec for sparse mat NOT USING THIS RN, BC SETTING stiffness == restLength for each spring
    
    # Assemble force vector
    f = zeros(length(PL),1)

    # Calculate resting length of springs
    L = Array{Float64, 1}(undef,77)  
    L .=  sqrt.(sum(eachrow((PL[:,EL[1,:]]-PL[:,EL[2,:]]).^2)))

    # Add pressure forces 
    for el = 1:size(CL,2)
        ea = PL[:, CL[2,el]] - PL[:, CL[1,el]]
        eb = PL[:, CL[3,el]] - PL[:, CL[1,el]]
        fe = (p/6)*cross(ea,eb)
        ne = [3*CL[1,el].-[2; 1; 0] 3*CL[2,el].-[2;1;0] 3*CL[3,el].-[2; 1; 0]]
        f[ne] = repeat(fe,3,1)    
    end
    
    # Add spring forces (should also calc damp forces here) 
    for ele = 1:size(EL,2)

        V = [ PL[:, EL[1,ele]]+d[3*EL[1,ele].-[2; 1; 0;]]   PL[:, EL[2,ele]]+d[3*EL[2,ele].-[2; 1; 0;]] ] # Coords of nodes belonging to this spring
        dL = sqrt(sum((V[:,2] - V[:,1]).^2)) # Length of deformed spring (after t = 0)
          
        ne = [3*EL[1,ele].-[2; 1; 0] 3*EL[2,ele].-[2; 1; 0]] # Indexing for each node within force vector 

        f[ne[:,1]] .+= -L[ele]*(dL - L[ele])*(V[:,2] - V[:,1])/dL # For node1 of edge, f = k*(norm(pt2-pt1) - restLength)*(pt2-pt1)/norm(pt2-pt1)
        f[ne[:,2]] .+= -L[ele]*(dL - L[ele])*(V[:,1] - V[:,2])/dL # Same for node2 of edge, except direction flipped, so that f = k*(norm(pt2-pt1) - restLength)*(pt1-pt2)/norm(pt2-pt1)

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
