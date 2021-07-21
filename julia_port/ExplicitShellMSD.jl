using GeometryBasics
using LinearAlgebra
using DifferentialEquations
import DelimitedFiles: readdlm

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
    
    # Initialize state vector
    d₀ = reshape(PL, length(PL), 1)
    v₀ = zeros(Float64, length(PL), 1)
    s₀ = vcat(d₀, v₀)

    
    # Initialize force vector
    f = zeros(length(PL),1)

    # Calculate resting length of springs
    L = Array{Float64, 1}(undef,77)  
    L .=  sqrt.(sum(eachrow((PL[:,EL[1,:]]-PL[:,EL[2,:]]).^2)))
    
    params = (PL, CL, EL, L, m, f, c)
    display(params)

    prob = ODEProblem(accel, s₀, (0.0,10.0), params)
    sol = solve(prob)


end

function importMesh()

    PL = readdlm("points.txt", ',', Float64)'    # Points list
    CL = readdlm("tris.txt", ',', Int8)'         # Connectivity list
    EL = readdlm("edges.txt", ',', Int8)'        # Edge List

    return PL, CL, EL
end

function accel(a, s, params, t) 
    PL, CL, EL, L, m, f, c = params

    # Add pressure forces 
    for el = 1:size(CL,2)
        ea = (PL[:, CL[2,el]]+s[3*CL[2,el].-[2; 1; 0;]]) - (PL[:, CL[1,el]]+s[3*CL[1,el].-[2; 1; 0;]])
        eb = (PL[:, CL[3,el]]+s[3*CL[3,el].-[2; 1; 0;]]) - (PL[:, CL[1,el]]+s[3*CL[1,el].-[2; 1; 0;]])
        fe = (p/6)*cross(ea,eb)
        ne = [3*CL[1,el].-[2; 1; 0] 3*CL[2,el].-[2;1;0] 3*CL[3,el].-[2; 1; 0]]
        f[ne] = repeat(fe,3,1)    
    end

    # Add spring & damping forces  
    for ele = 1:size(EL,2)

        P = [ PL[:, EL[1,ele]]+s[3*EL[1,ele].-[2; 1; 0;]]   PL[:, EL[2,ele]]+s[3*EL[2,ele].-[2; 1; 0;]] ] # Coords of nodes belonging to this spring
        V = [ s[(3*EL[1,ele].-[2; 1; 0;]).+96]     s[(3*EL[2,ele].-[2; 1; 0;]).+96] ] # Velocities of nodes belonging to this spring
        
        dL = sqrt(sum((P[:,2] - P[:,1]).^2)) # Length of deformed spring (after t = 0)
        
        disp = dL - L[ele] # Displacement from rest
        relV = (V[:,1] - V[:,2]) # Relative velocity
        n_hat = (P[:,2] - P[:,1])/dL # Directional unit vector
        
        spr = L[ele]*disp # Spring force
        damp = c[ele]*relV*n_hat # Damping force
          
        ne = [3*EL[1,ele].-[2; 1; 0] 3*EL[2,ele].-[2; 1; 0]] # Indexing for each node within force vector 

        f[ne[:,1]] .+= -(spr + damp)*n_hat # For node1 of edge, f = k*(norm(pt2-pt1) - restLength)*(pt2-pt1)/norm(pt2-pt1)
        f[ne[:,2]] .+= -(spr + damp)*-n_hat # Same for node2 of edge, except direction flipped, so that f = k*(norm(pt2-pt1) - restLength)*(pt1-pt2)/norm(pt2-pt1)

    end

    a[1:96] = s[1:96]

    a[97:192] .= f ./ m

    
    return
end
expShell()
