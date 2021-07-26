using GeometryBasics
using LinearAlgebra
using DifferentialEquations
using PrettyTables
import GLMakie as mk
import DelimitedFiles: readdlm

# Solution Parameters
global const p = 0.01              # applied pressure
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
        m[ne] .+= (th*rho*A/3)*10000  # mass is calc'd as th*density*area/3 - three is since 3 points on tri --- iterated for each triangle the mass belongs to
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

    # du = similar(s₀)
    # accel(du,s₀, params, 0.0)

    prob = ODEProblem(accel, s₀, (0.0,45.0), params)
    sol = solve(prob)

    
    times = sol.t
    states = sol.u

    plotMesh(states, EL, times)
    # @pt times
    # @pt states


    


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
        ea = (s[3*CL[2,el].-[2; 1; 0;]]) - (s[3*CL[1,el].-[2; 1; 0;]])
        eb = (s[3*CL[3,el].-[2; 1; 0;]]) - (s[3*CL[1,el].-[2; 1; 0;]])
        fe = (p/6)*cross(ea,eb)
        ne = [3*CL[1,el].-[2; 1; 0]; 3*CL[2,el].-[2;1;0]; 3*CL[3,el].-[2; 1; 0]]
        
        f[ne] += repeat(fe,3,1)
    end
    bkpt =1
    # Add spring & damping forces  
    for ele = 1:size(EL,2)

        P = [ (s[3*EL[1,ele].-[2; 1; 0;]])   (s[3*EL[2,ele].-[2; 1; 0;]]) ] # Coords of nodes belonging to this spring
        V = [ s[(3*EL[1,ele].-[2; 1; 0;]).+96]     s[(3*EL[2,ele].-[2; 1; 0;]).+96] ] # Velocities of nodes belonging to this spring

        dL = sqrt(sum((P[:,2] - P[:,1]).^2)) # Length of deformed spring (after t = 0)
        
        disp = dL - L[ele] # Displacement from rest
        relV = (V[:,2] - V[:,1]) # Relative velocity

        n_hat = (P[:,2] - P[:,1])/dL # Directional unit vector
        
        spr = L[ele]*disp # Spring force
        damp = dot((1/c[ele]).*relV, n_hat) # Damping force
          
        ne = [3*EL[1,ele].-[2; 1; 0] 3*EL[2,ele].-[2; 1; 0]] # Indexing for each node within force vector 

        neg = -n_hat
        
        f[ne[:,1]] += -(spr + damp).*n_hat 
        f[ne[:,2]] += -(spr + damp).*(-n_hat) 
        
    end
    bkpt = 1
    a[1:96] = s[97:192]
    a[97:192] = f ./ m

    return a
end

function plotMesh(s, EL, t)
    # fcs = [TriangleFace(val[1], val[2], val[3]) for val in eachcol(CL)]
    # pts = Array{Vector{Point{3, Float32}}, 1}(undef, size(states,1))
    # for i = 1:size(states,1)
    #     pts[i] = [Point3f0(val[1], val[2], val[3]) for val in eachrow(reshape(states[i][1:96],3,32)')]
    # end

    # msh = GeometryBasics.Mesh(pts, fcs)
    # msh = mk.wireframe!(msh)
    # scn = mk.mesh(msh, color = 1 : length(msh.position), camera = cam3d!)

    scene = mk.Scene()
    mk.Camera3D(scene)

    pts = []
    lines = Array{Float64, 2}(undef, length(EL),3)

    for i ∈ 1:size(s,1)
        pts = reshape(s[i][1:96], 3, 32)' 

        for j ∈ 1:size(EL,2)
            lines[2*j-1,:] = s[i][3*EL[1,j].-[2;1;0;]]
            lines[2*j,:] = s[i][3*EL[2,j].-[2;1;0]]
        end
        
        mk.scatter!(scene, pts, markersize=10, color = :red)
        mk.linesegments!(scene, lines, color = :red, linewidth = 2)
        
    end
    display(scene)
    # mk.record(scene, "test.mp4", 1:size(s,1); framerate = 30)
        
    # @pt lines
    # @pt pts

    # display(scn)
end
expShell()