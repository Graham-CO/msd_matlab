using GeometryBasics
using LinearAlgebra
using DifferentialEquations
import DelimitedFiles: readdlm

# Solution Parameters
global const tsim = 10          # total sim time
global const dt = 0.1           # timestep
global const p = 1              # applied pressure
global const fixed = 1:48       # fixed nodes (edges)
global const th = 0.1           # material thickness
global const rho = 10           # material density
global const beta = 0.25        # Newmark integration beta parameter
global const gamma = 0.5        # Newmark integration gamma parameter
global const damp = 1.5         # damping coefficient

function expShell()
    # Import pointslist, connectivitylist, and edgelist
    PL, CL, EL = importMesh()

    # Calculate spring rest lengths
    L = Array{Float64, 1}(undef, 77)
    L .= sqrt.(sum(eachrow((PL[:, EL[1,:]]-PL[:,EL[2,:]]).^2)))

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
    M = Diagonal(vec(m)) # Mass

    C = damp*M # damping
    
    # 0'd arrays for holding disp/vel
    d = zeros(Float64, length(PL), 1)
    # v = Array{Float64, 1}(undef, length(PL))
    v = zeros(Float64, length(PL), 1)
    

    i = reshape(Int[], 0, 2) # index vec for sparse mat
    k = Float64[] # accum vec for sparse mat
    f = zeros(length(PL),1)

    # Assemble stiffness matrix 
    
    for ele = 1:size(EL,2)

        # Store coords of end points
        V = [ PL[:, EL[1,ele]]+d[3*EL[1,ele].-[2; 1; 0;]]   PL[:, EL[2,ele]]+d[3*EL[2,ele].-[2; 1; 0;]] ]
        
        
        # Compute stiffness
        ke = SpringStiffnessMatrix(V,L[ele])
        

        # Nodal degrees of freedom
        ne = [3*EL[1, ele].-[2; 1; 0;]' 3*EL[2,ele].-[2;1;0]']

        je = [ne; ne; ne; ne; ne; ne;]

        ie = [ne'; ne'; ne'; ne'; ne'; ne']

        tmp = [ie je[:]]

        i = [i; tmp]
        k = [k; ke[:]]
    end
    
    K = accumarray(i,k)
    
    # Assemble force vector 
    for el = 1:size(CL,2)
        ea = PL[:, CL[2,el]] - PL[:, CL[1,el]]
        eb = PL[:, CL[3,el]] - PL[:, CL[1,el]]
        fe = (p/6)*cross(ea,eb)
        ne = [3*CL[1,el].-[2; 1; 0] 3*CL[2,el].-[2;1;0] 3*CL[3,el].-[2; 1; 0]]
        f[ne] = repeat(fe,3,1)    
    end

    params = (M, f, C, K)

    function accel(a, v, d, params, t) 
        M, f, C, K = params
        a .= inv(M) * (f - (C*v + K*d))
        return
    end


    prob = SecondOrderODEProblem(accel, v, d, params, (0.0,10.0))
    sol = solve(prob)

end

function importMesh()

    PL = readdlm("points.txt", ',', Float64)'    # Points list
    CL = readdlm("tris.txt", ',', Int8)'         # Connectivity list
    EL = readdlm("edges.txt", ',', Int8)'        # Edge List

    return PL, CL, EL
end

function SpringStiffnessMatrix(V,k)

    L = sqrt(sum((V[:,1] - V[:,2]).^2))

    Cx = (V[1,1] - V[1,2])/L
    Cy = (V[2,1] - V[2,2])/L
    Cz = (V[3,1] - V[3,2])/L

    A = [Cx^2 Cx*Cy Cx*Cz; Cy*Cx Cy^2 Cy*Cz; Cz*Cx Cz*Cy Cz^2]
    
    K = k*[A -A; -A A]

    return K
end

function accumarray(subs, val; sz=maximum(subs[:,1]))
    counts = Dict()
    
    for i = 1:size(subs,1)
         counts[subs[i,:]]=[get(counts,subs[i,:],[]);val[i...]]    
    end
    
    A = zeros(sz,sz)
    for j = keys(counts)
         A[j...] = sum(counts[j])
    end
    return A
 end

expShell()