## Graham Williams
# grwi2594@colorado.edu

# credit to Lawrence Smith for original MATLAB code
# lasm4254@colorado.edu


using BenchmarkTools
using Makie
using GeometryBasics
using LinearAlgebra
using PrettyTables
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
    # Import pointslist, connectivitylist, and edgelist from .txt DelimitedFiles
    PL, CL, EL = importMesh() 

    # Calculate spring rest lengths
    L = Array{Float64, 1}(undef,77)  
    L .=  sqrt.(sum(eachrow((PL[:,EL[1,:]]-PL[:,EL[2,:]]).^2)))

    
    # Calculate nodal masses, assemble mass & damping matrices
    m = zeros(Float64, 96, 1)
    for el ∈ 1:size(CL,2)
        ea = PL[:, CL[2,el]] - PL[:, CL[1,el]]
        eb = PL[:, CL[3,el]] - PL[:, CL[1,el]]
        A = norm(×(ea,eb)/2)
        ne = [3*CL[1,el].-[2; 1; 0] 3*CL[2,el].-[2;1;0] 3*CL[3,el].-[2; 1; 0]]
        m[ne] .+= th*rho*A/3 # mass is calculated as thickness*density*area/3 - three is since 3 points on tri --- iterated for each triangle the mass belongs to 
    end
    m[fixed] .= 1e4
    M = Diagonal(vec(m)) # Mass
    
    C = damp*M  # damping

    # 0'd arrays for holding disp/vel
    d = Array{Float64, 1}(undef, length(PL))
    v = Array{Float64, 1}(undef,length(PL))
    a = Matrix{Float64}(undef,96,1)
    v_predict = Matrix{Float64}(undef,96,1)
    d_predict = Matrix{Float64}(undef,96,1)
    # f = Array{Float64, 1}(undef, length(PL))
    
    
    for t = 0:dt:tsim

        i = reshape(Int[], 0, 2) # index vec for sparse mat
        k = Float64[] # accum vec for sparse mat
        f = zeros(length(PL),1) # 0'd vec for forces

        # Assemble stiffness matrix 
        for ele = 1:size(EL,2)

            # Store coords of end points
            V = [ PL[:, EL[1,ele]]+d[3*EL[1,ele].-[2; 1; 0;]]   PL[:, EL[2,ele]]+d[3*EL[2,ele].-[2; 1; 0;]] ]
            
            # Compute stiffness
            ke = SpringStiffnessMatrix(V,L[ele])
            
            # Nodal degrees of freedom
            ne = vec([3*EL[1, ele].-[2; 1; 0;]' 3*EL[2,ele].-[2;1;0]'])
            
            ts = time_ns()
            je, ie = meshgrid(ne)
            te = time_ns() - ts
            display(te*1e-9)
                     
            i= vcat(i, [ie[:] je[:]])
            k = vcat(k, ke[:])           
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
        
        
        
        # Solve for accelerations at first time step
        if t == 0
            a = inv(M) * (f - (C*v + K*d))
        end
        
        
        
        # Solve for predictors
        v_predict = v + dt * (1-gamma)*a
        d_predict = d + dt * v + dt^2/2 * (1-2*beta)*a 
        
        # Solve for accelerations
        a = inv(M + gamma*dt*C + beta*dt^2*K) * (f - (C*v_predict + K*d_predict))

        # Apply Correctors
        v = v_predict + gamma*dt*a
        d = d_predict + beta*dt^2*a   
  
    end
    
    
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
    A = zeros(sz,sz)
    for i ∈ 1:size(subs,1)
        @inbounds begin
        A[subs[i,1], subs[i,2]] += val[i]
        end
    end
    return A
 end

function meshgrid(x::AbstractVector{T}, y::AbstractVector{T}) where {T}
    m, n = length(y), length(x)
    x = reshape(x, 1, n)
    y = reshape(y, m, 1)
    (repeat(x, m, 1), repeat(y, 1, n))
end

meshgrid(v::AbstractVector) = meshgrid(v,v)

for i ∈ 1:20
    impShell()
end

