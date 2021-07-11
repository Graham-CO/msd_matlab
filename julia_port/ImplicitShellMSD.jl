using GeometryBasics: minimum
using Makie
using GeometryBasics
using Clustering
using Triangulate

tstart = time_ns()
function impShell() 
    # Construct a square
    outerpnts = [Point(0.,0.), Point(0., 0.25), Point(0., 0.5), Point(0., 0.75), Point(0.,1.), Point(0.25,0.), Point(0.5, 0.), Point(0.75, 0.), Point(1., 0.),  Point(1., 1.)]
    rect = Rect(Vec(0.0,0.0), Vec(1.0,1.0))

    # Generate points inside of the rectangle
    innerpnts = Point[]



    
    let
        n = 50 # Number of inner nodes
        rmin = (3/(4*π*n))^(1/3)
        maxiter = 10000
        outernodes = [getindex.(outerpnts, 1) getindex.(outerpnts, 2)]
        i = 0
        while length(innerpnts[:,1]) ≤ n && i < maxiter
            i += 1

            innernodes = [getindex.(innerpnts, 1) getindex.(innerpnts, 2)]
            p = Point(rand(2)...)

            r1 = minimum((sum((outernodes .- p').^2, dims=2)).^(1/2))
            r2 = isempty(innernodes) ? 1 : minimum((sum((innernodes .- p').^2, dims=2)).^(1/2))


            if r1 > rmin && r2 > rmin && p[1] > .15 && p[1] < .85 && p[2] > .15 && p[2] < .85
                push!(innerpnts, p)
            end
        end
    end

    outernodes = [getindex.(outerpnts, 1) getindex.(outerpnts, 2)]
    innernodes = [getindex.(innerpnts, 1) getindex.(innerpnts, 2)]
    # innernodes = collect(kmeans(innerpntsmat', 29).centers')



    # Triangulate the points
    triin = TriangulateIO()
    triin.pointlist = [outernodes; innernodes]'
    msh, _ = triangulate("vcDQ", triin)





    # Convert to GeometryBasics.Mesh
    pts = [Point(val[1], val[2]) for val in eachcol(msh.pointlist)]
    fcs = [TriangleFace(val[1], val[2], val[3]) for val in eachcol(msh.trianglelist)]
    msh = GeometryBasics.Mesh(pts, fcs)

    # Plot Mesh
    scn = Makie.mesh(msh, color = 1 : length(msh.position))
    Makie.wireframe!(msh)


    display(scn)
end


impShell()
tend = time_ns() - tstart
print(tend*1e-9)