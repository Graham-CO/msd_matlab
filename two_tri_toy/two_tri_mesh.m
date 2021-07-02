%% Two Triangle Mesh
% Graham Williams
% grwi2594@colorado.edu

function [mesh] = two_tri_mesh()
    [x,y] = meshgrid(0:1, 0:1);
    z = zeros(2);
    T = delaunay(x,y);
    mesh = triangulation(T, x(:), y(:), z(:));
    trisurf(mesh);
end