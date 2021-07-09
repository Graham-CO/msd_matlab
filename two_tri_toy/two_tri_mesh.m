%% Two Triangle Mesh
% Graham Williams
% grwi2594@colorado.edu

function [mesh] = two_tri_mesh()
    [x,y] = meshgrid(0:15, 0:15);
    z = zeros(16);
    T = delaunay(x,y);
    mesh = triangulation(T, x(:), y(:), z(:));
    trisurf(mesh);
end