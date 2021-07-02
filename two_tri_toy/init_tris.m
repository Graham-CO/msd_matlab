%% Triangle initialization function
% Graham Williams
% grwi2594@colorado.edu

function [tris] = init_tris(mesh, nodes)
    CL = mesh.ConnectivityList; % defines each triangle by 3 vertexIDs
                               % each row is corresponds to faceID
                               
    tris = Tri.empty(length(CL), 0); % initialize empty Tri array
    
    for i = 1:height(CL)
        tris(i) = Tri(nodes(CL(i,1)), nodes(CL(i,2)), nodes(CL(i,3)), i);
    end
end