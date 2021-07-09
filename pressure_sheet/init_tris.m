%% Triangle initialization function
% Graham Williams
% grwi2594@colorado.edu

function [tris] = init_tris(mesh, nodes)
    cl = mesh.ConnectivityList; % defines each triangle by 3 vertex IDs & faceID (for tracking textureVertexIDs)
    
    tris = Tri.empty(length(cl), 0); % initialize empty Tri array 
    
    for i = 1:length(cl)
        tris(i) = Tri(nodes(cl(i,1)), nodes(cl(i,2)), nodes(cl(i,3)), i); 
    end
        
end