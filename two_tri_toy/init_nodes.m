%% Node Initialization Function
% Graham Williams
% grwi2594@colorado.edu

function [nodes] = init_nodes(mesh, m, r)
    % initialize node_objs with mesh.Points data
    % fix nodes on shared edge
    
    nodes = Node.empty(length(mesh.Points),0);  % initialize empty Node array - for performance
    
    for i = 1:length(mesh.Points)
        nodes(i) = Node(mesh.Points(i,:), m, r); % init each node
        
        % Fix nodes on shared edge
        if nodes(i).x0(1) == 0 && nodes(i).x0(2) == 1
            nodes(i).fix = true; 
        elseif nodes(i).x0(1) == 1 && nodes(i).x0(2) == 0
            nodes(i).fix = true;      
        end
    end
end