%% Node Initialization Function
% Graham Williams
% grwi2594@colorado.edu



function [nodes] = init_nodes(mesh, m, r)  
    % initialize node_objs with mesh.Points data
    % fix outer edges
    
    nodes = Node.empty(length(mesh.Points),0); % initialize empty Node array - avg 0.2s faster over 1000 iterations

    for i = 1:length(mesh.Points)               
        nodes(i) = Node(mesh.Points(i,:), m, r, mesh.Points(i,1), mesh.Points(i,2)); % initialize each node

        if nodes(i).x0(1) == 0 || nodes(i).x0(1) == 1 || nodes(i).x0(2) == 0 || nodes(i).x0(2) == 1
            nodes(i).fix = true; % fix edge nodes
        end
    end      
end
