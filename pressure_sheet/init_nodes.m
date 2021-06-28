%% Node and Spring Initialization Function
% Graham Williams
% grwi2594@colorado.edu



function [] = init_nodes(mesh, m, r)  
    % initialize node_objs with mesh.Points data
    % fix outer edges
    
    node = Node.empty(length(mesh.Points),0); % initialize empty Node array - avg 0.2s faster over 1000 iterations

    for i = 1:length(mesh.Points)               
        node(i) = Node(mesh.Points(i,:), m, r); % initialize each node

        if node(i).x0(1) == 0 || node(i).x0(1) == 1 || node(i).x0(2) == 0 || node(i).x0(2) == 1
            node(i).fix = true; % fix edge nodes
        end
    end


        
end
