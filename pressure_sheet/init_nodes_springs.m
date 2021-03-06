%% Node and Spring Initialization Function
% Graham Williams
% grwi2594@colorado.edu



function [] = init_nodes_springs(mesh, m, r)  
    % initialize node_objs with mesh.Points data
    % fix outer edges
    
    node = Node.empty(length(mesh.Points),0); % initialize empty Node array - avg 0.2s faster over 1000 iterations

    for i = 1:length(mesh.Points)               
        node(i) = Node(mesh.Points(i,:), m, r); % initialize each node

        if node(i).x0(1) == 0 || node(i).x0(1) == 1 || node(i).x0(2) == 0 || node(i).x0(2) == 1
            node(i).fix = true; % fix edge nodes
        end
    end
    
    tic()
    % initialize spr_obj with node data & mesh.ConnectivityList data
    n_nodes = length(node);

    % each row in mesh.Points is a vertex ID
    % each element in mesh.ConnectivityList is a vertex ID

    % to generate springs, first sort connectivity list (for visual
    % debugging)
    sortlist1 = sortrows(mesh.ConnectivityList,1) % by first colmun
    sortlist2 = sortrows(mesh.ConnectivityList,2); % by second column
    sortlist3 = sortrows(mesh.ConnectivityList,3); % by third column

    % use a map (similar to dict) to gen key for each vertex ID, and add vertex
    % IDs to which it has already connected a spring to, checking against the
    % values each time to ensure no duplicates
    map = containers.Map('KeyType', 'uint32', 'ValueType', 'any'); % each key is 1 vertex ID mapped to the value of vertex ID that it attaches a spring to
                                                                   % can't initialize map with only keys and not values, need both or neither
                                                                     
    for i = 1:n_nodes % for every node
        vals = []; % initialize valset to check against, add the values to map at the end
        for j = 1:n_nodes % check every row        
            if i == sortlist1(j,1) %add vals from sort_list1 
                if ~ismember(sortlist1(j,2), vals) % if the vertex ID in position 2 hasn't been added to the valset
                    vals(end+1) = sortlist1(j,2);
                end
                if ~ismember(sortlist1(j,3),vals) % if the vertex ID in position 3 hasnt been added to the valset
                    vals(end+1) = sortlist1(j,3);
                end
            end
            
            if i == sortlist2(j,2) % add vals from sort_list2
                if ~ismember(sortlist2(j,1), vals) % if the vertex ID in position 1 hasn't been added to the valset
                    vals(end+1) = sortlist2(j,1);
                end
                if ~ismember(sortlist2(j,3), vals) % if the vertex ID in position 3 hasn't been added to the valset
                    vals(end+1) = sortlist2(j,3);
                end
            end
            
            if i == sortlist3(j,3) % add vals from sort_list3
                if ~ismember(sortlist3(j,1), vals) % if the vertex ID in position 1 hasn't been added to the valset
                    vals(end+1) = sortlist3(j,1);
                end
                if ~ismember(sortlist3(j,2), vals) % if the vertex ID in position 2 hasn't been added to the valset
                    vals(end+1) = sortlist3(j,2);
                end    
            end
        end
        map(i) = vals; % map now has 1 key for every node, each of which maps to the nodes to which that node is connected
    end
    toc()
    map(3)
    map(22)
    % initialize springs using the generated map
%     for i = 1:n_nodes
        
end
