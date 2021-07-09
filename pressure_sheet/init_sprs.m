function [b_sprs, s_sprs] = init_sprs(mesh, nodes, k, c) % not using k or c yet
    EL = edges(mesh);   % edge list
    
    b_sprs = Spring.empty(length(EL), 0); % initialize empty (bend)Spring array
    s_sprs = Spring.empty(length(EL), 0); % initialize empty (stretch)Spring arary
    
    for i = 1:length(EL) % for each edge
        % initialize bend springs
        bendbool = 1;
        stretchbool = 0;
        b_sprs(i) = Spring(nodes(EL(i,1)), nodes(EL(i,2)), bendbool, stretchbool);
        
        % initialize stretch springs
        bendbool = 0;
        stretchbool = 1;
        s_sprs(i) = Spring(nodes(EL(i,1)), nodes(EL(i,2)), bendbool, stretchbool);        
    end
end