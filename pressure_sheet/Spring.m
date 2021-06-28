%% Spring Class
% Graham Williams
% grwi2594@colorado.edu

classdef Spring < Tri
    properties
        n1;                 % Node 1
        n2;                 % Node 2
        k;                  % stiffness (N/m)
        c;                  % damping (Ns/m)
        bend;               % bending spring switch    
        stretch;            % stretch spring switch 
    end 
    methods
        function spr_obj = Spring(n1, n2, bendbool, stretchbool, k, c)    % Class Constructor
            if nargin < 5
                spr_obj.n1 = n1;                    
                spr_obj.n2 = n2;
                spr_obj.k = 1; % if k not specified, default to 1
                spr_obj.c = 1;
                spr_obj.bend = bendbool;
                spr_obj.stretch = stretchbool;
            else
                spr_obj.n1 = n1;                    
                spr_obj.n2 = n2;
                spr_obj.k = k; % need to have a good way of adjusting k
                spr_obj.c = c;
                spr_obj.bend = bendbool;
                spr_obj.stretch = stretchbool;
            end
        end
        
        function stretchC = stretch_constraint(n1, n2)
            w = n2(1) - n1(1); % along direction 1
        end         
    end
end
        
        