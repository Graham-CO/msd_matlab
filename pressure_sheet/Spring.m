%% Spring Class
% Graham Williams
% grwi2594@colorado.edu

classdef Spring 
    properties
        n1                 % Node 1
        n2                 % Node 2
        k                  % stiffness (N/m)
        c                  % damping (Ns/m)
        bend               % bending spring switch    
        stretch            % stretch spring switch 
    end 
    methods
<<<<<<< HEAD
        function spr_obj = Spring(n1, n2, k, bendbool, stretchbool)    % Class Constructor
            if nargin < 5 
=======
        function spr_obj = Spring(n1, n2, bendbool, stretchbool, k, c)    % Class Constructor
            if nargin < 5
>>>>>>> d48dd1fcece312424d7caffe183119663fe32981
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
        
        function stretchC = stretch_constraint(spr_obj)
            if stretchbool == 1 && isempty(spr_obj.n1.x) % if s_spr & simulation hasn't been performed
                w = spr_obj.n2(1) - spr_obj.n1(1); 
                
            end
        end         
    end
end
        
        