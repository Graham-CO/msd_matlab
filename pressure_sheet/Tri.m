%% Triangle Class
% Graham Williams
% grwi2594@colorado.edu

classdef Tri 
    properties
        n1     % node 1
        n2     % node 2
        n3     % node 3
        nor    % triangle's normal
        faceID % for tracking springs
    end
    
    properties (Dependent)
        Area   % triangle's area
    end
    
    methods
        function tri_obj = Tri(n1, n2, n3, faceID) % Class Constructor
            % accepts 4 arguments, the three constiuent nodes and an ID for
            % tracking
            if nargin > 0
                tri_obj.n1 = n1;
                tri_obj.n2 = n2;
                tri_obj.n3 = n3;
                tri_obj.faceID = faceID;
                tri_obj.nor = updateNor(tri_obj);
            end
        end
           
        function [nor] = updateNor(tri_obj) % update function for triangle normal
            if isempty(tri_obj.n1.x) && isempty(tri_obj.n2.x) && isempty(tri_obj.n3.x)  % if the forward integration has not taken place
                tri_obj.n1.x = tri_obj.n1.x0;   % for calculation of init nor, will clear at end 
                tri_obj.n2.x = tri_obj.n2.x0;   % so the forward integration has an empty array to dump into
                tri_obj.n3.x = tri_obj.n3.x0;
                
                v = tri_obj.n1.x - tri_obj.n2.x;    % first tri edge
                w = tri_obj.n1.x - tri_obj.n3.x;    % second tri edge
                
                nor = cross(v,w);   % compute unit normal by taking cross product of 2/3 constiuent edges
                
                tri_obj.n1.x = []; % clear xcaches
                tri_obj.n2.x = [];
                tri_obj.n3.x = [];
            else
                v = tri_obj.n1.x - tri_obj.n2.x;
                w = tri_obj.n1.x - tri_obj.n3.x;
                
                nor = cross(v,w);
            end
        end
        
        function a  = get.Area(tri_obj)  % get function for triangle area             
            a = 0.5 * norm(tri_obj.nor);
        end
    end
end
            
            