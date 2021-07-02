%% Triangle Class
% Graham Williams
% grwi2594@colorado.edu

classdef Tri
    properties
        n1;         % node 1
        n2;         % node 2
        n3;         % node 3
        k_stretch;  % stretch stiffness
        nor;        % triangle's normal
        faceID;     % for tracking purposes
    end
    
    properties (Dependent = true)
        Area;       % triangle's area in uv plane
    end
    
    methods
        % Class constructor
        function tri_obj = Tri(n1, n2, n3, faceID)
            if nargin > 0
                tri_obj.n1 = n1;
                tri_obj.n2 = n2;
                tri_obj.n3 = n3;
                tri_obj.k_stretch = 1; % can create new constructor that accepts this as an argument (for tuning/optimization)
                tri_obj.faceID = faceID;
                tri_obj.nor = updateNor(tri_obj);
            end
        end
        
        % Update function for triangle normal
        function [nor] = updateNor(tri_obj)
            % If simulation hasn't begun
            if isempty(tri_obj.n1.x) && isempty(tri_obj.n2.x) && isempty(tri_obj.n3.x)
                
                v = tri_obj.n1.x0 - tri_obj.n2.x0; % first tri edge
                w = tri_obj.n1.x0 - tri_obj.n3.x0; % second tri edge
                
                nor = cross(v,w)/norm(cross(v,w)); % compute unit normal vector
                
            else
                
                v = tri_obj.n1.x - tri_obj.n2.x;
                w = tri_obj.n1.x - tri_obj.n3.x;
                
                nor = cross(v,w)/norm(cross(v,w));
            
            end
        end
        
        % Get function for triangle's area in uv coords
        function [a] = get.Area(tri_obj)
            
            X = [tri_obj.n1.u, tri_obj.n2.u, tri_obj.n3.u];
            Y = [tri_obj.n1.v, tri_obj.n2.v, tri_obj.n3.v];
            
            a = polyarea(X,Y);
            
        end
        
        % Constraint function for stretch
        function [stretchF] = calcStretch(tri_obj)
            % Following section 4.2 of [Baraff, Witkin '98]
            u1 = tri_obj.n1.u;
            u2 = tri_obj.n2.u;
            u3 = tri_obj.n3.u;
            
            v1 = tri_obj.n1.v;
            v2 = tri_obj.n2.v;
            v3 = tri_obj.n3.v;
            
            if isempty(tri_obj.n1.x) % if simulation hasn't begun
                x1 = tri_obj.n1.x0;
                x2 = tri_obj.n2.x0;
                x3 = tri_obj.n3.x0;
            else
                x1 = tri_obj.n1.x;
                x2 = tri_obj.n2.x;
                x3 = tri_obj.n3.x;
            end
            
            dlta_x1 = x2 - x1;
            dlta_x2 = x3 - x1;
            
            dlta_u1 = u2 - u1;
            dlta_u2 = u3 - u1;
            dlta_v1 = v2 - v1;
            dlta_v2 = v3 - v1;
            
            % triangle's uv area
            a = tri_obj.Area;
            
            % stretch stiffness
            k = tri_obj.k_stretch;
            
            % determinant of inverse dlta_u, dlta_v matrix (eq 9)
            detinv = 1/(dlta_u1*dlta_v2 - dlta_u2*dlta_v1);
            
            % biases for controlling anisotropic stretch
            bu = 1; 
            bv = 1;
            
            % w is a (approximated) linear function that maps from plane coords to world space
            wu = (dlta_x1*dlta_v2 - dlta_x2*dlta_v1)*detinv; 
            wv = (dlta_x2*dlta_u1 - dlta_x1*dlta_u2)*detinv;
            
            % inverse magnitude of wu, wv
            wu_im = 1 / norm(wu);
            wv_im = 1 / norm(wv);
            
            % calc w_hat for first and second derivatives of constraint
            wu_hat = wu * wu_im;
            wv_hat = wv * wv_im;
            
            % calc dw/dx
            dwu_dx = zeros(1,3);
            dwv_dx = zeros(1,3);
            
            dwu_dx(1) = (dlta_v1 - dlta_v2) * detinv;
            dwu_dx(2) = dlta_v2 * detinv;
            dwu_dx(3) = -dlta_v1 * detinv;
            
            dwv_dx(1) = (dlta_u2 - dlta_u1) * detinv;
            dwv_dx(2) = -dlta_u2 * detinv;
            dwv_dx(3) = dlta_u1 * detinv;
            
            % calc constraint
            Cu = a * (length(wu) - bu);
            Cv = a * (length(wv) - bv);
            
            % calc first derivative of constraint
            dCu_dx = cell(3,1);
            dCv_dx = cell(3,1);
            for i = 1:3
                dCu_dx{i,1} = (a .* dwu_dx(i)) .* wu_hat;
                dCv_dx{i,1} = (a .* dwv_dx(i)) .* wv_hat;
            end
            
            % plz for the love of god just let me make arrs of vecs MATLAB
            % might need to be reshaped, idk yet...
            dCu_dx = cell2mat(dCu_dx);
            dCv_dx = cell2mat(dCv_dx);
        
            % add force to node force vector
            for i = 1:1
                val = -k * (dCu_dx(i,:)*Cu + dCv_dx(i,:)*Cv)
                
                % pass nodes into Tri as a vec to make this better
                if i == 1
                    tri_obj.n1.f0 = tri_obj.n1.f0 + val;
                end
                if i == 2
                    tri_obj.n2.f0 = tri_obj.n2.f0 + val;
                end
                if i == 3
                    tri_obj.n3.f0 = tri_obj.n3.f0 + val;
                end
            end    
        end
    end
end