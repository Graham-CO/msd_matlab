%% Node Class
% Graham Williams
% grwi2594@colorado.edu

classdef Node
    properties
        r               % radius - for collision handling (m)
        x0              % initial position (m)
        dx0 = [0,0,0]   % initial velocity (m/s)
        u               % texture coord 1
        v               % texture coord 2
        x               % position (m)
        dx              % velocity (m/s)
        m               % mass (kg)
        fix             % fixed point switch
        f0              % force vector
    end
    
    methods
        % Class constructor
        function node_obj = Node(x0, m, r)
            if nargin > 0 
                node_obj.r = r;
                node_obj.x0 = x0;
                node_obj.m = m;
                node_obj.u = x0(1);
                node_obj.v = x0(2);
                node_obj.f0 = 0;
            end
        end 
    end
end
    