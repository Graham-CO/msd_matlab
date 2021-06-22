%% Node Class
% Graham Williams
% grwi2594@colorado.edu

classdef Node
    properties
        r;              % radius - for collision handling (m)
        x0;             % initial position (m)
        v0 = [0,0,0];   % initial velocity (m/s)
        x;              % position (m)
        v;              % velocity (m/s)
        m               % mass (kg)
        fix;            % fixed point switch
    end
    methods
        function node_obj = Node(x0, m, r, fixbool)
            if nargin < 4
                node_obj.r = r;
                node_obj.x0 = x0;
                node_obj.m = m;
                node_obj.fix = false;
            end
        end
    end
end