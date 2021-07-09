%% Node Class
% Graham Williams
% grwi2594@colorado.edu

classdef Node
    properties
<<<<<<< HEAD
        r;              % radius - for collision handling (m)
        x0;             % initial position (m)
        v0 = [0,0,0];   % initial velocity (m/s)
        x;              % position (m)
        v;              % velocity (m/s)
        m;               % mass (kg)
        fix;            % fixed point switch
=======
        r               % radius - for collision handling (m)
        x0              % initial position (m)
        v0 = [0,0,0]    % initial velocity (m/s)
        u               % texture coord 1
        v               % texture coord 2
        x               % position (m)
        xd              % velocity (m/s)
        m               % mass (kg)
        fix             % fixed point switch
>>>>>>> d48dd1fcece312424d7caffe183119663fe32981
    end
    methods
        function node_obj = Node(x0, m, r, u, v, fixbool)
            if nargin < 6
                node_obj.r = r;
                node_obj.x0 = x0;
                node_obj.u = x0(1);
                node_obj.v = x0(2);
                node_obj.m = m;
                node_obj.fix = false;
            else                            % don't think is needed, but don't want the underline         
                node_obj.r = r; 
                node_obj.x0 = x0;
                node_obj.u = x0(1);
                node_obj.v = x0(2);
                node_obj.m = m;
                node_obj.fix = fixbool;
            end
        end
    end
end