%% 1D MSD
%Lawrence Smith | lasm4254@colorado.edu
clear; clc; close all

%derivative-free simplex optimization of dynamic system
k_opt = fminsearch(@blindOpt,1)

function C = blindOpt(k)

    %initialize physical parameters
    c = 0.5;    %[N*s/m]    damping coefficient
    m = 1.0;    %[kg]       mass

    tspan = [0 10];     %duration of simulation
    X0      = [0; 0;];  %intial conditions: zero position, zero velocity

    %integrate system forward in time
    [t, y] = ode23(@(t,X) forward(t,X,m,c,k), tspan, X0);

    %compute cost as squared distance from target displacement at t_final
    C = (pi - y(end,1))^2;
    
end

function [Xdot] = forward(t,X,m,c,k)

    %derivative of state vector
    Xdot = [ X(2);
            1/m * (-c*X(2) - k*X(1) + 10)];

end
