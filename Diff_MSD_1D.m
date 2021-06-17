%% 1D Differentiable MSD 
% Optimize values of k to reach target trajectory
% Utilizing backpropagation through time

% Graham Williams | grwi2594@colorado.edu
clear; clc; close all

C = diffMSD_offline();
semilogy(C/norm(C));

function [Ct] = diffMSD_offline()

    iter = 0;   % INITIALIZE ITERATION COUNTER
    
    % INITIALIZE LEARNING PARAMETERS 
    nu = 1;
    nu_decay = 0.8;
    
    % INITIALIZE COST & TRACKER FOR PLOTTING
    C = 1;
    Ct = [];

    % INITIALIZE PHYSICAL PARAMETERS
    c = 0.5; %[N*s/m] DAMPING COEFFICIENT
    m = 1.0; %[kg] MASS
    
    % INITIAL GUESS FOR PARAMETER OPTIMIZATION
    k = 1.0;

    tspan = [0 10]; % DURATION OF SIM - SIZE OF BATCH 
    
    S0 = [0;0]; % INITIAL STATE [y,Vy]
    e0 = 0; % INITIAL ERROR
    
    while C > 1e-5 && iter < 40  % STOPPING CONDS - TOL & ITER COUNT
        
        % FORWARD TIME INTEGRATION
        [t,y] = ode23(@(t,S) forward(t,S,m,c,k), tspan, S0);

        % COST FUNCTION - USED FOR PLOTTING, NOT ACTUAL OPTIMIZATION
        C = (pi - y(end,1))^2; % SUM (IF TRACKING VELOCITY) OF MEAN SQUARE ERRORS
        
        % APPEND COST TRACKER
        Ct = [Ct; C];
                                          
        % BACKWARD TIME INTEGRATION
        [s,e] = ode23(@(s,e) backward(s,e,m,c,k,t,y), tspan, e0);

        % DEFINE GRADIENT AS NEG. INTEGRAL FROM 0 --> T (e*J_1*K)
        % COMPUTE USING RECT INTEGRATION
        
        lambda_theta = 0;   % INITIALIZE GRADIENT
        for i = 1:length(s)-1
            [edot, J_0, J_1,K] = backward(s(i),e(i),m,c,k,t,y); % RETURN PARAMETERS THAT MAKE UP THE GRADIENT
            
            ds = s(i+1) - s(i); 
            lambda_theta = lambda_theta - ds*e(i)*inv(J_1)*K; % MINUS BECAUSE NEGATIVE INTEGRAL
        end
        
        % LEARNING PARAMETER DECAY - I'M NOT SURE WHY WE WOULD DECAY ON
        % FIRST ITERATION, BUT IT DOESN'T WORK OTHER WAY
        nu = nu*nu_decay;
        
        % PARAMETER UPDATE EQN
        k = k - nu * lambda_theta ./ norm(lambda_theta);
        
        iter = iter + 1;    % UPDATE ITER COUNT       
    end
end

function [Sdot] = forward(t,S,m,c,k)    % FORWARD DIFF EQ

    g = 9.8; %[m/s^2]

    Sdot = [ S(2);
             1/m*(-c*S(2) - k*S(1) + m*g)]; % DERIV OF STATE - [vy,ay]
         
end

function [edot,J_0,J_1,K] = backward(s,e,m,c,k,t,y) % BACKWARD DIFF EQ
    
    x = interp1(t,y(:,1),t(end)-s); % INTERPOLATE POSITION VALUE AT 
    
    J_0 = -k;   % FIRST JACOBIAN - df/dx
    J_1 = -c;   % SECOND JACOBIAN - df/dxdot
    K = -x;     % DERIV OF FUNC WRT PARAMETER - df/dk
    
    e_bar = 2*(x - pi); % DERIV OF COST WRT STATE - dc/dx
 
    edot = e_bar + e'*inv(J_1)*J_0; % DIFF EQ DESCRIBING e - SET UP FOR MATRIX
    
    edot = edot';   % ode23() REQ COLUMN VEC    
end