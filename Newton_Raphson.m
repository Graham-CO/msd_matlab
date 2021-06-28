%% NEWTON-RAPHSON METHOD
% Graham Williams
% grwi2594@colorado.edu
% MACLab Research

% BASED ON ASSUMPTION THAT SECOND DERIVATIVE IS CONTINUOUS

% SINCE FUNCTION IS QUADRATIC, THERE WILL BE 2 ROOTS
% HOWEVER, N.R. WILL RETURN ONLY ONE
% WHICH ROOT IS RETURNED IS DEPENDENT ON THE INITIAL GUESS 
% i.e. APPROACH FROM LEFT --> -3 ||| APPROACH FROM RIGHT --> -1

clear; clc; close all
%% PROBLEM STATEMENT
% FIND THE ROOT OF THE FUNCTION BELOW USING THE N.R. METHOD.
% GENERATE A PLOT OF THE FUNCTION
% PLOT THE INITIAL GUESS
% PLOT THE UPDATED GUESS AT EACH ITERATION

% f(x) = (x+2)^2 - 1
% df/dx = 2*(x+2) = 2x + 4

func = @(x) (x+2).^2 - 1;   % DEF THE FUNCTION ANONYMOUSLY TO MAKE PLOTTING
                            % GUESSES EASIER AND BUILD FAMILIARITY
                            
deriv = @(x) 2*x + 4;       % ANALYTICAL EXPRESSION FOR DERIVATIVE OF FUNCTION

x = -15:0.1:10;
y = func(x);


plot_2D(y,x,1)
newt_rap(func,deriv)

function plot_2D(func,inp,fig)
    %% 2D LINE PLOT

    figure(fig)
    
    plot(inp,func,'r')
    
    hold on;
end

function plot_guess(out,inp,fig,tol)
    figure(fig)
    init_guess = 0;
    
    if inp == init_guess     % IF INITIAL GUESS, PLOT AS CYAN
        plot(inp,out, 'c *', 'MarkerSize', 10)
    elseif out < tol         % IF CONVERGED, PLOT AS GREEN, ELSE BLUE
        plot(inp,out, 'g *', 'MarkerSize', 10)
    else
        plot(inp,out, 'b x')
        hold on;
    end
end

function newt_rap(func, deriv)
    init_guess = 0;     % INITIAL GUESS
    tol = 1e-5;         % TOLERANCE
    iter = 0;           % INITIALIZE ITERATION COUNTER
    max_iter = 100;      % STOP CONDITION BASED ON MAX ITER
    
    if iter == 0
        y = func(init_guess);
        plot_guess(y,init_guess,1,tol)
    end
    
    x = init_guess;
    while abs(y) > tol && iter <= max_iter
        
        y = func(x);
        dy = deriv(x);
        
        x = x - y/dy;    % UPDATE GUESS ACCORDING TO
                        % x_1 = x_0 - err - f(x_0-err)/f'(x_0-err)
                               
        
        plot_guess(func(x),x,1,tol)
        
        iter = iter + 1;
    end
        
    fprintf('Newton-Raphson method converged on a root value of %f in %i iterations.', x, iter)
end




    
    

    