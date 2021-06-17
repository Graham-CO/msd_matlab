% DETERMINE THE SPRING CONSTANT k THAT RESULTS IN A POSITIVE MASS
% DISPLACEMENT OF 3.14 after 1.0s


function one_mass_1D()
    clear
    clf
    
    s0 = [0,0,0,0]; %[x,y,vx,vy]
    t0 = linspace(0,50); 
    
    k = 3; %SPRING CONSTANT
    
    
    
    my_ode = @(k)ode23(@(t,s)one_dof_msd(t,s,k),t0,s0)
    
    
    function disp = Fy(expr, k)
        [t,x] = expr(k);
        
        disp = x(end,2)  %NET DISPLACEMENT
    end
     
    cost = @(k)sqrt(Fy(my_ode,k)^2 - 3.14^2);
    
    
     
    opt_k = fminsearch(cost, k)

    
    one_mass_1D_plot(t,x,1)

end

function ds = one_dof_msd(t,s,k)
    L = 1; %RESTING LENGTH (m) 
    c = 0.2; %DAMPING CONSTANT (N*s/m)

    m = 1; %MASS (kg)
    Fext = 10; %APPLIED FORCE (N)
    
    ds = zeros(4,1);
    
    Pm = [s(1), s(2)];
    Vm = [s(3), s(4)];
    
    fxpt = [0,1];
    
    Fm = k *(norm(fxpt - Pm) - L)*mydiv((fxpt - Pm), norm(fxpt - Pm))...
         + c * dot([0,0] - Vm, mydiv((fxpt - Pm), norm(fxpt - Pm))) * mydiv((fxpt - Pm), norm(fxpt - Pm))...
         + [ 0, -Fext];
     
     
    ds(1) = s(3);
    ds(2) = s(4);
    ds(3) = dot(Fm/m, [1,0]);
    ds(4) = dot(Fm/m, [0,1]);
%        
end

function one_mass_1D_plot(t,x,f)
    figure(f);
    p = animatedline('Color', [1 0 0], 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 30);
    xlim([-15 15])
    ylim([-15 15])
    xlabel('X position')
    ylabel('Y position')
    
    for i = 1:length(t)
        x1 = x(i,1);
        y1 = x(i,2);
        
        
        
        clearpoints(p)
        addpoints(p,x1,y1);
        addpoints(p,0,1);
        
        drawnow;
    end
end

function div = mydiv(num,denom)
    if denom < 1e-14
        denom = denom + 1e-10;
    end
    div = num/denom;
end



