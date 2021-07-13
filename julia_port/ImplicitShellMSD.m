%% Implicit Shell MSD Sim
%Using Linear Springs, Newmark Integration, no bending stiffness (yet)
%Lawrence Smith | lasm4254@colorado.edu
%1 July 2021

clear; clc; close all

load mesh.mat
PT = mesh.Points;
CM = mesh.ConnectivityList;

%Set up solution parameters
tsim = 10;                  %total simulation time
dt = 0.1;                   %timestep
p = 1;                      %applied pressure
fixed = 1:48;               %fixed nodes (the edges)
th = 0.1;                   %material thickness
rho = 10;                   %material density
beta = 0.25;                %Newmark Integration Beta Parameter
gamma = 0.5;                %Newmark Integration Gamma Parameter
damp = 1.5;                 %damping coefficient

%Extract edges from mesh
e = edges(mesh);

%Compute initial lengths of each spring
L = sqrt(sum((PT(e(:,1),:)-PT(e(:,2),:))'.^2))';

%Compute the mass of each particle
m = zeros(numel(PT),1);
tic
for el = 1:size(CM,1)                        %index thru triangles
    ea = PT(CM(el,2),:)-PT(CM(el,1),:);         %one edge vector
    eb = PT(CM(el,3),:)-PT(CM(el,1),:);         %another edge vector
    A = norm(cross(ea,eb)/2);                     %area of this triangle
    ne = [3*CM(el,1)-[2 1 0] 3*CM(el,2)-[2 1 0] 3*CM(el,3)-[2 1 0]];
    m(ne) = m(ne) + th*rho*A/3;   %account for this mass
    disp(m(ne))
end
toc


%Adjust the masses of the fixed nodes
m(fixed) = 1e4;


%Form global mass matrix

M = diag(m);

%Form global damping matrix
C = damp*M;

%Initialize a matrix of zeros to hold displacement/velocity solutions
d = zeros(numel(PT),1);
v = zeros(numel(PT),1);

tic
for t = 0:dt:tsim
    
    i = []; %initialize indexing vector for sparse matrix
    k = []; %initialize accumulate vector for sparse matrix
    f = zeros(numel(PT),1);
    
    %Assemble stiffness matrix
    
    for el = 1:size(e,1)
        
        V = [PT(e(el,1),:) + d(3*e(el,1)-[2 1 0])'; 
             PT(e(el,2),:) + d(3*e(el,2)-[2 1 0])';]; %store the coords of the endpoints  
        
        ke = SpringStiffnessMatrix(V,L(el)); %compute stiffness
        
        %identify the nodal degrees of freedom
        ne = [3*e(el,1)-[2 1 0] 3*e(el,2)-[2 1 0]];
        [je, ie] = meshgrid(ne);
        
        %Store the indices and stiffness values into arrays
        i = [i; [ie(:) je(:)]];
        k = [k; ke(:)];
        
    end
    
    %accumulate the stiffness array
%     K = sparse(i(:,1),i(:,2),k); %I wonder if sparse matrix is faster?

    K = accumarray(i,k);

   
    i = []; %initialize indexing vector for sparse matrix
    
    %Assemble the force vector
    
    for el = 1:size(CM,1)
        ea = PT(CM(el,2),:)-PT(CM(el,1),:);
        eb = PT(CM(el,3),:)-PT(CM(el,1),:);
        fe = (p/6)*cross(ea,eb);
        ne = [3*CM(el,1)-[2 1 0] 3*CM(el,2)-[2 1 0] 3*CM(el,3)-[2 1 0]];
        f(ne) = repmat(fe,1,3);
    end
    
    
    if t == 0 %Solve for accelerations at first timestep
        a = inv(M) * (f - (C*v + K*d));
    end
    
    %Solve for Predictors
    v_predict = v + dt * (1-gamma)*a;
    d_predict = d + dt * v + dt^2/2 * (1-2*beta)*a;
    
    
    %Solve for Acceleration
    a = inv(M + gamma*dt*C + beta*dt^2*K) * (f - (C*v_predict + K*d_predict));
    
    %Apply correctors
    v = v_predict + gamma*dt*a;
    d = d_predict + beta*dt^2*a;
      
    %display

    trimesh(triangulation(CM,PT + reshape(d,3,[])'),'facecolor','w','edgecolor','b')
    axis([0 1 0 1 0 1])
    drawnow()


end

toc