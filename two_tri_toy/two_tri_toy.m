%% Two Triangle Implementation of Baraff/Witkin Paper
% Graham Williams
% grwi2594@colorado.edu

% Generate mesh of two triangles
mesh = two_tri_mesh();

%% TODO: Calculate node mass based on triangle mass 
%%       node mass is sum of 1/3 of mass of each triangle it belongs to
%%       triangle mass is ~ material density * fixed area in uv plane

% Physical Parameters
m = 1;      % mass (kg)
r = 0.2;    % radius - for collision handling
k = 1;      % stiffness (N/m)
c = 1;      % damping (Ns/m)

% Initialize node objects
nodes = init_nodes(mesh, m, r);

% Initialize triangle objects
tris = init_tris(mesh, nodes);



calcStretch(tris(1))
calcStretch(tris(2))

nodes(1).f0 % these are 0 right now because constraint is satisfied...
nodes(2).f0 % will change upon application of pressure forces.. HOPEFULLY

% Ps = mesh.Points
% CL = mesh.ConnectivityList



