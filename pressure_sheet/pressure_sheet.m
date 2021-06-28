%% Pressure Force Sim
% Graham Williams
% grwi2594@colorado.edu

clc; clear; close all

% Simulate uniform pressure force to underside of triangulation() sheet
% Fixed edges, rest free

% Each spring has unique stiffness (k=1 for now)
% Each node has mass (m-1)
% Neglect gravity and damping

% Display deformed shape after uniform pressure of 0.1 for 4s
% plot of volume under sheet vs time

load mesh.mat % Load triangulation() mesh

triplot(mesh) % plot mesh
view([45 45]) % set viewing angle

% Physical Parameters
m = 1;      % mass (kg)
r = 0.2;    % radius (m)
% k = 1;      % stiffness (N/m)
% c = 0.5;    % damping (Ns/m)

% initialize Node class objects
init_nodes(mesh, m, r)

% initialize Tri class objects
init_tris(mesh)
% node(1)
% node(16)
% node(32)
