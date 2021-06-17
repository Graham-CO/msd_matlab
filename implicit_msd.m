% IMPLICIT EULER TOY PROBLEM

% STEP 1: DISPLAY A PIECE OF CLOTH

clear; clc; close all

% 4 NODES - DIAGONAL LINKS INCLUDED

% NODE IS A struct DATATYPE WITH FIELDS:
% x0, y0, m, c, k, s0, rest(only 2D for forces)
node = struct([]); % this isn't required, but I wish it was for clarity...


nnodes = 4; % number of nodes

%------------------------------------------------------------------------%
%                         ASSIGN INITIAL POSITIONS

% TOP LEFT
node(1).x0 = 0;
node(1).y0 = 0;

% TOP RIGHT
node(2).x0 = 0;
node(2).y0 = 1;

% BOTTOM RIGHT
node(3).x0 = 1;
node(3).y0 = 1;

% BOTTOM LEFT
node(4).x0 = 1;
node(4).y0 = 0;


%------------------------------------------------------------------------%
%                      ASSIGN PHYSICAL PARAMETERS

% SET THESE VALUES TO WHATEVER YOU WANT
% INITIALLY ASSUMING UNIFORM MASS FOR NODES, SO DON'T HAVE TO CALCULATE
% FROM TRIANGLE MASS - IF NOT UNIFORM, MASS OF EACH NODE IS 1/3 MASS OF
% TRIANGLE IT BELONGS TO
m = 1; % mass (kg)
k = 1; % stiffness (N/m)
c = 1; % damping (Ns/m)

for i = 1:nnodes
    node(i).m = m;
    node(i).k = k;
    node(i).c = c;
end

%------------------------------------------------------------------------%
%       ASSIGN INITIAL STATE VECTORS (s0) & 2D REST STATES (rest)

% [x position, y position, x velocity, y velocity]
% [x position, y position]

for i = 1:nnodes
    node(i).s0 = [node(i).x0, node(i).y0, 0, 0];
end

for i = 1:nnodes
    node(i).rest = [node(i).x0, node(i).y0];
end
%------------------------------------------------------------------------%
%                       DEFINE TRIANGLES FOR MESH

% tri IS A struct DATATYPE WITH FIELDS:
% n1, n2, n3, nor
tri = struct([]);

tri(1).n1 = node(1).rest;
tri(1).n2 = node(3).rest;
tri(1).n3 = node(4).rest;

tri(2).n1 = node(1).rest;
tri(2).n2 = node(3).rest;
tri(2).n3 = node(2).rest;


%------------------------------------------------------------------------%
%                       DEFINE CONSTRAINT FORCES
% STRETCH

%------------------------------------------------------------------------%
%                      IMPLICIT EULER IMPLEMENTATION
%------------------------------------------------------------------------%
%                     PLOT AS POINTS WITH LINES ATTACHING
[X,Y] = meshgrid(-5:1:5);
%------------------------------------------------------------------------%
%                        GIVE LATTICE TEXTURE
    