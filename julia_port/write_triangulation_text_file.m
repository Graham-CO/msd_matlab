%% Graham Williams
% grwi2594@colorado.edu

%% Credit to Jacob Haimes
% jacob.haimes@colorado.edu


clear; close all; clc;


%% Load previously formatted data
load('mesh.mat')

%% Choose points and connectivity list
% Alternatively input a triangulation
nodes = mesh.Points;
tris = mesh.ConnectivityList;
edges = mesh.edges;

%% Open a text file
% If you would like a different format, designate so here
file_name = 'points.txt';

% NOTE: The 'W' will flush the input file name if it already has
% contents in it
fileID = fopen(file_name,'W');

%% Define number of elements and nodes
num_nodes = size(nodes,1);
num_tris = size(tris,1);
num_edges = size(edges,1);

%% Write nodes in file
% Assuming 3 coordinates
% fprintf(fileID,'// NODES \n// node#, [spatial coords of this node]\n\n');

% Format specification for list of nodes
%       node#, [spatial coords of this node]
format =   '%f,\t%f,\t%f\n';

% Print node number and coordinates
for ii = 1:num_nodes
    contents = [nodes(ii,:)];
    fprintf(fileID,format,contents);
end

fclose(fileID);

file_name = 'tris.txt';

fileID = fopen(file_name, 'W');

%% Write elements (connectivity List) in file
% Assuming 4 nodes per element
% fprintf(fileID,'// CONNECTIVITY LIST\n// elem#, [node#s connected by this element]\n\n');


% Format specification for list of elements
%       elem#, [node#s connected by this element]
format =  '%i,\t%i,\t%i\n';

% Print element number and connectivity
for ii = 1:num_tris
    contents = [tris(ii,:)];
    fprintf(fileID,format,contents);
end

%% Close written input file
fclose(fileID);

file_name = 'edges.txt';

fileID = fopen(file_name, 'W');

format = '%i,\t%i\n';

for ii = 1:num_edges
    contents = [edges(ii,:)];
    fprintf(fileID,format,contents);
end

fclose(fileID)