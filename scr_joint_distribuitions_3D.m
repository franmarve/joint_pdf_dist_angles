% *************************************************************************
% Francisco J. Martin-Vega, fjmv@uma.es
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% This script generates a rectangle and computes the joint distribution of
% angles and distances of random points towards an arbitrary point u. 
% *************************************************************************
clear, clc, close all

% Add search path with the project files
% *************************************************************************
[~, oldPath] = addPaths();

%% Parameters
N = 1e5; % Number of points to place. Default 1e5
Lx = 200; % Length in x dimension(m)
Ly = 100; % Length in y dimension (m)

% Coordinates of reference point u = (ux, uy, uz)
u = [30, 25, 10];

% Heigh of antennas at reciver side (random locations)
vz = 150;

get_joint_dist_3D(N, Lx, Ly, u, vz);