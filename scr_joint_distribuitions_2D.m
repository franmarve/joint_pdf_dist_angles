% *************************************************************************
% Francisco J. Martin-Vega, fjmv@uma.es
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% This script generates a rectangle and computes the joint distribution of
% the azimuth angle and the distance between a random point and an 
% arbitrary point u. 
% *************************************************************************
clear, clc, close all

% Add search path with the project files
% *************************************************************************
[~, oldPath] = addPaths();

%% Parameters
N = 1e5; % Number of points to place
Lx = 10; % Length in x dimension(m)
Ly = 20; % Length in y dimension (m)

% Coordinates of reference point u = (ux, uy)
u = [-4.999, -9.99];

get_joint_dist_2D(N, Lx, Ly, u);