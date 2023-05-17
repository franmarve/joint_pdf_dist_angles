% *************************************************************************
% Francisco J. Martin-Vega, fjmv@uma.es
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% This script generates some graphs included in Fig. 8 (a) of [1]
% *************************************************************************
% REFERENCE:
% [1] Francisco J. Martin-Vega, Gerardo Gomez, David Morales-Jimenez, 
% F. Javier Lopez-Martinez and Mari Carmen Aguayo-Torres, "Joint 
% Distribution of Distance and Angles in Finite Wireless Networks", acepted
% for publication in IEEE Transactions on Vehicular Technology, 2023.

clear, clc, close all

% Add search path with the project files
% *************************************************************************
[~, oldPath] = addPaths();

%% Parameters
N = 1e5; % Number of points to place
Lx = 200; % Length in x dimension(m)
Ly = 25; % Length in y dimension (m)

% Coordinates of reference point u = (ux, uy)
u = [0, 0];

get_joint_dist_2D(N, Lx, Ly, u);