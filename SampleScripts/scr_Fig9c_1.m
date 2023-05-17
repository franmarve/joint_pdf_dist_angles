% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega 
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This script generates some graphs included in Fig. 9 (c) of [1]
% *************************************************************************
% REFERENCE:
% [1] Francisco J. Martin-Vega, Gerardo Gomez, David Morales-Jimenez, 
% F. Javier Lopez-Martinez and Mari Carmen Aguayo-Torres, "Joint 
% Distribution of Distance and Angles in Finite Wireless Networks", acepted
% for publication in IEEE Transactions on Vehicular Technology, 2023.

clc; clear; 

% Add search path with the project files
% *************************************************************************
[~, oldPath] = addPaths();

% LOAD DEFAULT PARAMETERS
% *************************************************************************
defaultParms;

% CUSTOMIZE PARAMETERS
% *************************************************************************
% *) Result folder
p.resultFolder = 'scr_Fig9c_1_Lx_Ly_10_10_u_0tom4_m4';

% *) Number of elements of the vectorial parameters. Number of points of
% graphs showing key performance metrics.
p.nVal = 20;

% *) Number of subframes to be simulated at each of the nVal simulations
p.nReal = 1e1;

% *) x-label for the performance figures
p.xLabel = 'u_x';

% *) x-axis for the performance figures
p.xVct = num2cell(sort(linspace(-4, 0, p.nVal), 'descend'));

% SPECIFIC PARAMETERS
% *************************************************************************
p.antennaStr = 'singleHorizontal3GPP';
vP.Lx = 10; % Length in x dimension(m)
vP.Ly = 10; % Length in y dimension (m)
vP.fcHz = 2e9;

% mehods: 'avoidDiscont', '2dint'
p.ccdfIntmethod = 'avoidDiscont';

% Coordinates of reference point u = (ux, uy)
vP.ux = cell2mat(p.xVct);
vP.uy = -4;

% Transmit power density dBm/Hz
vP.rhot_dBmpHz = -95; 

% Noise power spectral density dBm/Hz
vP.N0_dBmpHz = -165;

% Vector of SNR threshold values
p.tVctdB = linspace(-30, +30, 18);

% MAIN FUNCTION
% *************************************************************************
getNumResults(p, vP);

% Restore search paths
% *************************************************************************
path(oldPath);
