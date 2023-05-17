% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega
% *************************************************************************
% GROUP: Room 1.3.7., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This script generates some graphs included in Fig. 10 of [1]
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
p.resultFolder = 'scr_Fig10_1_Lx_Ly_0p01To10_100_u_0_0_xi_pid2';

% *) Number of elements of the vectorial parameters. Number of points of
% graphs showing key performance metrics.
p.nVal = 4;

% *) Number of subframes to be simulated at each of the nVal simulations
p.nReal = 1e5;

% *) x-label for the performance figures
p.xLabel = 'L_x';

% *) x-axis for the performance figures
p.xVct = num2cell([0.01, 0.1, 1, 10]);

% SPECIFIC PARAMETERS
% *************************************************************************
p.antennaStr = 'singleHorizontal3GPP';
vP.Ly = 100; % Length in y dimension (m)
vP.Lx = cell2mat(p.xVct);
vP.fcHz = 2e9;

% mehods: 'avoidDiscont', '2dint', 'degeneratedLy', 'approximationLy'
p.ccdfIntmethod = 'approximationLy';

% Coordinates of reference point u = (ux, uy)
vP.ux = 0;
vP.uy = 0;

vP.varphi = pi/2;

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
