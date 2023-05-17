% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega
% *************************************************************************
% GROUP: Room 1.3.7., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This script calls the main function to obtain the CCDF of the SNR and the 
% p-th percentile of the SNR. The scenario considers a base station located
% at an arbritary location, u, and a randomly deplyed user, placed within Add
% finite rectangular region.  
% *************************************************************************
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
p.resultFolder = 'Ly20_Lx10_u_3_m2_ccdfSNR';

% *) Number of elements of the vectorial parameters. Number of points of
% graphs showing key performance metrics.
p.nVal = 1;

% *) Number of subframes to be simulated at each of the nVal simulations
p.nReal = 1e5;

% *) x-label for the performance figures
p.xLabel = 'L_x';

% *) x-axis for the performance figures
p.xVct = num2cell([10]);

% SPECIFIC PARAMETERS
% *************************************************************************
p.antennaStr = 'singleHorizontal3GPP';
vP.Ly = 20; % Length in y dimension (m)
vP.Lx = cell2mat(p.xVct);
vP.fcHz = 2e9;

% mehods: 'avoidDiscont', '2dint', 'degeneratedLy', 'approximationLy'
p.ccdfIntmethod = 'avoidDiscont';

% Coordinates of reference point u = (ux, uy)
vP.ux = 3;
vP.uy = -2;

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
