% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega 
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION: This script contains the default parameters.
% *************************************************************************

%% STATIC PARAMS
% *************************************************************************
% *) Result folder
p.resultFolder = 'test_case_1';

% *) Number of elements of the vectorial parameters. Number of points of
% graphs showing key performance metrics.
p.nVal = 10;

% *) Number of subframes to be simulated at each of the nVal simulations
p.nReal = 1e2;

% *) x-label for the performance figures. It represents the name of the
% parameter that is varying in the drawn figures
% 
p.xLabel = '';

% *) x-axis for the performance figures
p.xVct = {};

% Vector of SNR threshold values
p.tVctdB = linspace(-30, +30, 18);

% Probability value to compute percentiles of the CDF. Interval [0, 1].
p.percentile = 0.1;

% Number of UEs to place
p.nUE = 1e3; 

% Boolean variable to indicate whether to plot debug figures, i.e., 
% the beams,  UE locations, etc.
p.plotRealizations = true;

% Antenna array structure
p.antennaStr = 'ULA'; % 'ULA', 'singleHorizontal3GPP', 'sectorized'

% Codebook design method
% 'equidistant': The beams are equidistant in the angular domain between 0
%   and 2pi
% 'equidensity': The beams serve the same user density. The steering
% directions of each beam is computed based on the marginal CDF of the user
% density in the angle domain. 
% Note: Only used for p.antennaStr == 'ULA'
p.method = 'equidistant'; 

% Path loss model: 
% - 'explicit' It uses the path loss slope (PLS), vP.tau, and path loss
%    exponent (PLE), vP.alpha to compute the path loss. 
% - 'InH-Office LOS' Computes the PLS, sP.tau, and PLE, sP.alpha
%   according to the InH-Office model (3GPP TR 38.913 v17.0.0) for a given 
%   carrier frequency, vP.fc. 
p.pathLossModel = 'InH-Office LOS';

% Implementation of serving beam selection function
% Note: Only used for p.antennaStr == 'ULA'
p.straightforward_imp = true;

% Method to integrate the CCDF of the SNR
% mehods: 'avoidDiscont', '2dint'
p.ccdfIntmethod = '2dint';

%% VECTOR PARAMS
% *************************************************************************
% Shape of the region
vP.Lx = 40; % Length in x dimension(m)
vP.Ly = 25; % Length in y dimension (m)

% Coordinates of reference point u = (ux, uy)
vP.ux = 0;
vP.uy = 0;

% Orientation azimith angle of the antenna array
vP.varphi = 0;

% Transmit power density dBm/Hz
vP.rhot_dBmpHz = -105; 

% Noise power spectral density dBm/Hz
vP.N0_dBmpHz = -165;

% Verification of the sinc radiation pattern
thetaVct = linspace(0, 2*pi, 1e3);

% Carrier frequency in Hz
vP.fcHz = 4e9;

% Path loss slope
% Note: It is only used if p.pathLossModel == explicit
vP.tau = 1;

% Path loss exponent
% Note: It is only used if p.pathLossModel == explicit
vP.alpha = 3.8;

% PARAMS for p.antennaStr == 'ULA'
% *************************************************************************
% Number of transmit antennas
% Note: Only used for p.antennaStr == 'ULA'
vP.Nt = 8;

% Distance between antenna elements normalized by the wavelength
% Note: Only used for p.antennaStr == 'ULA'
vP.dlambda = 0.25;

% Number of beams
% Note: Only used for p.antennaStr == 'ULA'
vP.Nb = 15;

% Max and min steering angles
% Note: Only used for p.antennaStr == 'ULA'
vP.phiMax = 2*pi;
vP.phiMin = 0;

% PARAMS for p.antennaStr == 'singleHorizontal3GPP'
% *************************************************************************
% Maximum directional gain of an antenna element
vP.gMax_dBi = 8;

% Half power beamwidth (HPBW) in radians
vP.HPBW_rad = 65*pi/180;

% side-lobe attenuation
vP.AMax_dB = 30;
