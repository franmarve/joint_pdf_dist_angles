% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin-Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function computes the antenna gain of the serving RF beam, as
% g_t(theta - varphi), being theta the azimuth angle where the antenna gain
% is evaluated and varphi the bearing angle of the physical antenna. 
% *************************************************************************
% REFERENCE:
% [1] Francisco J. Martin-Vega, Gerardo Gomez, David Morales-Jimenez, 
% F. Javier Lopez-Martinez and Mari Carmen Aguayo-Torres, "Joint 
% Distribution of Distance and Angles in Finite Wireless Networks", acepted
% for publication in IEEE Transactions on Vehicular Technology, 2023.

function g = gServ(theta, p, sP)

% Apply rotation of the antenna array in azimuth direction
theta = theta - sP.varphi;

% Convert angles to the interval ranging from -pi to pi
theta = angle_mpi_to_pi(theta);

if strcmp(p.antennaStr, 'singleHorizontal3GPP') 
    % Antenna gain as per [1], eq. (38) 
    g_dB = sP.gMax_dBi - min( min(12*(theta/sP.HPBW_rad).^2, sP.AMax_dB), ...
        sP.AMax_dB );
    g = db2pow(g_dB);

else
    error('gServ: Unkown p.antennaStr value');
end