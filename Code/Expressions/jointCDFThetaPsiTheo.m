% *************************************************************************
% Francisco J. Martin-Vega, fjmv@uma.es
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% It obtains the joint CDF of azimuth and zenith angles between a 
% reference point and a set of randomly placed points. 
% *************************************************************************
function outVct = jointCDFThetaPsiTheo(epsilon, hi, theta, psi, Lx, Ly, u, vz)

T1 = double(u(3) >= vz).*(  double(pi/2 <= psi & psi <= pi).*...
    pp(cdfThetaTheo(epsilon, hi, theta, Lx, Ly, u) - ...
    jointCDFTheo((u(3) - vz).*tan(pi-psi), theta, Lx, Ly, u)) + ...
    double(psi <= pi/2).*cdfThetaTheo(epsilon, hi, theta, Lx, Ly, u) );

T2 = double(u(3) < vz).*double(psi <= pi/2).*...
    jointCDFTheo( (vz - u(3)).*tan(psi), theta, Lx, Ly, u);

outVct = T1 + T2;