% *************************************************************************
% Francisco J. Martin-Vega, fjmv@uma.es
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% It obtains the joint CDF of distance and angle between a reference point
% and a set of randomly placed points with uniform distribution. 
% *************************************************************************

function outVct = jointCDF3DTheo(d, theta, psi, u, vz, Lx, Ly)

E11 = double(pi/2 <= psi & psi <= pi).*...
    pp(jointCDFTheo(sqrt(d.^2 - (u(3) - vz).^2), theta, Lx, Ly, u) - ...
    jointCDFTheo((u(3) - vz).*tan(pi-psi), theta, Lx, Ly, u));

E12 = double(psi < pi/2).*...
    jointCDFTheo(sqrt(d.^2 - (u(3) - vz).^2), theta, Lx, Ly, u);

E1 = double(d > abs(u(3) - vz)).*double(u(3) >= vz).*(E11 + E12);

arg_CDF_E2 = min(sqrt(d.^2 - (u(3) - vz).^2), (vz - u(3)).*tan(psi));
E2 = double(d > abs(u(3) - vz)).*double(u(3) < vz).*...
    jointCDFTheo(arg_CDF_E2, theta, Lx, Ly, u);

outVct = E1 + E2;
