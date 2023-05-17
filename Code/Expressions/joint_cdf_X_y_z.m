% *************************************************************************
% Francisco J. Martin-Vega, frmvega@gmail.com
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function estimates the cdf of the RV X for a given value 
% Y. F_{x,y0,z0} = Pr(X <=x, Y <= y0, Z <= <z0). 
% *************************************************************************

function cdfXyz = joint_cdf_X_y_z(X, Y, Z, xVct, y0, z0)

nX = length(X);
cdfXyz = zeros(length(xVct), 1);
for iX = 1:length(xVct)
    cdfXyz(iX) = sum(X <= xVct(iX) & Y <= y0 & Z <= z0)/nX;
end