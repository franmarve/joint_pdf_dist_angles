% *************************************************************************
% Francisco J. Martin-Vega, frmvega@gmail.com
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function estimates the cdf of the RV X for a given value 
% Y. F_{x,y0} = Pr(X <=x, Y <= y0). 
% *************************************************************************

function cdfXy = joint_cdf_X_y(X, Y, xVct, y)

nX = length(X);
cdfXy = zeros(length(xVct), 1);
for iX = 1:length(xVct)
    cdfXy(iX) = sum(X <= xVct(iX) & Y <= y)/nX;
end