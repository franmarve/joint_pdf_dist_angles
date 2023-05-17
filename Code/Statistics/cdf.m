% *************************************************************************
% Francisco J. Martin-Vega, frmvega@gmail.com
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function estimates the cdf of the RV X. 
% *************************************************************************

function cdfX = cdf(X, xVct)

nX = length(X);
cdfX = zeros(length(xVct), 1);
for iX = 1:length(xVct)
    cdfX(iX) = sum(X <= xVct(iX))/nX;
end