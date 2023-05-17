% *************************************************************************
% Francisco J. Martin-Vega, frmvega@gmail.com
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function estimates the ccdf of the RV X. 
% *************************************************************************

function ccdfX = ccdf(X, xVct)

nX = length(X);
ccdfX = zeros(length(xVct), 1);
for iX = 1:length(xVct)
    ccdfX(iX) = sum(X > xVct(iX))/nX;
end