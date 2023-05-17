% *************************************************************************
% Francisco J. Martin-Vega, frmvega@gmail.com
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function estimates the pdf of the RV X. It is used nBins to estimate
% the histogram. Then with the histogram it is obtained the estimated pdf.
% *************************************************************************

function LTX = laplaceTransformRV(X, sVct)

nVal = length(sVct);
LTX = zeros(nVal, 1);
for iVal = 1:nVal 
    LTX(iVal) = mean(exp(-sVct(iVal)*X));
end

