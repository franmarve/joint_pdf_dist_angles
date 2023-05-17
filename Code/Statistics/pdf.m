% *************************************************************************
% Francisco J. Martin-Vega, frmvega@gmail.com
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function estimates the pdf of the RV X. It is used nBins to estimate
% the histogram. Then with the histogram it is obtained the estimated pdf.
% *************************************************************************

function [pdfX, x] = pdf(X, xVct)

nReal = length(X);
[n, x] = hist(X, xVct);
deltaX = diff(x);
pdfX =  n(1:end-1)./(nReal*deltaX);
x = x(1:end-1);

