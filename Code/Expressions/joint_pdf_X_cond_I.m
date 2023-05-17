% *************************************************************************
% Francisco J. Martin-Vega, frmvega@gmail.com
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function estimates the pdf of the RV X conditioned on event I, 
% which is defined as the RV Y being within the interval I = [y_min, y_max]
% *************************************************************************

function [pdfX_cond_I, x] = joint_pdf_X_cond_I(X, Y, xVct, y_min, y_max)

ind_in_I = (Y >= y_min) & (Y <= y_max);
X = X(ind_in_I);

if sum(ind_in_I) == 0
    pdfX_cond_I = zeros(size(xVct));
    x = xVct;
else
    [pdfX_cond_I, x] = pdf(X, xVct);
end

% Sanity check
ind_NaN = isnan(pdfX_cond_I);

if sum(ind_NaN) > 1
    warning('joint_pdf_X_cond_I: NaN');
end