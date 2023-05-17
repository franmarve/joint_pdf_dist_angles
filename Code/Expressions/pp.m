% *************************************************************************
% Francisco J. Martin-Vega, fjmv@uma.es
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% This functions returns the positive part of a vector or 0 otherwise 
% *************************************************************************
function y = pp(x)

% Initialization
y = zeros(size(x));

% Conditional assigment
y(x > 0) = x(x > 0);