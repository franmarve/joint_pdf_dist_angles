% *************************************************************************
% Francisco J. Martin-Vega, fjmv@uma.es
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function represents a conditional assigment. The function returns an 
% x if the condition C is true and it returns y otherwise
% *************************************************************************
function z = cA(x, y, C)

% Sanity check
if ~(size(x, 1) == size(C, 1) && size(x, 2) == size(C, 2) ...
        && size(x, 1) == size(y, 1) && size(x, 2) == size(y, 2))
    error('cA: Input matrices must have the same length');
end

% Initialization
z = zeros(size(x));

% Conditional assigment
z(C) = x(C);
z(~C) = y(~C);