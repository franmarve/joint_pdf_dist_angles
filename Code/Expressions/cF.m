% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin-Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% Conditional definite integral. This function computes the definite
% integral of the function, fun, being its integration limits, a(x), and 
% b(x), functions of the independent variable x, if and only if a(x) <
% b(x). The function returns 0 otherwise. 
% *************************************************************************

function z = cF(fun, x, a, b)

aa = a(x);
bb = b(x);

% Resize in case that one of the inputs is an scalar
if numel(aa) == 1 && numel(bb) > 1
    aa = repmat(aa, size(bb));
elseif numel(bb) == 1 && numel(aa) > 1
    bb = repmat(bb, size(aa));
end

% Sanity check
if ~(size(aa, 1) == size(bb, 1) && size(aa, 2) == size(bb, 2))
    error('cF: Input matrices must have the same length');
end

% Initialization
z = zeros(size(aa));

% Conditional assigment
bga = bb > aa;
zfun = @(xx) integral(@(y) fun(xx,y), a(xx), b(xx));

if (sum(bga) > 1)
    z(bga) = arrayfun(zfun, x(bga));
end
