% *************************************************************************
% Francisco J. Martin-Vega, fjmv@uma.es
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function represents a conditional function. 
% *************************************************************************
function y = F(f, a, b)

% Resize in case that one of the inputs is an scalar
if numel(a) == 1 && numel(b) > 1
    a = repmat(a, size(b));
elseif numel(b) == 1 && numel(a) > 1
    b = repmat(b, size(a));
end

% Sanity check
if ~(size(a, 1) == size(b, 1) && size(a, 2) == size(b, 2))
    error('F: Input matrices must have the same length');
end

% Initialization
y = zeros(size(a));

% Conditional assigment
bga = b > a;
y(bga) = f(b(bga)) - f(a(bga));

% Equivalent code: 
% N = numel(a);
% for n = 1:N
%     if b(n) > a(n)
%         y(n) = f(b(n)) - f(a(n));
%     end
% end