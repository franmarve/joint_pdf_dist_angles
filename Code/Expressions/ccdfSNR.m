% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin-Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function computes the analytical ccdf of the SNR for a vector of
% threshold according to [1]
% *************************************************************************
% REFERENCE:
% [1] Francisco J. Martin-Vega, Gerardo Gomez, David Morales-Jimenez, 
% F. Javier Lopez-Martinez and Mari Carmen Aguayo-Torres, "Joint 
% Distribution of Distance and Angles in Finite Wireless Networks", acepted
% for publication in IEEE Transactions on Vehicular Technology, 2023.

function [pSNR, percSNRdB] = ccdfSNR(p, sP, method)

if nargin < 3
    method = avoidDiscont;
end

% Anonymous functions
dist = @(x) sqrt(x(:,1).^2 + x(:,2).^2); 
[Chil, Chig, gfun, qfun, hi] = anonymousFunctions(sP.Lx, sP.Ly, sP.u);

% vertex of the rectangle
sP.hxp = sP.Lx/2 - sP.u(1);  % Latex notation: h_x^{+}
sP.hxm = -sP.Lx/2 - sP.u(1); % Latex notation: h_x^{-}
if sP.hxm == 0, sP.hxm = -0; end
sP.hyp = sP.Ly/2 - sP.u(2);  % Latex notation: h_y^{+}
sP.hym = -sP.Ly/2 - sP.u(2); % Latex notation: h_y^{-}
if sP.hym == 0, sP.hym = -0; end

V(:, 1) = [sP.hxm, sP.hxp,  sP.hxp, sP.hxm];
V(:, 2) = [sP.hym, sP.hym,  sP.hyp, sP.hyp];

% Maximum distance
rMax = max(dist(V)) + 1;

% Compute the CCDF of the SNR
pSNR = arrayfun(@(x) ccdfSNRfun(x, hi, Chil, Chig, p, sP, rMax, ...
        method), sP.tVct);

% Compute the percentile
percSNRdB = pow2db(fzero(@(x) 1 - ccdfSNRfun(x, hi, Chil, Chig, p, sP, ...
        rMax, method) - p.percentile, [db2pow(-30), db2pow(+30)]));
