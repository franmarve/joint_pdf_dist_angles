% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function computes dependant variables
% *************************************************************************

function sP = dependantVar(p, sP)

% Compute location of the reference access point/base station
sP.u = [sP.ux, sP.uy];

% Sanity check: arbritary location must be within the open set R(u)
if abs(sP.u(1)) >= sP.Lx/2 || abs(sP.u(2)) >= sP.Ly/2
    error(['dependantVar: The arbritary location, u, must lie '....
        'within the open set R(o)']);
end

% Compute the PLE and PLS
sP = pathLoss(p, sP);

sP.straightforward_imp = p.straightforward_imp;

% Variable conversion
sP.tVct = db2pow(p.tVctdB);
sP.rhot = db2pow(sP.rhot_dBmpHz - 30);
sP.N0 = db2pow(sP.N0_dBmpHz - 30);