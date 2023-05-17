% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function fill the result struc for iVal
% *************************************************************************

function sP = pathLoss(p, sP)

if strcmp(p.pathLossModel, 'explicit')
    if ~(isfield(sP,'tau') && isfield(sP,'alpha'))
        error('pathLoss: non existent fields tau and alpha');
    end

elseif strcmp(p.pathLossModel, 'InH-Office LOS')
    tau = @(fc) ((fc/1e9)^2*10^3.24)^(1/1.73);
    sP.tau = tau(sP.fcHz);
    sP.alpha = 1.73;
    
else
    error('pathLoss: Unkown path loss model (p.pathLossModel)')
    
end