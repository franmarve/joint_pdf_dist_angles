% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function fill the result struc for iVal
% *************************************************************************    

function [vR, vF, uP] = resultStruc(p, sR, vR, sF, vF, sP, uP, iVal)

% Figures
vF.hf{iVal} = sF.hf;

uP{iVal} = sP;

% Numerical evaluation of empirical CCDF of the SNR
vR.ccdf_simSNR.x{iVal} = p.tVctdB;
vR.ccdf_simSNR.y{iVal} = sR.ccdf_simSNR.y;

% Numerical evaluation of theoretical CCDF of the SNR
vR.ccdf_theorSNR.x{iVal} = p.tVctdB;
vR.ccdf_theorSNR.y{iVal} = sR.theoCcdfSNR;

vR.percSNRdB(iVal) = sR.percSNRdB;

