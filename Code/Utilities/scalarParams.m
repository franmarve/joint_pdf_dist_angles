% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function fixes the length of the vectorial params
% *************************************************************************

function sP = scalarParams(vP, iVal)

vPFields = fieldnames(vP);

nFields = length(vPFields);

for iField = 1:nFields
    eval(['sP.' vPFields{iField} '= vP.' vPFields{iField} '(iVal);'])
end


