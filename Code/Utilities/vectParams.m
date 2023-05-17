% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function fixes the length of the vectorial params
% *************************************************************************
% INPUTS:
% *) p: Struct with the statics values. It is used by the eval function to
% get the number of simulations p.nVal
% *) vP: Struct with the vectorial parameters. These are those parameters
% that might vary between a simulation realted to iVal=k, and a different
% simulation related to iVal=q, with q!=k. 
% *************************************************************************

function vP = vectParams(p, vP)

vPFields = fieldnames(vP);
nFields = length(vPFields);

for iField = 1:nFields
    eval(['x = ' 'vP.' vPFields{iField} ';'])
    if length(x) == 1
        eval(['vP.' vPFields{iField} ' = repmat(x, 1, p.nVal);'])
    end
end


