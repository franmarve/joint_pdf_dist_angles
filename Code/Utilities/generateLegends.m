% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function returns the string of the legend of figures that compares
% theoretical vs simulation results
% *************************************************************************

function strLegend = generateLegends(baseStr, p)

strLegend = cell(2*p.nVal, 1);
nTypes = length(baseStr);
for j = 1:nTypes
    for iVal = 1:p.nVal
        if isnumeric(p.xVct{iVal})
            str_legend_iVal = num2str(p.xVct{iVal});
        
        elseif ischar(p.xVct{iVal})
            str_legend_iVal = p.xVct{iVal};
        
        else
            error('iValLegend: p.xVct(iVal) does not have a proper type');
        end
        strLegend{iVal + (j-1)*p.nVal} = ...
            sprintf('%s %s = %s', baseStr{j}, p.xLabel, str_legend_iVal);
    end
end