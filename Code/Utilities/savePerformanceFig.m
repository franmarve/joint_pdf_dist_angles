% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function save figures related to key performance metrics
% *************************************************************************

function savePerformanceFig(p, vR, vF)

% Save vF figures
% Figure names
if p.plotRealizations
    strFigNames = {'polarPlotBeams', 'polarPlotServBeam', 'plotBeams', ...
        'pdfTheta', 'cdfTheta', 'Locations'};
    nF = length(vF.hf{1});
    for iF = 1:nF
        for iVal = 1:p.nVal
            saveas(vF.hf{iVal}(iF), [p.folderName '\' strFigNames{iF} '_iVal' ...
                num2str(iVal) '.fig'])
            close(vF.hf{iVal}(iF));
        end
    end
end

% # This is sample code. Remove with your own code
hf = figure;
% Plot theoretical results
for iVal = 1:p.nVal
    plot(vR.ccdf_theorSNR.x{iVal}, vR.ccdf_theorSNR.y{iVal});
    if iVal == 1, hold on; end
end
ylabel('F_{SNR}(t)'); xlabel('t (dB)'); grid on;

% Plot simulation results with markers
for iVal = 1:p.nVal
    plot(vR.ccdf_simSNR.x{iVal}, vR.ccdf_simSNR.y{iVal}, 'o');
end
strLegend = generateLegends({'Theor', 'Sim'}, p);
legend(strLegend);
saveas(hf, [p.folderName '\ccdfSNR.fig'])
close(hf)

hf = figure;
if isnumeric(p.xVct{1})
    plot(cell2mat(p.xVct), vR.percSNRdB, '-x');
elseif ischar(p.xVct{1})
    plot(1:p.nVal, vR.percSNRdB);
else
    error('savePerformanceFig: incorrect type of p.xVct');
end
title(sprintf('Percentile %s', num2str(p.percentile)))
ylabel('SNR_p(dB)'); xlabel(p.xLabel); grid on;
saveas(hf, [p.folderName '\percentileSNR.fig'])
close(hf)




