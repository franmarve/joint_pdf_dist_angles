% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This script performs numerical evaluation of theoretical expressions and
% monte carlo simulations about the CCDF of the 
% *************************************************************************

function getNumResults(p, vP)

vR = []; vF = []; uP = [];

% VECTORIAL PARAMETERS
% *************************************************************************
vP = vectParams(p, vP);

% CREATION OF A RESULT FOLDER
% *************************************************************************
p = createResultsFolder(p);

% MAIN BODY
% *************************************************************************
hwin = waitbar(0, 'Running simulation points...');
count = 0;
nIters = p.nVal;

t0 = clock;

for iVal = 1:p.nVal
    % OBTAINING SCALAR VALUES
    % *********************************************************************
    sP = scalarParams(vP, iVal);

    % COMPUTE DEPENDANT VARIABLES
    % *********************************************************************
    sP = dependantVar(p, sP);
    
    % MONTE CARLO SIMULATION
    % *********************************************************************
    [sR, sF] = ranSimMetrics(p, sP);
    
    % NUMERICAL EVALUATION
    % *********************************************************************
    [theoCcdfSNR, percSNRdB] = ccdfSNR(p, sP, p.ccdfIntmethod);
    [sR.theoCcdfSNR, sR.percSNRdB] = deal(theoCcdfSNR, percSNRdB);
    
    % SAVE CURRENT RESULTS AND FILL RESULTS STRUCT
    % *********************************************************************
    [vR, vF, uP] = resultStruc(p, sR, vR, sF, vF, sP, uP, iVal);

    save([p.folderName '\sR_' num2str(iVal) '.mat'], 'p', 'sR', 'sP')

    count = count + 1;
    waitbar(count/nIters, hwin);
end

% SAVING RESULTS AND KEY PERFORMANCE FIGURES
% *********************************************************************
% Elapsed time
t1 = clock; 
t0min = t0(4)*60 + t0(5);
t1min = t1(4)*60 + t1(5);
vR.eHours = (t1min - t0min)/60;
% Saving results
save([p.folderName '\results.mat'], 'vR', 'p', 'vP', 'uP')

savePerformanceFig(p, vR, vF);

close(hwin);


