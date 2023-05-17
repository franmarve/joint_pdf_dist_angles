% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function obtains the key performance metrics of a given set of
% system parameters stored in the struct sP and p
% *************************************************************************

function [sR, sF] = ranSimMetrics(p, sP)

% Struct to store metrics 
sR = [];

% struct of figures
sF = [];

% Anonymus function definitions
dist = @(x) sqrt(x(:,1).^2 + x(:,2).^2); 

% REPRESENTATION OF ANALOG CODEBOOK
% *************************************************************************
nfig = 0;
if p.plotRealizations && strcmp(p.antennaStr, 'ULA')
    thetaVct = linspace(0, 2*pi, 1e3);
     nfig = nfig + 1; hf(nfig) = figure;
    for n = 1:sP.Nb
        polarplot(thetaVct, gULA(thetaVct, sP, n)); hold on
        polarplot([sP.phin(n)], linspace(0,1,10), '.k');
        polarplot([sP.etan(n)], linspace(0,1,10), 'or');
    end

    nfig = nfig + 1; hf(nfig) = figure;
    polarplot(thetaVct, gServ(thetaVct, p, sP))

    nfig = nfig + 1; hf(nfig) = figure;
    for n = 1:sP.Nb
        plot(thetaVct/pi, gULA(thetaVct, sP, n)); hold on
        plot([sP.phin(n)/pi, sP.phin(n)/pi], [0, 1], '--k')
    end
    
    % Plot marginal PDF of the azimuth angle
    thetaVct = linspace(0, 2*pi, 100);
    nfig = nfig + 1; hf(nfig) = figure;
    pdfThetaVct = pdfThetaTheo(thetaVct, sP.Lx, sP.Ly, sP.u);
    plot(thetaVct/pi, pdfThetaVct, 'b'); hold on;
    for n = 1:sP.Nb
        plot([sP.phin(n)/pi, sP.phin(n)/pi], [0, max(pdfThetaVct)], '--k')
    end
    grid on; ylabel('f_\Theta(\theta)'); xlabel('\theta/\pi')
    
    nfig = nfig + 1; hf(nfig) = figure;
    % Plot marginal CDF of theta
    cdfThetaVct = cdfThetaTheo(thetaVct, sP.Lx, sP.Ly, sP.u);
    plot(thetaVct/pi, cdfThetaVct , 'b'); hold on;
    for n = 1:sP.Nb
        plot([sP.phin(n)/pi, sP.phin(n)/pi], [0, 1], '--k')
    end
    grid on; ylabel('F_\Theta(\theta)'); xlabel('\theta/\pi')
    
    sF.hf = hf;
end

% REALIZATIONS OF KEY PERFORMANCE METRICS
% *************************************************************************
hwin = waitbar(0,sprintf('Running simulation of %d slots...', p.nReal));

% Place p.nReal UEs randomly within the rectangular region
[~, Theta, ~, R, V, hf] = randomCoord(p.plotRealizations, sP.Lx, sP.Ly, ...
    p.nReal, sP.u);
if p.plotRealizations
    % Plot radiation pattern
    thetaVct = linspace(0, 2*pi, 1e3);
    gdBVct = pow2db(gServ(thetaVct, p, sP));

    % Maximum and miinimum values for graphical representation
    gMax = max(gdBVct); gMin = min(gdBVct);

    % Normalize gdBVct for graphical representation (gain->rho > 0)
    % Minimum desired gain on the plot: gMin
    gDmin = 0.5;
    % Maximum desired gain on the plot: gMax
    gDmax = 5;
    % Displacement eta0 and normalization factors eta1
    eta0 = gDmin*(gMin - gMax)/(gDmin - gDmax) - gMin;
    eta1 = (gMin - gMax)/(gDmin - gDmax);

    gdBVct = (gdBVct + eta0)/eta1;

    % Convert to cartesian coordinates
    [X, Y] = pol2cart(thetaVct, gdBVct);
    % Gain in dB and normalized
    figure(hf); hold on; plot(X, Y, 'g', 'LineWidth', 2);

    % Gain in linear form
    gVct = gServ(thetaVct, p, sP);

    % Convert to cartesian coordinates
    [X, Y] = pol2cart(thetaVct, gVct);
    plot(X, Y, 'r', 'LineWidth', 2);
    nfig = nfig + 1; sF.hf(nfig) = hf; 
end

% Generate fading realizations
ExpMean = 1;
B = random('exp', ExpMean, p.nReal, 1);

% Compute the SNR for each UE
SNR = (B.*(sP.tau.*R).^-sP.alpha.*sP.rhot.*gServ(Theta, p, sP))./sP.N0;

% for iReal = 1:p.nReal
%     waitbar(iReal/p.nReal, hwin);
% end % end iSpatReal

close(hwin)

% CCDFs, CDFs and PDFs
% *************************************************************************
sR.ccdf_simSNR.x = sP.tVct;
sR.ccdf_simSNR.y = ccdf(SNR, sR.ccdf_simSNR.x);

% MOMENTS
% *************************************************************************
% sR.AvMetric_1 = mean(samples_1);