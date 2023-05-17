% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin-Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function draws the theoretical and estimated distribuitions of 
% distance, zenith and azimuth angle for the 3D case
% *************************************************************************

function get_joint_dist_3D(N, Lx, Ly, u, vz)

% Sanity check: arbritary location must be within the open set R(u)
if abs(u(1)) >= Lx/2 || abs(u(2)) >= Ly/2
    error(['get_joint_dist_3D: The arbritary location, u, must lie '....
        'within the open set R(o)']);
end

font = 'Times new roman';
fontSize = 16;
cell_latex = {'fontsize',fontSize,'interpreter','latex'}; % 'Fontname','Times new roman',

% PDF estimation: MATLAB ksdensity or histogram based: 'ksdensity', 'hist'
pdfEst = 'hist';

[Chil, Chig, gfun, qfun, hi, epsilon] = anonymousFunctions(Lx, Ly, u);

%% SIMULATION results: Place in polygon
% Euclidean distance
dist = @(x) sqrt(x(:,1).^2 + x(:,2).^2); 
psif = @(rVct, uz, vz) cA(pi - atan(rVct./(uz - vz)), ...
    atan(rVct./(vz - uz)), (uz > vz));

[D, Theta, Psi, R, V, hf] = randomCoord(true, Lx, Ly, N, u, vz, psif);
figure(hf); set(gca,'fontsize',fontSize)

%% Estimate joint CDF and marginals for distance (d), azimuth (theta) 
% and zenith (psi) angles
nValSim = 50;
thetaMax = 2*pi;
zeroPlus = 1e-3;
ThetaVct = linspace(0, thetaMax, nValSim);
rMax = max(dist(V)) + 1;
dMax = sqrt(rMax^2 + (u(3) - vz)^2);
DVct = linspace(abs(u(3) - vz)*0.9, dMax, nValSim);
psiMax = cA(pi, psif(rMax, u(3), vz), (u(3) >= vz));
psiMin = cA(psif(rMax, u(3), vz), zeroPlus, (u(3) >= vz));
PsiVct = linspace(psiMin, psiMax, nValSim);

theta0_per_quad_deg = [40, 60, 70; 120, 135, 150; 210, 225, 240; ...
    300, 315, 330];

nPsi = 2;
[nQuad, nTheta] = size(theta0_per_quad_deg);
cdfR_theta = cell(nQuad, nTheta, nPsi);

% Take for instance two equally spaced
psiMin_iQ = zeros(nQuad, 1);
psiMax_iQ = zeros(nQuad, 1);
psi0_per_quad_deg = zeros(nQuad, nPsi);
for iQ = 1:nQuad
    rMax_iQ = sqrt(V(iQ,1)^2 + V(iQ,2)^2);
    psiMax_iQ(iQ) = cA(pi, psif(rMax_iQ, u(3), vz), (u(3) >= vz));
    psiMin_iQ(iQ) = cA(psif(rMax_iQ, u(3), vz), 0, (u(3) >= vz));
    
    psi0_per_quad_deg(iQ, :) = [(psiMax-psiMin_iQ(iQ))/3 + psiMin_iQ(iQ), ...
    2*(psiMax-psiMin_iQ(iQ))/3 + psiMin_iQ(iQ)];
end

psi0_per_quad_deg = round(180*psi0_per_quad_deg/pi);

for iQ = 1:nQuad
    for iTheta0 = 1:nTheta
        for iPsi0 = 1:nPsi
            cdfR_theta{iQ, iTheta0, iPsi0} = ...
                joint_cdf_X_y_z(D, Theta, Psi, ...
                DVct, theta0_per_quad_deg(iQ, iTheta0)/180*pi, ...
                psi0_per_quad_deg(iQ, iPsi0)/180*pi); 
        end
    end
end

% Estimate joint CDF of angles
data = [Theta, Psi];
stR = gkde2(data);

% Marginal PDF zenith angle
if (strcmp(pdfEst, 'ksdensity'))
    [f_Psi, xi_Psi, bw_Psi] = ksdensity(Psi);
    
elseif (strcmp(pdfEst, 'hist'))

    [f_Psi, xi_Psi] = pdf(Psi, PsiVct);
    
else
    error('pdfEst: Invalid value');
end

% Marginal CDF zenith angle
cdfPsi = cdf(Psi, PsiVct);

%% Plot statistics

hf_Q = cell(nQuad);
% Different colors for different theta_0
color_theta = {'blue', 'red', 'green'};
% Marker only for simulation results
marker_psi = {'x', 'o'};
% Line style none for sim. - theor. psi_0(1) -- theor. psi_0(2)
lineStyle_sim_theo = {'none', '-', '--'}; 

legend_str = cell(nPsi, nTheta, nQuad);
for iQ = 1:nQuad
    hf_Q{iQ} = figure('DefaultTextFontName', font, 'DefaultAxesFontName', font); 
    for iTheta0 = 1:nTheta
        for iPsi0 = 1:nPsi
            plot(DVct, cdfR_theta{iQ, iTheta0, iPsi0}, 'LineStyle', ...
                lineStyle_sim_theo{1}, 'Color', color_theta{iTheta0}, ...
                'Marker', marker_psi{iPsi0})         
            if (iTheta0 == 1 && iPsi0 == 1), hold on; end
        end
    end
    
    title(sprintf(['Joint CDF for \\theta = (%s) and \\psi = (%s) '...
        'degrees. Quadrant: %d'], ...
        num2str(theta0_per_quad_deg(iQ,:)),...
        num2str(psi0_per_quad_deg(iQ,:)),iQ));
    xlabel('$d (m)$', cell_latex{:}); grid on;
    ylabel('$F_{D, \Theta, \Psi}(D, \theta = \theta_0, \Psi = \psi_0)$', cell_latex{:});
    set(gca,'fontsize',fontSize)
end 

% Plot joint CDF distance, angle
hf_theta_psi_PDF = figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
surf(stR.x/pi, stR.y/pi, stR.pdf);
xlabel('$\theta/\pi$', cell_latex{:}); ylabel('$\psi/\pi$', cell_latex{:}); 
zlabel('$f_{\Theta, \Psi}(\theta, \psi)$', cell_latex{:}); grid on;
title(['Estimated joint PDF of distance azimuth and zenith angles '...
    '$(\Theta, \Psi)$'], cell_latex{:})
set(gca,'fontsize',fontSize)

hf_theta_psi_CDF = figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
surf(stR.x/pi, stR.y/pi, stR.cdf);
xlabel('$\theta/\pi$', cell_latex{:}); 
ylabel('$\psi/\pi$', cell_latex{:}); 
zlabel('$F_{\Theta, \Psi}(\theta, \psi)$', cell_latex{:}); grid on;
title(['Estimated joint CDF of distance azimuth and zenith angles '...
    '$(\Theta, \Psi)$'], cell_latex{:})
set(gca,'fontsize',fontSize)

%% Add Theoretical results
% Joint CDF for \Theta = \theta_0 and \Psi = \psi_0
nValTheo = 200;
ThetaTheoVct = linspace(zeroPlus, thetaMax, nValTheo);
DTheoVct = linspace(abs(u(3)-vz)*0.9, dMax, nValTheo);

cdfR_theta_Theo = cell(nQuad, nTheta, nPsi);

for iQ = 1:nQuad
    figure(hf_Q{iQ}); 
    for iTheta0 = 1:nTheta
        for iPsi0 = 1:nPsi
            cdfR_theta_Theo{iQ, iTheta0, iPsi0} = ...
                jointCDF3DTheo(DTheoVct, ones(size(ThetaTheoVct))...
                *theta0_per_quad_deg(iQ, iTheta0)/180*pi, ...
                ones(size(ThetaTheoVct))*...
                psi0_per_quad_deg(iQ, iPsi0)/180*pi, u, vz, Lx, Ly);
            plot(DTheoVct, cdfR_theta_Theo{iQ, iTheta0, iPsi0}, ...
                'LineStyle', lineStyle_sim_theo{iPsi0+1}, ...
                'Color', color_theta{iTheta0}); 
            if iTheta0 == 1, hold on; end
            legend_str{iPsi0, iTheta0, iQ} = sprintf(...
                '\\theta = %dº, \\psi = %dº', ...
                (theta0_per_quad_deg(iQ, iTheta0)), ...
                (psi0_per_quad_deg(iQ, iPsi0)));

        end
    end
    legend_str_iQ = legend_str(:,:,iQ);
    legend(legend_str_iQ(:))
    set(gca,'fontsize',fontSize)
end 

% Joint CDF and PDF of zenith and azimuth angles
ThetaTheoVct = linspace(zeroPlus, thetaMax, nValTheo);
PsiTheoVct = linspace(psiMin, psiMax, nValTheo);
[TTheta, PPsi] = meshgrid(ThetaTheoVct, PsiTheoVct);

jointCDFAngTheoVct = jointCDFThetaPsiTheo(epsilon, hi, TTheta, PPsi, Lx, Ly, u, vz);
hf_jointCDFAng = figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
surf(TTheta/pi, PPsi/pi, jointCDFAngTheoVct); 
xlabel('$\theta/\pi$', cell_latex{:}); ylabel('$\psi/\pi$', cell_latex{:}); 
zlabel('$F_{\Theta, \Psi}(\theta, \psi)$', cell_latex{:}); grid on;
title(['Theoretical joint CDF of distance azimuth and zenith angles '...
    '$(\Theta, \Psi)$'], cell_latex{:})
set(gca,'fontsize',fontSize)

jointPDFAngTheoVct = jointPDFThetaPsiTheo(Chil, Chig, TTheta, PPsi, Lx, Ly, u, vz);
hf_jointPDFAng = figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
surf(TTheta/pi, PPsi/pi, jointPDFAngTheoVct); 
xlabel('$\theta/\pi$', cell_latex{:}); ylabel('$\psi/\pi$', cell_latex{:}); 
zlabel('$f_{\Theta, \Psi}(\theta, \psi)$', cell_latex{:}); grid on;
title(['Theoretical joint PDF of distance azimuth and zenith angles '...
    '$(\Theta, \Psi)$'], cell_latex{:})
set(gca,'fontsize',fontSize)