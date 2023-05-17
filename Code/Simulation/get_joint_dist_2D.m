% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin-Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function draws the theoretical and estimated distribuitions of 
% distance and azimuth angle for the 2D case
% *************************************************************************

function get_joint_dist_2D(N, Lx, Ly, u)

% Sanity check: arbritary location must be within the open set R(u)
if abs(u(1)) >= Lx/2 || abs(u(2)) >= Ly/2
    error(['get_joint_dist_2D: The arbritary location, u, must lie '....
        'within the open set R(o)']);
end

font = 'Times new roman';
fontSize = 16;
cell_latex = {'fontsize',fontSize,'interpreter','latex'}; % 'Fontname','Times new roman',

% PDF estimation: MATLAB ksdensity or histogram based: 'ksdensity', 'hist'
pdfEst = 'hist';

% Debug booleans to plot different figures
debug.region = true;
debug.jointPDF = true;
debug.jointCDF = true;
debug.conditionedJointCDF = true;
debug.conditionedJointPDF = true;
debug.marginalCDFTheta = true;
debug.marginalPDFTheta = true;
debug.marginalPDF_R = true;

version_R_cond_O = 2;

% Angles for conditional distribuitions
angles_per_quadrant_deg = [40, 60, 70; 120, 135, 150; 210, 225, 240; ...
    300, 315, 330];

limit_angles_per_quadrant_deg = [0, 90; 90, 180; 180, 270; 270, 360];

[Chil, Chig, gfun, qfun, hi, epsilon] = anonymousFunctions(Lx, Ly, u);

%% SIMULATION results: Place in polygon
% Euclidean distance
dist = @(x) sqrt(x(:,1).^2 + x(:,2).^2); 

[~, Theta, ~, R, V, hf] = randomCoord(true, Lx, Ly, N, u);
figure(hf); set(gca,'fontsize',fontSize)
if ~debug.region
    close(hf);
end

%% Estimate joint PDF and marginals for distance and angle
data = [R, Theta];
stR = gkde2(data);

thetaMax = 2*pi;
ThetaVct = linspace(0, thetaMax, 50);
rMax = max(dist(V)) + 1;
RVct = linspace(0, rMax, 50);

% Marginal PDF angle
if (strcmp(pdfEst, 'ksdensity'))
    [f_Thetha, xi_Thetha, bw_Thetha] = ksdensity(Theta);
    
elseif (strcmp(pdfEst, 'hist'))

    [f_Thetha, xi_Thetha] = pdf(Theta, ThetaVct);
    
else
    error('pdfEst: Invalid value');
end

% CDF angle
cdfTheta = cdf(Theta, ThetaVct);

[nQuad, nAng] = size(angles_per_quadrant_deg);
cdfR_theta = cell(nQuad, nAng);

for iQ = 1:nQuad
    for iAng = 1:nAng
        cdfR_theta{iQ, iAng} = joint_cdf_X_y(R, Theta, RVct, ...
            angles_per_quadrant_deg(iQ, iAng)/180*pi); 
    end
end

% Conditional PDF
pdfR_theta = cell(nQuad, 1);
for iQ = 1:nQuad
    theta_min = min(limit_angles_per_quadrant_deg(iQ, :))/180*pi;
    theta_max = max(limit_angles_per_quadrant_deg(iQ, :))/180*pi;

    [pdfX_cond_I_iQ, x] = joint_pdf_X_cond_I(R, Theta, RVct, theta_min, ...
        theta_max);

    pdfR_theta{iQ}.pdf = pdfX_cond_I_iQ; 
    pdfR_theta{iQ}.x = x; 
end

% Marginal PDF R
if debug.marginalPDF_R
    [pdfR, r] = pdf(R, RVct);
    
    pdf_R.pdf = pdfR; 
    pdf_R.r = r; 
end

% Sanity CHECKs: compute the definite double integral of the bivariate PDF
% It must be exactly 1.0.
sanity_1 = jointCDFTheo(rMax, 2*pi, Lx, Ly, u);
disp(['Sanity check 1: the CDF of the bivariate PDF evaluated at r=rMax & \theta = 2*\pi is:', ...
    num2str(sanity_1)])

sanity_2 = intJointPDF_R_Theta(epsilon, hi, qfun, Chil, Chig, Lx, Ly, u, rMax, 'option1');
disp(['Sanity check 2: the integral of the bivariate PDF is (method 1):', ...
    num2str(sanity_2)])

sanity_3 = intJointPDF_R_Theta(epsilon, hi, qfun, Chil, Chig, Lx, Ly, u, rMax, 'option2');
disp(['Sanity check 3: the integral of the bivariate PDF is (method 2):', ...
    num2str(sanity_3)])



%% Plot statistics
if debug.marginalPDF_R
    hf_PDF_R = figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
    plot(pdf_R.r, pdf_R.pdf, 'bo'); hold on; 
    title(sprintf('Marginal PDF R'));
    xlabel('$r (m)$', cell_latex{:}); ylabel('$f_{R}(r)$', cell_latex{:});grid on;
    set(gca,'fontsize',fontSize)
end

if debug.conditionedJointPDF
    hf_condPDF = figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
    linspec_sim = {'bo', 'ro', 'go', 'co'};
    legend_str_conditionedJointPDF = cell(nQuad*2, 1);
    for iQ = 1:nQuad
        plot(pdfR_theta{iQ}.x, pdfR_theta{iQ}.pdf, linspec_sim{iQ});
        if iQ == 1, hold on; end
        theta_min_deg = min(limit_angles_per_quadrant_deg(iQ, :));
        theta_max_deg = max(limit_angles_per_quadrant_deg(iQ, :));
        legend_str_conditionedJointPDF{iQ, 1} = sprintf('Sim I = [%d, %d]º', theta_min_deg, ...
            theta_max_deg);
    end
    title(sprintf('Conditional PDF'));
    xlabel('$r (m)$', cell_latex{:}); ylabel('$f_{R}(R| \Theta \in I)$', cell_latex{:});grid on;

    set(gca,'fontsize',fontSize)
end

if debug.jointPDF
    % Plot joint PDF distance, angle
    hf2a = figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
    surf(stR.x, stR.y/pi, stR.pdf);
    xlabel('$r (m)$', cell_latex{:}); ylabel('$\theta/pi$', cell_latex{:}); 
    zlabel('$f_{R,\Theta}(r,\theta)$', cell_latex{:})
    grid; title('Estimated joint PDF of distance and angle $(R, \Theta)$', cell_latex{:})
    set(gca,'fontsize',fontSize)
end

if debug.jointCDF
    % Plot joint CDF distance, angle
    hf2b = figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
    surf(stR.x, stR.y/pi, stR.cdf);
    xlabel('$r (m)$', cell_latex{:}); ylabel('$\theta/\pi$', cell_latex{:}); 
    zlabel('$F_{R,\Theta}(r,\theta)$', cell_latex{:})
    grid; title('Estimated joint CDF of distance and angle $(R, \Theta)$', cell_latex{:})
    set(gca,'fontsize',fontSize)
end

if debug.marginalPDFTheta
    % Plot marginal PDF of the angle on radians
    hf3 = figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
    plot(xi_Thetha/pi, f_Thetha, 'o'); 
    xlabel('$\theta/\pi$', cell_latex{:}); ylabel('$f_{\Theta}(\theta)$', cell_latex{:}); 
    hold on; grid on; title('Estimated marginal PDF of the angle (radians)', cell_latex{:})
    set(gca,'fontsize',fontSize)
end

if debug.marginalCDFTheta
    % Plot marginal CDF of the angle 
    hf6 = figure('DefaultTextFontName', font, 'DefaultAxesFontName', font); 
    plot(ThetaVct/pi, cdfTheta, 'o'); 
    xlabel('$\theta/\pi$', cell_latex{:}); ylabel('$F_{\Theta}(\theta)$', cell_latex{:}); 
    hold on; grid on; title('Marginal CDF of the angle', cell_latex{:})
    set(gca,'fontsize',fontSize)
end

if debug.conditionedJointCDF
    % Plot joint CDF R and Theta for each quadrant
    hf_Q = cell(nQuad);
    linspec_sim = {'o', 'ro', 'go'};
    legend_str = cell(nAng, 2);
    for iQ = 1:nQuad
        hf_Q{iQ} = figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
        for iAng = 1:nAng
            plot(RVct, cdfR_theta{iQ, iAng}, linspec_sim{iAng}); 
            if iAng == 1, hold on; end
            legend_str{iAng, 1} = ['sim. \theta_0 = ', ...
                num2str(angles_per_quadrant_deg(iQ, iAng)), 'º'];
        end
        
        title(sprintf('Joint CDF for \\theta (degrees) = (%s). Quadrant: %d', ...
            num2str(angles_per_quadrant_deg(iQ, :)), iQ))
        xlabel('$r (m)$', cell_latex{:}); ylabel('$F_{R, \Theta}(R, \theta = \theta_0)$', cell_latex{:});grid on;
    end 
    set(gca,'fontsize',fontSize)
end

%% Theoretical results
ThetaTheoVct = linspace(0+1e-3, thetaMax, 100);
RTheoVct = linspace(0+1e-3, rMax, 100);
[RR,TTheta] = meshgrid(RTheoVct, ThetaTheoVct);

if debug.marginalPDF_R
    figure(hf_PDF_R);
    pdf_R_Theo = joint_pdf_R_cond_O(epsilon, hi, Chil, Chig, RTheoVct, 0, 2*pi, ...
        Lx, Ly, u, version_R_cond_O);
    plot(RTheoVct, pdf_R_Theo, 'b');
    legend('Simulation', 'Analysis');
    title('Marginal PDF of the distance')
end

if debug.conditionedJointPDF
    figure(hf_condPDF);
    linspec_sim = {'b', 'r', 'g', 'c'};
    for iQ = 1:nQuad
        % Compute theoretical PDF
        theta_min = min(limit_angles_per_quadrant_deg(iQ, :))/180*pi;
        theta_max = max(limit_angles_per_quadrant_deg(iQ, :))/180*pi;
        pdfR_cond_O = joint_pdf_R_cond_O(epsilon, hi, Chil, Chig, RTheoVct, ...
            theta_min, theta_max, Lx, Ly, u, version_R_cond_O);
        plot(RTheoVct, pdfR_cond_O, linspec_sim{iQ});
        if iQ == 1, hold on; end
        theta_min_deg = min(limit_angles_per_quadrant_deg(iQ, :));
        theta_max_deg = max(limit_angles_per_quadrant_deg(iQ, :));
        legend_str_conditionedJointPDF{iQ+nQuad, 1} = sprintf('Theor I = [%d, %d]º', theta_min_deg, ...
            theta_max_deg);
    end

    title(sprintf('Conditional PDF'));
    legend(legend_str_conditionedJointPDF(:));
end

% https://www.math24.net/basic-trigonometric-inequalities
% Look in this link how solving trigonometric inequalities
if debug.jointCDF
    jointCDFTheoVct = jointCDFTheo(RR, TTheta, Lx, Ly, u);
    hf_Teo1 = figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
    surf(RR, TTheta/pi, jointCDFTheoVct); 
    title('Theoretical joint CDF $(R, \Theta)$', cell_latex{:})
    xlabel('$r (m)$', cell_latex{:}); ylabel('$\theta/\pi$', cell_latex{:}); zlabel('$F_{R,\Theta}(r,\theta)$', cell_latex{:})
    set(gca,'fontsize',fontSize)
end

if debug.conditionedJointCDF
    cdfR_theta_Theo = cell(nQuad, nAng);
    linspec_theo = {'b', 'r', 'g'};
    for iQ = 1:nQuad
        figure(hf_Q{iQ}); 
        for iAng = 1:nAng
            cdfR_theta_Theo{iQ, iAng} = ...
                jointCDFTheo(RTheoVct, ones(size(ThetaTheoVct))...
                *angles_per_quadrant_deg(iQ, iAng)/180*pi, Lx, Ly, u);
            plot(RTheoVct, cdfR_theta_Theo{iQ, iAng}, linspec_theo{iAng}); 
            if iAng == 1, hold on; end
            legend_str{iAng, 2} = ['theor. \theta_0 = ', ...
                num2str(angles_per_quadrant_deg(iQ, iAng)), 'º'];
        end
        legend(legend_str(:))
    end 
    set(gca,'fontsize',fontSize)
end

if debug.marginalCDFTheta
    % Plot marginal CDF of theta
    vs = 2;
    marginalCDFThetaVct = cdfThetaTheo(epsilon, hi, ThetaTheoVct, Lx, Ly, u, vs);
    figure(hf6); plot(ThetaTheoVct/pi, marginalCDFThetaVct, 'b');
    set(gca,'fontsize',fontSize)
end 

if debug.jointPDF
    jointPDFTheoVct = jointPDFTheo(Chil, Chig, RR, TTheta, Lx, Ly, u);
    hf_Teo2 = figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
    surf(RR, TTheta/pi, jointPDFTheoVct); 
    xlabel('$r (m)$', cell_latex{:}); ylabel('$\theta/\pi$', cell_latex{:}); zlabel('$f_{R,\Theta}(r,\theta)$', cell_latex{:})
    grid on; title('Theoretical joint PDF of distance and angle $(R, \Theta)$', cell_latex{:})
    set(gca,'fontsize',fontSize)
end

if debug.marginalPDFTheta
    % Marginal PDF of angle
    marginalPDFThetaVct = pdfThetaTheo(ThetaTheoVct, Lx, Ly, u);
    figure(hf3); plot(ThetaTheoVct/pi, marginalPDFThetaVct, 'b');
    set(gca,'fontsize',fontSize)
end