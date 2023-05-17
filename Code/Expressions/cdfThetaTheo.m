% *************************************************************************
% Francisco J. Martin-Vega, fjmv@uma.es
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% It obtains the marginal CDF angle between a reference point
% and a set of randomly placed points with uniform distribution. 
% *************************************************************************
function cdfThetaVct = cdfThetaTheo(epsilon, hi, theta, Lx, Ly, u, version)

if nargin < 7
    version = 2;
end

% Function handle for inverse tan(x), 1/tan(x)
itan = @(x) 1./tan(x);

if version == 1
    hxp = Lx/2 - u(1);  % Latex notation: h_x^{+}
    hxm = -Lx/2 - u(1); % Latex notation: h_x^{-}
    if hxm == 0, hxm = -0; end
    hyp = Ly/2 - u(2);  % Latex notation: h_y^{+}
    hym = -Ly/2 - u(2); % Latex notation: h_y^{-}
    if hym == 0, hym = -0; end    
    
    I_rToInf_xp = hxp^2/2.*( tan(min(theta, atan(hyp/hxp))) ...
        + F(@tan, atan(hym/hxp) +2*pi, theta)  );
    
    I_rToInf_xm = hxm^2/2.*( F(@tan, pi, min(theta, atan(hym/hxm) + pi)) ...
        + F(@tan, atan(hyp/hxm) + pi, min(theta, pi)) );
    
    I_rToInf_yp = -hyp^2/2.*( ...
        F(itan, atan(hyp/hxp), min(theta, pi/2)) ...
        + F(itan, pi/2, min(theta, atan(hyp/hxm) +pi)) );
    
    I_rToInf_ym = -hym^2/2.*( F(itan, atan(hym/hxm) +pi, min(theta, 3*pi/2))...
        + F(itan,3*pi/2,min(theta, atan(hym/hxp) +2*pi)) );
    
    cdfThetaVct = 1/(Lx*Ly)*(I_rToInf_xp + I_rToInf_xm + I_rToInf_yp + I_rToInf_ym);
else
    % Memory allocation
    cdfThetaVct = zeros(size(theta));
    
    % Margianl CDF of the azimuth angle theta as per Corollary 4
    for i = 1:8
	    if (i < 5) % i = 1:4
            cdfThetaVct = cdfThetaVct + hi(i).^2 ...
            .*F(@tan, epsilon{i,1}(theta), epsilon{i,2}(theta));
        else % i = 5:8
            cdfThetaVct = cdfThetaVct - hi(i).^2 ...
            .*F(itan, epsilon{i,1}(theta), epsilon{i,2}(theta));
        end
    end
    cdfThetaVct = cdfThetaVct/(2*Lx*Ly);
end
end