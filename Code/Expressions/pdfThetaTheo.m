% *************************************************************************
% Francisco J. Martin-Vega, fjmv@uma.es
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% It obtains the marginal CDF angle between a reference point
% and a set of randomly placed points with uniform distribution. 
% *************************************************************************
function outVct = pdfThetaTheo(theta, Lx, Ly, u)

hxp = Lx/2 - u(1);  % Latex notation: h_x^{+}
hxm = -Lx/2 - u(1); % Latex notation: h_x^{-}
if hxm == 0, hxm = -0; end
hyp = Ly/2 - u(2);  % Latex notation: h_y^{+}
hym = -Ly/2 - u(2); % Latex notation: h_y^{-}
if hym == 0, hym = -0; end

Ipp_rToInf_xp = hxp^2./(2*(cos(theta)).^2).*( ...
    double(0 <= theta & theta < atan(hyp/hxp)) ...
    + double(atan(hym/hxp) + 2*pi < theta & theta <= 2*pi) );

Ipp_rToInf_xm = hxm^2./(2*(cos(theta)).^2).*( ...
    double(pi < theta & theta < atan(hym/hxm) +pi) ...
    + double(atan(hyp/hxm) +pi < theta & theta < pi) );

Ipp_rToInf_yp = hyp^2./(2*(sin(theta)).^2).*( ...
    double (atan(hyp/hxp) < theta & theta < pi/2) ...
    + double(pi/2 < theta & theta < atan(hyp/hxm) + pi) );

Ipp_rToInf_ym = hym^2./(2*(sin(theta)).^2).*( ...
    double (atan(hym/hxm) +pi < theta & theta <= 3*pi/2) ...
    + double(3*pi/2 < theta & theta < atan(hym/hxp) +2*pi) );

outVct = 1/(Lx*Ly)*(Ipp_rToInf_xp + Ipp_rToInf_xm ...
    + Ipp_rToInf_yp + Ipp_rToInf_ym);