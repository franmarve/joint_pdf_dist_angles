% *************************************************************************
% Francisco J. Martin-Vega, fjmv@uma.es
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% It obtains the joint PDF of distance and angle between a reference point
% and a set of randomly placed points with uniform distribution. 
% *************************************************************************

function outVct = jointPDFTheo(Chil, Chig, r, theta, Lx, Ly, u)


hxp = Lx/2 - u(1);  % Latex notation: h_x^{+}
hxm = -Lx/2 - u(1); % Latex notation: h_x^{-}
if hxm == 0, hxm = -0; end
hyp = Ly/2 - u(2);  % Latex notation: h_y^{+}
hym = -Ly/2 - u(2); % Latex notation: h_y^{-}
if hym == 0, hym = -0; end

hi = [hxp, hxp, -hxm, -hxm, hyp, hyp, -hym, -hym];

outVct = double(r < hi(1)).*...
    double(Chil{1,1}(r) <= theta & theta < Chil{1,2}(r)) ...
    + double(r >= hi(1)).*...
    double(Chig{1,1}(r) <= theta & theta < Chig{1,2}(r));

for i = 2:8
    outVct = outVct + double(r < hi(i)).*...
        double(Chil{i,1}(r) <= theta & theta < Chil{i,2}(r)) ...
        + double(r >= hi(i)).*...
        double(Chig{i,1}(r) <= theta & theta < Chig{i,2}(r));
end

outVct = outVct.*r/(Lx*Ly);
