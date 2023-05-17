% *************************************************************************
% Francisco J. Martin-Vega, fjmv@uma.es
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function integrates the joint PDF of distance and angle. 
% *************************************************************************

function result = intJointPDF_R_Theta(epsilon, hi, qfun, Chil, Chig, Lx, Ly, u, rMax, version)

hxp = Lx/2 - u(1);  % Latex notation: h_x^{+}
hxm = -Lx/2 - u(1); % Latex notation: h_x^{-}
if hxm == 0, hxm = -0; end
hyp = Ly/2 - u(2);  % Latex notation: h_y^{+}
hym = -Ly/2 - u(2); % Latex notation: h_y^{-}
if hym == 0, hym = -0; end

hi = [hxp, hxp, -hxm, -hxm, hyp, hyp, -hym, -hym];

if strcmp(version, 'option1')
    result = 0;
    for i = 1:8
        I1 = integral(qfun{i,1}, 0, hi(i));
        I2 = integral(qfun{i,2}, hi(i), rMax);
        result = result + I1 + I2;
    end
    
    result = result/(Lx*Ly);
elseif strcmp(version, 'option2')
    pdfR = @(r) joint_pdf_R_cond_O(epsilon, hi, Chil, Chig, r, 0, 2*pi, Lx, Ly, u, 2);
    result = integral(pdfR, 0, rMax);
else
    error('intJointPDF_R_Theta: Error unkonwn version')
end
