% *************************************************************************
% Francisco J. Martin-Vega, frmvega@gmail.com
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% 
% *************************************************************************

function pdfR_cond_O = joint_pdf_R_cond_O(epsilon, hi, Chil, Chig, RVtc, ...
    theta_min, theta_max, Lx, Ly, u, version)

if nargin < 9
    version = 1;
end

nElm = length(RVtc);
pdfR_cond_O = zeros(nElm, 1);

PrO = cdfThetaTheo(epsilon, hi, theta_max, Lx, Ly, u) ...
    - cdfThetaTheo(epsilon, hi, theta_min, Lx, Ly, u);

if version == 1
    for i = 1:nElm
        fun = @(theta) jointPDFTheo(Chil, Chig, RVtc(i), theta, Lx, Ly, u);
        pdfR_cond_O(i) = integral(fun, theta_min, theta_max)/PrO;
    end
else
    % version 2
    hxp = Lx/2 - u(1);  % Latex notation: h_x^{+}
    hxm = -Lx/2 - u(1); % Latex notation: h_x^{-}
    if hxm == 0, hxm = -0; end
    hyp = Ly/2 - u(2);  % Latex notation: h_y^{+}
    hym = -Ly/2 - u(2); % Latex notation: h_y^{-}
    if hym == 0, hym = -0; end
    
    hi = [hxp, hxp, -hxm, -hxm, hyp, hyp, -hym, -hym];

    % Anonymus function
    tgfun = cell(8, 2);
    theta = @(x) x;
    for i = 1:8
        % Index j = 1 stands for superscrip '<' (less) whereas
        % j = 2 is for '>=' (greater or equal)
        for j = 1:2 
            if j == 1
                tgfun {i, j} = @(r) ...
                    F(theta, max(theta_min, Chil{i,1}(r)), min(theta_max, Chil{i,2}(r)));
            else
                tgfun {i, j} = @(r) ...
                    F(theta, max(theta_min, Chig{i,1}(r)), min(theta_max, Chig{i,2}(r)));
            end
        end
    end
    
    pdfR_cond_O = double(RVtc < hi(1)).*tgfun {1,1}(RVtc) ...
        + double(RVtc >= hi(1)).*tgfun {1,2}(RVtc);
    
    for i = 2:8
	    pdfR_cond_O = pdfR_cond_O + double(RVtc < hi(i)).*tgfun {i,1}(RVtc) ...
        + double(RVtc >= hi(i)).*tgfun {i,2}(RVtc);
    end
    pdfR_cond_O = pdfR_cond_O.*RVtc/(Lx*Ly);
    pdfR_cond_O = pdfR_cond_O/PrO;
    % pdfR_cond_O = pdfR_cond_O.*RVtc/(Lx*Ly)/PrO;

end