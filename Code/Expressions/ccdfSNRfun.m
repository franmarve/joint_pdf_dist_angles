% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin-Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function computes the analytical ccdf of the SNR for a given (scalar)
% threshold value, t, according to Corollary 9 and Approximation 2 of [1]
% *************************************************************************
% REFERENCE:
% [1] Francisco J. Martin-Vega, Gerardo Gomez, David Morales-Jimenez, 
% F. Javier Lopez-Martinez and Mari Carmen Aguayo-Torres, "Joint 
% Distribution of Distance and Angles in Finite Wireless Networks", acepted
% for publication in IEEE Transactions on Vehicular Technology, 2023.

function pSNR = ccdfSNRfun(t, hi, Chil, Chig, p, sP, rMax, method)

sl = cell(8, 1);
sg = cell(8, 1);

if strcmp(method, '2dint')
    % This block of code represents the straightforward approach of 
    % computing the CCDF of the SNR as a double integral over the domain of
    % the joint PDF of distances and azimuth angle. This method leads to 
    % inaccurate results since the numerical integration method does not 
    % perform well in the discontinuities. The option 'avoidDiscont' should
    % be used instead.  
    fun = @(r, theta) ...
        jointPDFTheo(Chil, Chig, r, theta, sP.Lx, sP.Ly, sP.u, 2) ...
            .*exp( -(t.*sP.N0.*(sP.tau.*r).^sP.alpha) ./ ...
            (sP.rhot.*gServ(theta, p, sP)) );
    pSNR = integral2(fun, 0, rMax, 0, 2*pi);

elseif strcmp(method, 'avoidDiscont')
    % This block of code computes the CCDF of the SNR as per [1], 
    % Corollary 9
    pSNR = 0;
    for ii = 1:8
        fun = @(r, theta) ...
            exp( -(t.*sP.N0.*(sP.tau.*r).^sP.alpha)...
            ./(sP.rhot.*gServ(theta, p, sP)));
        slScalar = @(rr) integral(@(theta) fun(rr, theta), ...
            Chil{ii,1}(rr), Chil{ii,2}(rr));
        sl{ii} = @(r) arrayfun(slScalar, r);
        sg{ii} = @(r) cF(fun, r, Chig{ii,1}, Chig{ii,2});
        pSNR = pSNR ...
            + integral(@(r) sl{ii}(r).*r, 0, hi(ii)) ...
            + integral(@(r) sg{ii}(r).*r, hi(ii), rMax);
    end
    pSNR = pSNR/(sP.Lx*sP.Ly);

elseif strcmp(method, 'degeneratedLy')
    fun1 = @(r, theta) exp( -(t.*sP.N0.*(sP.tau.*r).^sP.alpha) ./ ...
            (sP.rhot.*gServ(pi/2, p, sP)) )/sP.Ly;
    fun2 = @(r, theta) exp( -(t.*sP.N0.*(sP.tau.*r).^sP.alpha) ./ ...
            (sP.rhot.*gServ(3*pi/2, p, sP)) )/sP.Ly;
    pSNR = integral(fun1, 0, sP.Ly/2) + integral(fun2, 0, sP.Ly/2);

elseif strcmp(method, 'approximationLy')
    lambda = @(phi) t.*sP.N0./(sP.rhot.*gServ(phi, p, sP));
    pSNR = (gamma(1/(sP.alpha)) ...
        -gamma(1/(sP.alpha))*...
        gammainc(2^(-sP.alpha)*lambda(pi/2)*(sP.Ly*sP.tau)^(sP.alpha), ...
        1/(sP.alpha), 'upper'))/...
        (sP.Ly*sP.alpha*sP.tau*(lambda(pi/2)^(1/(sP.alpha)))) ...
        ...
        + (gamma(1/(sP.alpha)) ...
        -gamma(1/(sP.alpha))*...
        gammainc(2^(-sP.alpha)*lambda(3*pi/2)*(sP.Ly*sP.tau)^(sP.alpha), ...
        1/(sP.alpha), 'upper'))/...
        (sP.Ly*sP.alpha*sP.tau*(lambda(3*pi/2)^(1/(sP.alpha))));


else
    error('ccdfSNR: Invalid method');
end