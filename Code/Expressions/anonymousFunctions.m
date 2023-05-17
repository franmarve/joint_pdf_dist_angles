% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega 
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function creates some intermediate functions defined in [1] 
% that are needed to compute the distribuitions of distance and angles. 
% *************************************************************************
% REFERENCES:
% [1] Francisco J. Martin-Vega, Gerardo Gomez, David Morales-Jimenez, 
% F. Javier Lopez-Martinez and Mari Carmen Aguayo-Torres, "Joint 
% Distribution of Distance and Angles in Finite Wireless Networks", acepted
% for publication in IEEE Transactions on Vehicular Technology, 2023.

function [Chil, Chig, gfun, qfun, hi, epsilon] = anonymousFunctions(Lx, Ly, u)

hxp = Lx/2 - u(1);  % Latex notation: h_x^{+}
hxm = -Lx/2 - u(1); % Latex notation: h_x^{-}
if hxm == 0, hxm = -0; end
hyp = Ly/2 - u(2);  % Latex notation: h_y^{+}
hym = -Ly/2 - u(2); % Latex notation: h_y^{-}
if hym == 0, hym = -0; end

hi = [hxp, hxp, -hxm, -hxm, hyp, hyp, -hym, -hym];

% Eq. (22) of [1]
Chil{1,1} = @(r) zeros(size(r));
Chil{1,2} = @(r) atan(hyp/hxp)*ones(size(r));
Chil{2,1} = @(r) (atan(hym/hxp) + 2*pi)*ones(size(r));
Chil{2,2} = @(r) 2*pi*ones(size(r));
Chil{3,1} = @(r) pi*ones(size(r));
Chil{3,2} = @(r) (atan(hym/hxm) + pi)*ones(size(r));
Chil{4,1} = @(r) (atan(hyp/hxm) + pi)*ones(size(r));
Chil{4,2} = @(r) pi*ones(size(r));
Chil{5,1} = @(r) atan(hyp/hxp)*ones(size(r));
Chil{5,2} = @(r) pi/2*ones(size(r));
Chil{6,1} = @(r) pi/2*ones(size(r));
Chil{6,2} = @(r) (atan(hyp/hxm) + pi)*ones(size(r));
Chil{7,1} = @(r) (atan(hym/hxm) + pi)*ones(size(r));
Chil{7,2} = @(r) (3*pi)/2*ones(size(r));
Chil{8,1} = @(r) (3*pi)/2*ones(size(r));
Chil{8,2} = @(r) (atan(hym/hxp) + 2*pi)*ones(size(r));

Chig{1,1} = @(r) acos(hxp./r);
Chig{1,2} = @(r) atan(hyp/hxp)*ones(size(r));
Chig{2,1} = @(r) (atan(hym/hxp) + 2*pi)*ones(size(r));
Chig{2,2} = @(r) 2*pi - acos(hxp./r);
Chig{3,1} = @(r) 2*pi - acos(hxm./r);
Chig{3,2} = @(r) (atan(hym/hxm) + pi)*ones(size(r));
Chig{4,1} = @(r) (atan(hyp/hxm) + pi)*ones(size(r));
Chig{4,2} = @(r) acos(hxm./r);
Chig{5,1} = @(r) atan(hyp/hxp)*ones(size(r));
Chig{5,2} = @(r) asin(hyp./r);
Chig{6,1} = @(r) (pi - asin(hyp./r));
Chig{6,2} = @(r) (atan(hyp/hxm) + pi)*ones(size(r));
Chig{7,1} = @(r) (atan(hym/hxm) + pi)*ones(size(r));
Chig{7,2} = @(r) pi - asin(hym./r);
Chig{8,1} = @(r) 2*pi + asin(hym./r);
Chig{8,2} = @(r) (atan(hym/hxp) + 2*pi)*ones(size(r));

% Eq. defined with Corollary 4
epsilon{1,1} = @(theta) zeros(size(theta));
epsilon{1,2} = @(theta) min(theta, atan(hyp/hxp));
epsilon{2,1} = @(theta) (atan(hym/hxp) + 2*pi)*ones(size(theta));
epsilon{2,2} = @(theta) theta;
epsilon{3,1} = @(theta) pi*ones(size(theta));
epsilon{3,2} = @(theta) min(theta, (atan(hym/hxm) + pi)*ones(size(theta)));
epsilon{4,1} = @(theta) (atan(hyp/hxm) + pi)*ones(size(theta));
epsilon{4,2} = @(theta) min(theta, pi*ones(size(theta)));
epsilon{5,1} = @(theta) atan(hyp/hxp)*ones(size(theta));
epsilon{5,2} = @(theta) min(theta, pi/2*ones(size(theta)));
epsilon{6,1} = @(theta) pi/2*ones(size(theta));
epsilon{6,2} = @(theta) min(theta, (atan(hyp/hxm) + pi)*ones(size(theta)));
epsilon{7,1} = @(theta) (atan(hym/hxm) + pi)*ones(size(theta));
epsilon{7,2} = @(theta) min(theta, (3*pi)/2*ones(size(theta)));
epsilon{8,1} = @(theta) (3*pi)/2*ones(size(theta));
epsilon{8,2} = @(theta) min(theta, (atan(hym/hxp) + 2*pi)*ones(size(theta)));

gfun = cell(8, 2);
qfun = cell(8, 2);
theta = @(x) x;
for i = 1:8
    % Index j = 1 stands for superscrip '<' (less) whereas
    % j = 2 is for '>=' (greater or equal)
    for j = 1:2 
        if j == 1
            gfun{i, j} = @(r) ...
                F(theta, max(0, Chil{i,1}(r)), min(2*pi, Chil{i,2}(r)));
        else
            gfun{i, j} = @(r) ...
                F(theta, max(0, Chig{i,1}(r)), min(2*pi, Chig{i,2}(r)));
        end
        qfun{i, j} = @(r) r.*gfun{i, j}(r);
    end
end


