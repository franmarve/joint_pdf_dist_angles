% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin-Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function places and gets the distances of a BPP towards a reference
% 3D location u, in cartesian coordinates. It is considered the heigh of 
% the reference node, u, and of the rest of randon nodes of the BPP. 
% *************************************************************************
function [D, Theta, Psi, R, V, hf] = ...
    randomCoord(plotFig, Lx, Ly, N, u, vz, psif)

font = 'Times new roman';

% Sanity checks: input arguments 
if nargin == 5 % 2D case
    % Empty variables in 2D case
    D = []; 
    Psi = [];
    if length(u) ~= 2
        error(['randomCoord: The location, u, must have 2 cartesian '...
            'coordinates'])
    end
elseif nargin == 7 % 3D case
    if length(u) ~= 3
        error(['randomCoord: The location, u, must have 3 cartesian'...
            ' coordinates'])
    end
else
    error(['randomCoord: wrong number of input arguments'])
end

hxp = Lx/2 - u(1);  % Latex notation: h_x^{+}
hxm = -Lx/2 - u(1); % Latex notation: h_x^{-}
hyp = Ly/2 - u(2);  % Latex notation: h_y^{+}
hym = -Ly/2 - u(2); % Latex notation: h_y^{-}

% Construct the vertex of the rectangle
V(:, 1) = [hxm, hxp,  hxp, hxm];
V(:, 2) = [hym, hym,  hyp, hyp];

% Place in polygon with vertexs V
Phi = placeInPolygon(N, V);

% Number of vertices
[nV, ~] = size(V);

% Plot points
if plotFig
    if nargin == 5 % 2D case
        hf = figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
        plot(Phi(:, 1), Phi(:, 2), 'x'); hold on; grid; 
        xlabel('x (m)'); ylabel('y (m)'); % axis square;
        plot([V(:, 1); V(1, 1)], [V(:, 2); V(1, 2)], 'r')
        axis([min(V(:, 1)) - 1, max(V(:, 1)) + 1, min(V(:, 2)) ...
            - 1, max(V(:, 2)) + 1])
        plot([hxm hxp], [0, 0], 'k'); plot([0, 0], [hym hyp], 'k');
        axis equal
    else % 3D case
        hf = figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
        plot3(Phi(:, 1), Phi(:, 2), vz*ones(N, 1), 'x'); hold on; grid; 
        xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)')
        plot3([V(:, 1); V(1, 1)], [V(:, 2); V(1, 2)], zeros(nV + 1, 1),'r')
        axis([min(V(:,1))-1, max(V(:,1))+1, min(V(:,2))-1, max(V(:,2))+1])
        plot3([hxm hxp], [0, 0], [0, 0], 'k'); 
        plot3([0, 0], [hym hyp], [0, 0], 'k');
        plot3([0, 0], [0, 0], [0, u(3)], '-*k','LineWidth',5, ...
            'MarkerSize',10)
        axis equal
    end
else
    hf = [];
end 

% Convert to polar coordinates
R = sqrt(Phi(:, 1).^2 + Phi(:, 2).^2);
Theta = cart2pol(Phi(:, 1), Phi(:, 2));
Theta(Theta < 0) = Theta(Theta < 0) + 2*pi;

if nargin == 7 % 3D case
    % Conversion into spherical coordinates
    D = sqrt(R.^2 + (u(3) - vz).^2);
    Psi = psif(R, u(3)*ones(size(R)), vz*ones(size(R)));
end
