% *************************************************************************
% Francisco J. Martin-Vega, fjmv@uma.es
% Lab 1.3.5., Dpto. of Ingenieria de Comunicaciones. University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function places N points within a given polygon defined by its 
% vertex. 
% INPUTS:
% +) N: Number of points to be placed
% +) V: Vextex that define the polygon. Each vertex [xv, yv] is a row
% vector of the matrix V [nV x 2], which is defined by nV vertexs. 
% *************************************************************************

function Phi = placeInPolygon(N, V)

%% Compute the rectangle that contains this polygon
xMax = max(V(:, 1));
xMin = min(V(:, 1));
yMax = max(V(:, 2));
yMin = min(V(:, 2));


%% Place 2*N points within that rectangle
Phi = rand(2*N, 2).*[xMax - xMin, yMax - yMin] + [xMin, yMin];

% Check the points that have fallen within the polygon
in = inpolygon(Phi(:,1) ,Phi(:,2), V(:,1), V(:,2));

% Discard points outside the sector
Phi = Phi(in, :);

% Get the number of points already within the sector
points_inside = sum(in);

%% Generate more points until there are more than N points within the sector

% Scale parameter on the number of points generated per iteration
scale = 5;
while points_inside < N
    % Add floor(N/scale) random points within the rectangle
    Phi_added = rand(ceil(N/scale), 2).*[xMax - xMin, yMax - yMin] + [xMin, xMin]; 
    
    % Check how many of the added points are within the polygon
    in = inpolygon(Phi_added(:,1) ,Phi_added(:,2), V(:,1), V(:,2));
    
    % Add to the point set only the points that have fallen within the 
    % polygon in the current iteration
    Phi = [Phi; Phi_added(in, :)];
    
    % Update the number of points inside the sector
    points_inside = points_inside + sum(in);
end

% Select only the first N points within the sector
if points_inside > N
    Phi = Phi(1:N, :);
end



