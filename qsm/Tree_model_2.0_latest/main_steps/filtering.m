function Pass = filtering(P0,r1,n1,r2,n2,Scaling,Comp)

% ---------------------------------------------------------------------
% FILTERING.M       Filters noise from point clouds.
%
% Version 1.3
% Author        Pasi Raumonen
% Created       21 February 2014
%
% This software cannot be distributed without the permission of P. Raumonen
% and it is only for non-commercial use. These restrictions holds also for
% the modules or subprograms of the software:
% PARTITION_OF_POINT_CLOUD.M, COVER_SETS_AND_NEIGHBORS.M,
% CONNECTED_COMPONENTS.M
% ---------------------------------------------------------------------

% Performs an initial filtering of the point cloud.
% 
% At first the possible multiples of points or NaNs are removed.
% 
% Secondly those points which belong to r1-balls that have at least
% n1 points them are included. Comprehensive filtering here means that
% the r1-ball is defined for every point whereas non-comprehensive means
% that only a cover of point cloud with r1-balls is defined
%
% Lastly small components (less than n2 r2-balls) are removed.
%
% For both filtering procedures the given minimum number (n1 and n2) 
% applies for the first three meters above the ground and after that it is
% modified for each meter by the change in the mean point density. That is,
% if the point density decreases, which usually is the case, then the
% minimum number is also decreased accordingly.
%
% Inputs:
% P0        Point cloud
% r1        Radius of the balls used in the first filtering
% n1        Minimum number of points in the accepted balls of the first filtering
% r2        Radius of the balls used in the second filtering
% n2        Minimum number of balls in the components passing the second filtering
% Scaling   If true, the first filtering threshold "n1" is scaled along the 
%               height with average point density
% Comp      If true, does the first filtering process for every point
%
% Outputs:
% Pass      Logical vector indicating points passing the filtering


% Default firts filtering is not comprehensive and with scaling
if nargin == 5
    Scaling = true;
    Comp = false; 
elseif nargin == 6
    Comp = false;
end


%% Remove the possible multiples of points and NaNs
P = unique(P0,'rows');
I = all(~isnan(P),2);
P = P(I,:);
np = size(P,1);
np0 = np;


%% Partition the point cloud into cubes
[partition,CC] = partition_of_point_cloud(P,r1);

if ~Comp
    % Do the first filtering with a cover, not comprehensive filtering
    %% Generate a cover and determine the largest (number of points) ball for each point
    NotInspected = true(np,1);
    NumOfPoints = zeros(np,1);
    r1 = r1^2;
    for i = 1:np
        if NotInspected(i)
            points = partition(CC(i,1)-1:CC(i,1)+1,CC(i,2)-1:CC(i,2)+1,CC(i,3)-1:CC(i,3)+1);
            points = vertcat(points{:,:});
            cube = P(points,:);
            dist = (P(i,1)-cube(:,1)).^2+(P(i,2)-cube(:,2)).^2+(P(i,3)-cube(:,3)).^2;
            J = dist < r1;
            I = points(J);
            NotInspected(I) = false;
            N = NumOfPoints(I);
            n = length(I);
            K = N < n;
            NumOfPoints(I(K)) = n;
        end
    end
else
    % Do the first filtering comprehensively by defining the neighborhoods
    % for all points
    
    %% Generate the balls and determine the number of points for each point
    NumOfPoints = zeros(np,1);
    r1 = r1^2;
    for i = 1:np
        points = partition(CC(i,1)-1:CC(i,1)+1,CC(i,2)-1:CC(i,2)+1,CC(i,3)-1:CC(i,3)+1);
        points = vertcat(points{:,:});
        cube = P(points,:);
        dist = (P(i,1)-cube(:,1)).^2+(P(i,2)-cube(:,2)).^2+(P(i,3)-cube(:,3)).^2;
        J = dist < r1;
        NumOfPoints(i) = nnz(J);
    end
end
clearvars partition CC


%% First filtering
if Scaling
    % Use smaller treshold for upper parts of the tree
    % Change the treshold "n1" according the average points in the cover
    % sets for every meter
    hmin = min(P(:,3));
    hmax = max(P(:,3));
    H = ceil(hmax-hmin);
    D = zeros(H,1);
    J = false(np,1);
    Pass = false(np,1);
    A = (1:1:np)';
    for i = 1:H
        I = P(:,3) < hmin+i;
        K = I&~J;
        J = I;
        D(i) = ceil(mean(NumOfPoints(K)));
        if i <= 2
            M = A(K);
            N = NumOfPoints(K) >= n1;
            M = M(N);
            Pass(M) = true;
        else
            M = A(K);
            m = max(ceil(n1*D(i)/D(2)),ceil(n1/3));
            N = NumOfPoints(K) >= m;
            M = M(N);
            Pass(M) = true;
        end
    end
else
    Pass = NumOfPoints >= n1;
end
clearvars NumOfPoints

% Display filtering results
np = nnz(Pass);
nf = np0-np;
str = ['    All points: ',num2str(np0),', First filtering: ',num2str(nf),', Points left: ',num2str(np)];
disp(str)


%% Cover the point cloud with r2-balls for the second filtering
[Balls,Neighbors] = cover_sets(P(Pass,:),r2,1.5*r2,0);


%% Determine the separate components
Components = connected_components(Neighbors,0,n2);
clearvars Neighbors


%% Second filtering
B = vertcat(Components{:});
I = vertcat(Balls{B});
clearvars Balls Neighbors B components
J = false(np,1);
J(I) = true;
I = (1:1:np0)';
I = I(Pass);
I = I(~J);
Pass(I) = false;

% Display filtering results
nf = np-nnz(Pass);
npl = nnz(Pass);
str = ['    All points: ',num2str(np),', Second filtering: ',num2str(nf),', Points left: ',num2str(npl)];
disp(str)

nf = np0-npl;
str = ['    All points: ',num2str(np0),', All filtered points: ',num2str(nf),', Points left: ',num2str(npl)];
disp(str)


%% Plot the filtered and unfiltered point clouds
comparison_plot(P(Pass,:),P,1,8)
%point_cloud_plotting(P,2,4)