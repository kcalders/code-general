function [Bal,Nei,Cen,BoP] = cover_sets(P,dmin,rcov,ncov,RS)

% ---------------------------------------------------------------------
% COVER_SETS.M          Creates cover sets (surface patches), their
%                       neighbor-relations, and geometric characteristics 
%                       for a point cloud
%
% Version 1.3
% Author        Pasi Raumonen
% Created       21 February 2014
%
% This software cannot be distributed without the permission of P. Raumonen
% and it is only for non-commersial use. This restriction holds also for 
% the modules or subprograms of the software:
% PARTITION_OF_POINT_CLOUD.M
% ---------------------------------------------------------------------

% Covers the point cloud with small sets, which are along the surface, 
% such that each point belongs at most one cover set; i.e. the cover is 
% a partition of the point cloud. 
% 
% The cover is generated such that at first the point cloud is covered 
% with balls with radius RCOV. This first cover is such that 
% 1) the minimum distance between the centers is DMIN, and 
% 2) the maximum distance from any point to nearest center is also DMIN.
% Then the first cover of RCOV-balls is used to define a second cover:
% each RCOV-ball A defines corresponding cover set B in the second cover
% such that B contains those points of A that are nearer to the center of
% A than any other center of RCOV-balls. The RCOV-balls also define 
% the neighbors for the second cover: Let CA and CB denote cover sets in 
% the second cover, and BA and BB their RCOV-balls. Then CB is 
% a neighbor of CA, and vice versa, if BA and CB intersect or 
% BC and CA intersect. However if the gap between CA and CB is large, then
% thei are not neighbors.
%
% Inputs: 
% P         Point cloud
% dmin      Minimum distance between centers of cover sets; i.e. the
%               minimum diameter of a cover set
% rcov      Radius of the balls used to generate the cover sets, these 
%               balls are also used to determine the neighbors and the 
%               cover set characteristics              
% nmin      Minimum number of points in a rcov-ball
%
% Outputs:
% Bal       Cover sets, (n_sets x 1)-cell
% Cen       Center points of the cover sets, (n_sets x 1)-vector
% Nei       Neighboring cover sets of each cover set, (n_sets x 1)-cell


%% Large balls and the centers
np = size(P,1);
s = RandStream.create('mt19937ar','seed',sum(100*clock));
RandStream.setGlobalStream(s);
R = int32(randperm(np)); % random permutation of points, results different covers for same input
Ball = cell(np,1);  % the large balls which are used to generate the cover sets, their neighbors, and characteristics
Cen = zeros(np,1,'int32');  % the center points of the balls/cover sets
NotExa = true(np,1); % Points not yet examined
Dist = 1e8*ones(np,1,'single');  % distance of point to the closest center 
BoP = zeros(np,1,'int32');  % the balls/cover sets the points belong
Vec = zeros(np,3,'single');  % the lines/vectors between points and the closest centers
nb = 0;             % number of sets generated

if nargin == 4
    % normal approach, size of the cover sets is same everywhere
    % Partition the point cloud into cubes
    [partition,CC] = partition_of_point_cloud(P,rcov);
    
    Radius = rcov^2;
    MaxDist = dmin^2;
    for i = 1:np
        if NotExa(R(i))
            Q = R(i);
            points = partition(CC(Q,1)-1:CC(Q,1)+1,CC(Q,2)-1:CC(Q,2)+1,CC(Q,3)-1:CC(Q,3)+1);
            points = vertcat(points{:});
            V = [P(points,1)-P(Q,1) P(points,2)-P(Q,2) P(points,3)-P(Q,3)];
            dist = sum(V.*V,2);
            J = dist < Radius;
            if nnz(J) >= ncov
                I = points(J);
                d = dist(J);
                V = V(J,:);
                J = (dist < MaxDist);
                NotExa(points(J)) = false;
                nb = nb+1;
                Ball{nb} = I;
                Cen(nb) = Q;
                D = Dist(I);
                L = d < D;
                I = I(L);
                Dist(I) = d(L);
                BoP(I) = nb;
                Vec(I,:) = V(L,:);
            end
        end
    end
    Ball = Ball(1:nb,:);
    Cen = Cen(1:nb);
else
    % Use relative sizes, in which case the size varies
    % Partition the point cloud into cubes
    r = (double(min(RS))/256*0.5+0.5)*rcov+1e-4;
    [partition,CC] = partition_of_point_cloud(P,r);
    
    e = rcov-dmin;
    for i = 1:np
        if NotExa(R(i))
            Q = R(i);
            rs = double(RS(i))/256*0.5+0.5;
            Dmin = dmin*rs;
            Radius = Dmin+e;
            a = ceil(Radius/r);
            points = partition(CC(Q,1)-a:CC(Q,1)+a,CC(Q,2)-a:CC(Q,2)+a,CC(Q,3)-a:CC(Q,3)+a);
            points = vertcat(points{:});
            V = [P(points,1)-P(Q,1) P(points,2)-P(Q,2) P(points,3)-P(Q,3)];
            dist = sum(V.*V,2);
            J = dist < Radius^2;
            if nnz(J) >= ncov
                I = points(J);
                d = dist(J);
                V = V(J,:);
                J = (dist < Dmin^2);
                NotExa(points(J)) = false;
                nb = nb+1;
                Ball{nb} = I;
                Cen(nb) = Q;
                D = Dist(I);
                L = d < D;
                I = I(L);
                Dist(I) = d(L);
                BoP(I) = nb;
                Vec(I,:) = V(L,:);
            end
        end
    end
    Ball = Ball(1:nb,:);
    Cen = Cen(1:nb);
end


%% Cover sets
% Number of points in each ball and index of each point in its ball
Num = zeros(nb,1,'int32');
Ind = zeros(np,1,'int32');
for i = 1:np
    if BoP(i) > 0
        Num(BoP(i)) = Num(BoP(i))+1;
        Ind(i) = Num(BoP(i));
    end
end

% Initialization of the "Bal"
Bal = cell(nb,1);
for i = 1:nb
    Bal{i} = zeros(Num(i),1,'int32');
end

% Define the "Bal"
for i = 1:np
    if BoP(i) > 0
        Bal{BoP(i),1}(Ind(i)) = i;
    end
end


%% Neighbors
% Define neighbors. 
% Simple rule: Sets A and B are neighbors if the large ball of A 
% contains points of B. 
% However, this may contain neighbor's neighbors or large gaps between sets
% and these cases are removed.
% Notice that this is not a symmetric relation.
Nei = cell(nb,1);
for i = 1:nb
    B = Ball{i};        % the points in the big ball of cover set "i" 
    I = (BoP(B) ~= i);  
    N = B(I);           % the points of B not in the cover set "i"
    N = unique(BoP(N)); % the neighboring cover sets, "simple rule"
    n = length(N);     % number of possible neighbors
    if n > 1
        % if there are more than one possible neighbor, check which are
        % "true" neighbors and not neighbor's neighbors
        V = Vec(Bal{i},:);    % "lines" of the points in the cover set "i"
        T = zeros(n,3);      % unit direction vectors of "N" as seen from "i"
        D = zeros(n,2);      % Distances between "i" and "N" (1. column)
                          % 2. column: sum of extents of "i" and "N" in the T-lines
                          % the width of the gaps between "i" and "N" are D(:,1)-D(:,2) 
        for j = 1:n    % determine T and D for each possible neighbor
            W = Vec(Bal{N(j)},:);  % "lines" of ball "n(j)"
            U = P(Cen(N(j)),:)-P(Cen(i),:); % line joining center's of "i" and "N(j)"
            L = norm(U);   % distance between the cover sets "i" and "N(j)"
            U = U/L;       % normalize
            T(j,:) = U;
            d1 = V*U';     % projected distances of points of "i"
            d2 = W*U';     % projected distances of points of "N(j)"
            D(j,:) = [L max(d1)+abs(min(d2))];  % distance and sum of T-extents
        end
        I = D(:,2) > 0.8*D(:,1);  % small gaps compared to the distances
        J = D(:,1)-D(:,2) < dmin/2; % small gaps, under half of the minimum set diameter
        K = I|J;   % the neighbors, sets with small gaps to the "i" 
        if ~any(I) % if no small gaps, try larger gaps
            I = D(:,2) > 0.6*D(:,1);
            J = D(:,1)-D(:,2) < rcov/2;
            K = I|J;
        end
        if any(K) && ~all(K) && nnz(K) <= 2
            % if there are 1-2 neighbors, check the ones at the opposite
            % direction of the closest neighbor and add the closest if close enough 
            [~,J] = min(D(:,1));  % N(J) is the closest neighbor
            a = T*T(J,:)';    % cosines of the T-vectors compared to T(J)-vector
            J = a < -0.3;     % sets that are opposite direction
            if ~any(K(J)) && ~isempty(K(J))
                % opposite sets not yet neighbors exists
                D(~J,1) = 10*rcov;
                [~,J] = min(D(:,1));  % the closest opposite set
                if D(J,1)-D(J,2) < rcov
                    % accept the closest opposite set if the distance is
                    % smaller than rcov, the maximum distance of a possible
                    % neighbor
                    K(J) = true;
                end
            end
        end
        Nei{i} = int32(N(K));
    else
        Nei{i} = int32(N); 
    end
end

% Make the relation symmetric by adding, if needed, A as B's neighbor 
% in the case B is A's neighbor
for i = 1:nb
    N = Nei{i};
    for j = 1:length(N)
        K = (Nei{N(j)} == i);
        if ~any(K)
            Nei{N(j)} = int32([Nei{N(j)}; i]);
        end
    end
end


%% Display statistics
str = ['    ',num2str(nb),' cover sets, points not covered: ',num2str(np-nnz(BoP))];
disp(str)
