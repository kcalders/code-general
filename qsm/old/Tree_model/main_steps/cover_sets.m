function [Bal,Cen,Nei,Bale,Dim,Dir,Prin] = cover_sets(P,dmin,rcov,ncov)

% ---------------------------------------------------------------------
% COVER_SETS.M          Creates cover sets (surface patches), their
%                       neighbor-relations, and geometric characteristics 
%                       for a point cloud
%
% Version 1.0
% Author        Pasi Raumonen
% Created       14 June 2013
%
% This software cannot be distributed without the permission of P. Raumonen
% and it is only for non-commersial use. This restriction holds also for 
% the modules or subprograms of the software:
% PARTITION_OF_POINT_CLOUD.M, DIMENSIONS.M, DIMENSIONALITY.M, 
% OPTIMAL_ORTHOGONAL_VECTOR.M
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
% Bale      Extended cover sets, (n_sets x 1)-cell
% Dim       Box-dimensions and dimensionality values of the sets, (n_sets x 6)-vector
% Dir       Direction of the underlying branch at each cover set, (n_sets x 3)-vector
% Prin      Principal components of the sets, (n_sets x 9)-vector

%% Partition the point cloud into cubes
[partition,CC,CN] = partition_of_point_cloud(P,rcov);


%% Large balls and the centers
np = size(P,1);
R = randperm(np);   % use random permutation of points, results different covers for same input
%R = 1:1:np;
Ball = cell(np,1);  % the large balls which are used to generate the cover sets, their neighbors, and characteristics
Cen = zeros(np,1);  % the center points of the balls/cover sets
NotExa = true(np,1); % Points not yet examined
Dist = ones(np,1);  % distance of point to the closest center 
BoP = zeros(np,1);  % the balls/cover sets the points belong
Vec = zeros(np,3);  % the lines/vectors between points and the closest centers
Dim = zeros(np,6);  % box-dimensions and dimensionality values of the cover sets
Prin = zeros(np,9); % principal components of the cover sets
nc = 0;             % number of sets generated

Radius = rcov^2;
MaxDist = dmin^2;
if ~iscell(CN)
    % Uses normal fast approach when possible
    for i = 1:np
        if NotExa(R(i))
            points = partition(CC(R(i),1)-1:CC(R(i),1)+1,CC(R(i),2)-1:CC(R(i),2)+1,CC(R(i),3)-1:CC(R(i),3)+1);
            points = vertcat(points{:});
            V = mat_vec_subtraction(P(points,:),P(R(i),:));
            dist = sum(V.*V,2);
            J = dist < Radius;
            if nnz(J) >= ncov
                I = points(J);
                d = dist(J);
                V = V(J,:);
                J = (dist < MaxDist);
                NotExa(points(J)) = false;
                nc = nc+1;
                Ball{nc} = I;
                if length(I) > 5
                    %[D,W] = dimensionality(P(I,:));
                    [D,W] = dimensions(P(I,:));
                    Dim(nc,:) = D;
                    Prin(nc,:) = W;
                end
                Cen(nc) = R(i);
                D = Dist(I);
                L = d < D;
                I = I(L);
                Dist(I) = d(L);
                BoP(I) = nc;
                Vec(I,:) = V(L,:);
            end
        end
    end
    
else
    % Uses "sparse" approach because the memory requirements of the normal 
    % approach
    for i = 1:np
        if NotExa(R(i))
            N = CN{CC(R(i))};
            points = vertcat(partition{N});
            V = mat_vec_subtraction(P(points,:),P(R(i),:));
            dist = sum(V.*V,2);
            J = dist < Radius;
            if nnz(J) >= ncov
                I = points(J);
                d = dist(J);
                V = V(J,:);
                J = (dist < MaxDist);
                NotExa(points(J)) = false;
                nc = nc+1;
                Ball{nc} = I;
                if length(I) > 5
                    %[D,W] = dimensionality(P(I,:));
                    [D,W] = dimensions(P(I,:));
                    Dim(nc,:) = D;
                    Prin(nc,:) = W;
                end
                Cen(nc) = R(i);
                D = Dist(I);
                L = d < D;
                I = I(L);
                Dist(I) = d(L);
                BoP(I) = nc;
                Vec(I,:) = V(L,:);
            end
        end
    end
end
Ball = Ball(1:nc,:);
Cen = Cen(1:nc);
Dim = Dim(1:nc,:);
Prin = Prin(1:nc,:);


%% Number of points in each ball and index of each point in its ball
Num = zeros(nc,1);
Ind = zeros(np,1);
for i = 1:np
    if BoP(i) > 0
        Num(BoP(i)) = Num(BoP(i))+1;
        Ind(i) = Num(BoP(i));
    end
end


%% Cover sets
% Initialization of the "Bal"
Bal = cell(nc,1);
for i = 1:nc
    Bal{i} = zeros(Num(i),1);
end

% Define the "Bal"
for i = 1:np
    if BoP(i) > 0
        Bal{BoP(i),1}(Ind(i)) = i;
    end
end


%% Neighbors and extended cover sets
% Define neighbors. 
% Simple rule: Sets A and B are neighbors if the large ball of A 
% contains points of B. 
% However, this may contain neighbor's neighbors or large gaps between sets
% and these cases are removed.
% Notice that this is not a symmetric relation.
Nei = cell(nc,1);
for i = 1:nc
    B = Ball{i};        % the points in the big ball of cover set "i" 
    I = (BoP(B) ~= i);  
    N = B(I);           % the points of B not in the cover set "i"
    N = unique(BoP(N)); % the neighboring cover sets, "simple rule"
    nb = length(N);     % number of possible neighbors
    if nb > 1
        % if there are more than one possible neighbor, check which are
        % "true" neighbors and not neighbor's neighbors
        V = Vec(Bal{i},:);    % "lines" of the points in the cover set "i"
        T = zeros(nb,3);      % unit direction vectors of "N" as seen from "i"
        D = zeros(nb,2);      % Distances between "i" and "N" (1. column)
                          % 2. column: sum of extents of "i" and "N" in the T-lines
                          % the width of the gaps between "i" and "N" are D(:,1)-D(:,2) 
        for j = 1:nb    % determine T and D for each possible neighbor
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
        Nei{i} = N(K);
    else
        Nei{i} = N; 
    end
end

% Make the relation symmetric by adding, if needed, A as B's neighbor 
% in the case B is A's neighbor
for i = 1:nc
    N = Nei{i};
    for j = 1:length(N)
        K = (Nei{N(j)} == i);
        if ~any(K)
            Nei{N(j)} = [Nei{N(j)}; i];
        end
    end
end

if nargout > 3
    %% Extended cover sets and directions
    Bale = cell(nc,1);
    Dir = zeros(nc,3);
    I = true(nc,1);
    for i = 1:nc
        Bale{i} = [i; Nei{i}];
        
        if Dim(i,1) > 0.85 % elongated sets
            Dir(i,:) = Prin(i,1:3);
            I(i) = false;
        elseif Dim(i,2) < 0.9
            Dir(i,:) = optimal_orthogonal_vector(Prin(Bale{i},7:9));
            I(i) = false;
        end
    end
    for i = 1:nc
        if I(i)
            % very planar
            B = unique(vertcat(Bale{Bale{i}}));
            Dir(i,:) = optimal_orthogonal_vector(Prin(B,7:9));
        end
    end
end

%% Display statistics
str = ['    ',num2str(nc),' cover sets, points not covered: ',num2str(np-nnz(BoP))];
disp(str)
