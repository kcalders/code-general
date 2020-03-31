function [Bal,Cen,Nei,Bale,Ball] = cover_sets_and_neighbors(P,dmin,rcov,ncov)

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
% the and second cover, and BA and BB their RCOV-balls. Then CB is 
% a neighbor of CA, and vice versa, if BA and CB intersect or 
% BC and CA intersect.
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
% Cen       Center points of the cover sets, (n_sets x 1)-vectot
% Nei       Neighboring cover sets of each cover set, (n_sets x 1)-cell
% Bale      Extended cover sets, (n_sets x 1)-cell
% Ball      rcov-balls generated at first


%% Partition the point cloud into cubes
[partition,CC,CN] = partition_of_point_cloud(P,rcov);


%% Large balls and the centers
np = length(P(:,1));
%R = randperm(np);   % use random permutation of points, results different covers for same inputs
R = 1:1:np;
Ball = cell(np,1);
Cen = zeros(np,1);
NotInspected = true(np,1);
Dist = ones(np,1);
BallOfPoint = zeros(np,1);
t = 0;

if ~iscell(CN)
    % Uses normal fast approach when possible
    for i = 1:np
        if NotInspected(R(i))
            points = partition(CC(R(i),1)-1:CC(R(i),1)+1,CC(R(i),2)-1:CC(R(i),2)+1,CC(R(i),3)-1:CC(R(i),3)+1);
            points = vertcat(points{:});
            cube = P(points,:);
            dist = (P(R(i),1)-cube(:,1)).^2+(P(R(i),2)-cube(:,2)).^2+(P(R(i),3)-cube(:,3)).^2;
            J = (dist < rcov^2);
            I = points(J);
            d = dist(J);
            if length(I) >= ncov
                J = (dist < dmin^2);
                NotInspected(points(J)) = false;
                t = t+1;
                Ball{t} = I;
                Cen(t) = R(i);
                D = Dist(I);
                L = d < D;
                Dist(I(L)) = d(L);
                BallOfPoint(I(L)) = t;
            end
        end
    end
    
else
    % Uses "sparse" approach because the memory requirements of the normal 
    % approach
    disp('Sparce approach')
    for i = 1:np
        if NotInspected(R(i))
            N = CN{CC(R(i))};
            points = vertcat(partition{N});
            %points = vertcat(points{:,:});
            cube = P(points,:);
            dist = (P(R(i),1)-cube(:,1)).^2+(P(R(i),2)-cube(:,2)).^2+(P(R(i),3)-cube(:,3)).^2;
            J = (dist < rcov^2);
            I = points(J);
            d = dist(J);
            if length(I) >= ncov
                J = (dist < dmin^2);
                NotInspected(points(J)) = false;
                t = t+1;
                Ball{t} = I;
                Cen(t) = R(i);
                D = Dist(I);
                L = d < D;
                Dist(I(L)) = d(L);
                BallOfPoint(I(L)) = t;
            end
        end
    end
end
Ball = Ball(1:t,:);
Cen = Cen(1:t);


%% Number of points in each ball and index of each point in its ball
Num = zeros(t,1);
Ind = zeros(np,1);
for i = 1:np
    if BallOfPoint(i) > 0
        Num(BallOfPoint(i)) = Num(BallOfPoint(i))+1;
        Ind(i) = Num(BallOfPoint(i));
    end
end


%% Cover sets
% Initialization of Bal
Bal = cell(t,1);
for i = 1:t
    Bal{i} = zeros(Num(i),1);
end

% Define the Bal
for i = 1:np
    if BallOfPoint(i) > 0
        Bal{BallOfPoint(i),1}(Ind(i)) = i;
    end
end


%% Neighbors and extended cover sets
% Define neighbors. Sets A and B are neighbors if the large ball of A 
% contains points of B. This is not a symmetric relation.
Nei = cell(t,1);
for i = 1:t
    B = Ball{i};
    I = (BallOfPoint(B) ~= i);
    n = B(I);
    N = BallOfPoint(n);
    Nei{i} = unique(N);
end

% Make the relation symmetric by adding, if needed, A as B's neighbor 
% in the case B is A's neighbor
for i = 1:t
    N = Nei{i};
    for j = 1:length(N)
        K = (Nei{N(j)} == i);
        if ~any(K)
            Nei{N(j)} = [Nei{N(j)}; i];
        end
    end
end

if nargout > 3
    %% Extended cover sets
    Bale = cell(t,1);
    for i = 1:t
        Bale{i} = [i; Nei{i}];
    end
end


%% Display statistics
%n = nnz(NotInspected);
n = nnz(BallOfPoint);
str = ['    ',num2str(t),' cover sets, points not covered: ',num2str(np-n)];
disp(str)
