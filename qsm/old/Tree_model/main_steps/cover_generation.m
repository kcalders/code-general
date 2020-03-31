function [Bal,Cen,Nei,Bale,Ball] = cover_generation(P,dmin,rcov,ncov)

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


x_min = min(P(:,1));
x_max = max(P(:,1));
y_min = min(P(:,2));
y_max = max(P(:,2));
z_min = min(P(:,3));
z_max = max(P(:,3));
dx = x_max-x_min;
dy = y_max-y_min;
dz = z_max-z_min;
d = max([dx dy dz]);
e = floor(log10(d));
r = 10^e/200;


%% Partition the point cloud into cubes
[partition,CC,~,info] = partition_of_point_cloud(P,r);

%disp([info(4)*info(7) info(5)*info(7) info(6)*info(7)])

np = length(P(:,1));
R = randperm(np);
s = max(100,ceil(np/1000));
%r = info(7);
N = zeros(s,5);
for j = 1:s
    q = r/0.75;
    points = partition(CC(R(j),1)-1:CC(R(j),1)+1,CC(R(j),2)-1:CC(R(j),2)+1,CC(R(j),3)-1:CC(R(j),3)+1);
    points = vertcat(points{:,:});
    cube = P(points,:);
    dist = (P(R(j),1)-cube(:,1)).^2+(P(R(j),2)-cube(:,2)).^2+(P(R(j),3)-cube(:,3)).^2;
    for i = 1:5
        %q = r*2^(-i+1);
        %q = r-(i-1)*10^(e-3);
        q = 0.75*q;
        J = (dist < q^2);
        dist = dist(J);
        N(j,i) = nnz(J);
    end
end
disp([r 0.75*r 0.75^2*r 0.75^3*r 0.75^4*r])
disp(max(N))
disp(round(mean(N)))
disp(min(N))

M = round(mean(N));
I = M > 5;
n = nnz(I)

R = randperm(np);
for j = 1:s
    q = r/0.75;
    points = partition(CC(R(j),1)-1:CC(R(j),1)+1,CC(R(j),2)-1:CC(R(j),2)+1,CC(R(j),3)-1:CC(R(j),3)+1);
    points = vertcat(points{:,:});
    cube = P(points,:);
    V = mat_vec_subtraction(cube,P(R(j),:));
    L = sum(V.*V,2);
    l = sqrt(L);
    v = [V(:,1)./l V(:,2)./l V(:,3)./l];
    dist = (P(R(j),1)-cube(:,1)).^2+(P(R(j),2)-cube(:,2)).^2+(P(R(j),3)-cube(:,3)).^2;
    for i = 1:n
        q = 0.75*q;
        J = (dist < q^2);
        cube = cube(J,:);
        dist = dist(J);
        V = V(J,:);
        v = v(J,:);
        l = l(J);
        
        [D,C] = dimensions(cube);
        A = V*[C(1:3)' C(4:6)'];
        K = atan2(A(:,1),A(:,2));
        k = acos(V*C(7:9)');
        disp(nnz(J))
        disp([q D])
        
        figure(2)
        subplot(3,1,1)
        plot(A(:,1),A(:,2),'.b','Markersize',10)
        axis equal
        subplot(3,1,2)
        plot(l,K,'.b','Markersize',10)
        subplot(3,1,3)
        plot(l,k,'.b','Markersize',10)
        
        figure(1)
        plot3(P(R(j),1),P(R(j),2),P(R(j),3),'.g','Markersize',20)
        hold on
        comparison_plot(cube,P(points,:),1,10)
        hold off
        pause
    end
end


%% Large balls and the centers
np = length(P(:,1));
R = randperm(np);   % use random permutation of points, results different covers for same inputs
%R = 1:1:np;
Ball = cell(np,1);
Cen = zeros(np,1);
NotInspected = true(np,1);
Dist = ones(np,1);
BallOfPoint = zeros(np,1);
t = 0;

for i = 1:4
for j = 1:np
    if NotInspected(R(j))
        points = partition(CC(R(j),1)-1:CC(R(j),1)+1,CC(R(j),2)-1:CC(R(j),2)+1,CC(R(j),3)-1:CC(R(j),3)+1);
        points = vertcat(points{:,:});
        cube = P(points,:);
        dist = (P(R(j),1)-cube(:,1)).^2+(P(R(j),2)-cube(:,2)).^2+(P(R(j),3)-cube(:,3)).^2;
        J = (dist < rcov^2);
        I = points(J);
        d = dist(J);
        if length(I) >= ncov
            J = (dist < dmin^2);
            NotInspected(points(J)) = false;
            t = t+1;
            Ball{t} = I;
            Cen(t) = R(j);
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
n = nnz(NotInspected);
str = ['    ',num2str(t),' cover sets, points not covered: ',num2str(n)];
disp(str)
