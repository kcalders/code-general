function [Base, Forb, Bale, Nei] = tree_sets(P, Cen, Bale, Nei, Dir, Dim, Prin, dmin, Bal)

% ---------------------------------------------------------------------
% TREE_SETS.M       Determines the base of the trunk and the cover sets 
%                   belonging to the tree, updates the neighbor-relation
%
% Version 1.0
% Author        Pasi Raumonen
% Created       14 June 2013
%
% This software cannot be distributed without the permission of P. Raumonen
% and it is only for non-commercial use. These restrictions holds also for 
% the modules or subprograms of the software:
% CONNECTED_COMPONENTS.M, PARTITION_OF_POINT_CLOUD.M.
% ---------------------------------------------------------------------

% Determines the cover sets that belong to the tree. Determines also the
% base of the tree and updates the neighbor-relation such that all of the
% tree is connected, i.e., the cover sets belonging to the tree form a
% single connected component.

% Inputs:
% P         Point cloud
% Cen       Center points of the cover sets, (n_sets x 1)-vector
% Bale      Extended cover sets, (n_sets x 1)-cell
% Nei       Neighboring cover sets, (n_sets x 1)-cell
% Dir       Direction vectors of the cover sets (n_sets x 3)-maxtrix
% Dim       Dimensionality values of the sets, (n_sets x 3)-matrix
% dmin      Minimum diameter of the cover sets
%
% Outputs:
% Base      Base of the trunk (vector)
% Forn      Cover sets not part of the tree
% Bale      Updated extended cover sets
% Nei       Updated neigbors


nb = max(size(Bale));   % number of cover sets
Ind = (1:1:nb)';
Fal = false(nb,1);
Bottom = min(P(Cen,3));
[Top,top] = max(P(Cen,3));
Height = Top-Bottom;


%% Initial classification of possible trunk cover sets
% Cover sets that are 1) fairly planar, and 2) rather parallel to the trunk
%X = cov(P(Cen,:));
%[U,~,~] = svd(X);
%ref = U(:,1)  % approximate direction of the trunk
ref = [0 0 1]';
%ref = [-1 0 1]';
ref = ref/norm(ref);
Par = abs(Dir*ref);

a = 0.2;
b = 1.8*dmin;
c = 0.95;
I = (Prin(:,9) < a);  % normals are nearly horizontal
J = (Dim(:,2) > b);  % sets are wide enough
K = Par > c;  % directions are nearly vertical
trunk = I&J&K;
trunk(vertcat(Bale{trunk})) = true; % expand "trunk" to make it more connnected
trunk(vertcat(Bale{trunk})) = true;
trunk(vertcat(Bale{trunk})) = true;
I = P(Cen,3)-Bottom < 0.3*Height; % define bottom region "bot"
J = P(Cen,3)-Bottom > 0.1*Height;
bot = I&J;
N = nnz(bot);
%plot_segments(P,Bal,3,Ind(trunk),Ind)
%disp([nnz(I) nnz(J) nnz(K) nnz(trunk)])
%pause
while c > 0.6 && nnz(trunk & bot) < 0.66*N
    b = 0.9*b;
    c = c-0.05;
    J = (Dim(:,2) > b);
    K = Par > c;
    trunk = J & K;
    trunk(vertcat(Bale{trunk})) = true;
    trunk(vertcat(Bale{trunk})) = true;
    trunk(vertcat(Bale{trunk})) = true;
    %plot_segments(P,Bal,3,Ind(trunk),Ind)
    %disp([nnz(J) nnz(K) nnz(trunk)])
    %pause
end

%plot_segments(P,Bal,3,Ind(trunk),Ind)
%pause


%% Redefine the classification of trunk sets
% Remove small parts of branches included in the initial classification by 
% including only the largest components of the trunk sets.
if Height > 20
    B = 0.5; % The base of the trunk/tree (50 cm slice at the bottom of the trunk)
    MinComp = ceil(0.5/dmin)^2;  % minimum number of cover sets in a trunk component
elseif Height > 5
    B = 0.5;
    MinComp = ceil(0.3/dmin)^2;
else
    B = 0.1;
    MinComp = ceil(0.1/dmin)^2;
end
[trunk,compsize] = connected_components(Nei,trunk,MinComp,Fal);

T = vertcat(trunk{:});
T = mean(P(Cen(T),:));

%plot_segments(P,Bal,2,vertcat(trunk{:}),1:1:nb)
%plot_segs(P,trunk,4,Bal)
%pause

% Find the lowest component
n = length(compsize);
CompHeight = zeros(n,1);
Dist = zeros(n,1);
for i = 1:n
    d = mat_vec_subtraction(P(Cen(trunk{i}),1:2),T(1:2));
    Dist(i) = min(sqrt(sum(d.*d,2)));
    CompHeight(i) = min(P(Cen(trunk{i}),3));
end
I = Dist < 0.5;
trunk = trunk(I);
CompHeight = CompHeight(I,:);
n = nnz(I);
[hmin,I] = min(CompHeight);

% Expand components and define Trunk again from these components
h = min(0.05,2*dmin);
for i = 1:n
    if i ~= I
        C = trunk{i};
        C = vertcat(Bale{C});
        C = vertcat(Bale{C});
        trunk{i} = unique(C);
    else
        C = trunk{i};
        D = vertcat(Bale{C});
        D = vertcat(Bale{D});
        J = P(Cen(D),3) >= hmin+h; 
        trunk{i} = unique(D(J));
    end
end
Trunk = Fal;
Trunk(vertcat(trunk{:})) = true;

%plot_segs(P,trunk,4,Bal)
%pause


% The base of the trunk/tree
Foot = trunk{I(1)};  % The lowest large trunk component
Bottom = min(P(Cen(Foot),3));
I = (P(Cen(Foot),3) < Bottom+B);
Base = Foot(I);    % The cover sets forming the base of the trunk
%plot_segments(P,Bal,3,Base,Foot,Ind(Trunk),Ind)
%pause




%% The ground
% The ground component connected to the foot of the trunk
GL = min(P(Cen(Base),3));   % the ground level
G = setdiff(vertcat(Bale{Base}),Foot);
I = P(Cen(G),3) < GL+2*dmin;
G = G(I);
if ~isempty(G)
    Forb = Fal;
    Forb(G) = true;
    Forb(Base) = false;
    Added = Forb;
    while any(Added)
        Added(vertcat(Bale{Added})) = true;
        Added(Forb) = false;
        Added(Trunk) = false;
        Added(Base) = false;
        Forb(Added) = true;
    end
else
    Forb = Fal;
end

%J = P(Cen,3) < GL+20*dmin;
%Forb(J) = true;
%Forb(Base) = false;
%Forb(Trunk) = false;

%plot_segments(P,Bal,1,Ind(Forb),Base,Ind(Trunk),Ind)
%pause


%% Update the neighbor-relation near the Trunk
% Select only the cover sets nearest to the "Trunk"
T = Ind(Trunk);
n = length(T);
Fast = true;
EdgeLength = min(0.15,5*dmin);
[Par,CC,~,info] = partition_of_point_cloud(P(Cen,:),EdgeLength,Fast);

I = false(info(4),info(5),info(6)); % cover sets near the "Trunk"
for j = 1:n
    J = T(j);
    I(CC(J,1)-1:CC(J,1)+1,CC(J,2)-1:CC(J,2)+1,CC(J,3)-1:CC(J,3)+1) = true;
end
N = Fal;
N(vertcat(Par{I})) = true;
N(T) = false;   % Only the non-trunk cover sets
N = Ind(N);


% Determine separate components of "N", cover sets near the trunk
[Comps,CompSize] = connected_components(Nei,N,1,Fal);
nc = length(CompSize);

% Expand each component
for i = 1:nc
    C = Comps{i};
    C = unique(vertcat(Bale{C}));
    C = unique(vertcat(Bale{C}));
    C = unique(vertcat(Bale{C}));
    C = unique(vertcat(Bale{C}));
    C = unique(vertcat(Bale{C}));
    Comps{i} = C;
end


N = Fal;
N(vertcat(Comps{:})) = true;
N(T) = false;    % Only the non-trunk cover sets
N = Ind(N); 

%plot_segments(P,Bal,1,Ind(Trunk),N,Ind(Forb))
%pause

% Determine separate components of "N", cover sets near the trunk
[Comps,CompSize] = connected_components(Nei,N,1,Fal);
nc = length(CompSize);

% Check the components and possibly update the neighbor-relation
K = false(nc,1);
for i = 1:nc
    C = Comps{i};
    J = Trunk(C);
    L = Trunk(vertcat(Nei{C}));
    
    if any(J) || any(L)
        %disp('yhteydessa')
    else
        NC = length(C);
        % Select only the cover sets the nearest to the component
        I = false(info(4),info(5),info(6));
        for j = 1:NC
            J = C(j);
            I(CC(J,1)-1:CC(J,1)+1,CC(J,2)-1:CC(J,2)+1,CC(J,3)-1:CC(J,3)+1) = true;
        end
        
        M = Fal;
        M(vertcat(Par{I})) = true;
        M(~Trunk) = false; % The nearest "Trunk" cover sets
        M = Ind(M);
        
        d = pdist2(P(Cen(C),:),P(Cen(M),:));
        if NC == 1 && length(M) == 1
            dt = d;
            I = C;
            J = M;
        elseif NC == 1
            [dt,J] = min(d);
            I = C;
            J = M(J);
        elseif length(M) == 1
            [dt,I] = min(d);
            J = M;
            I = C(I);
        else
            [d,I] = min(d);
            [dt,J] = min(d);
            I = I(J);
            I = C(I);
            J = M(J);
        end
        
        if (NC > 5) || ((NC <= 5) && (dt < 3*dmin))
            
            K(i) = true;
            Nei{I} = [Nei{I}; J];
            Nei{J} = [Nei{J}; I];
            Bale{I} = [I; Nei{I}];
            Bale{J} = [J; Nei{J}];
        end
        
    end
end


%% Expand Trunk as much as possible
Trunk(Forb) = false;
Trunk(Base) = false;
Added = Trunk;
while any(Added)
    Added(vertcat(Bale{Added})) = true;
    Added(Trunk) = false;
    Added(Forb) = false;
    Added(Base) = false;
    Trunk(Added) = true;
end

%plot_segments(P,Bal,1,Base,Ind(Forb),Ind(Trunk),Ind)
%pause


%% Update the neighbor-relation near the expanded Trunk
% Select only the cover sets nearest to the "Trunk"
T = Ind(Trunk);
n = length(T);
Fast = true;
EdgeLength = min(0.1,5*dmin);
[Par,CC,~,info] = partition_of_point_cloud(P(Cen,:),EdgeLength,Fast);

I = false(info(4),info(5),info(6)); % cover sets near the "Trunk"
for j = 1:n
    J = T(j);
    I(CC(J,1)-1:CC(J,1)+1,CC(J,2)-1:CC(J,2)+1,CC(J,3)-1:CC(J,3)+1) = true;
end
N = Fal;
N(vertcat(Par{I})) = true;
N(T) = false;   % Only the non-trunk cover sets
N = Ind(N);


% Determine separate components of "N", cover sets near the trunk
[Comps,CompSize] = connected_components(Nei,N,1,Fal);
nc = length(CompSize);

% Expand each component
for i = 1:nc
    C = Comps{i};
    C = unique(vertcat(Bale{C}));
    C = unique(vertcat(Bale{C}));
    C = unique(vertcat(Bale{C}));
    C = unique(vertcat(Bale{C}));
    C = unique(vertcat(Bale{C}));
    Comps{i} = C;
end


N = Fal;
N(vertcat(Comps{:})) = true;
N(T) = false;    % Only the non-trunk cover sets
N = Ind(N); 

%plot_segments(P,Bal,1,Ind(Trunk),N)
%pause

% Determine separate components of "N", cover sets near the trunk
[Comps,CompSize] = connected_components(Nei,N,1,Fal);
nc = length(CompSize);

% Check the components and possibly update the neighbor-relation
K = false(nc,1);
for i = 1:nc
    C = Comps{i};
    J = Trunk(C);
    L = Trunk(vertcat(Nei{C}));
    
    if any(J) || any(L)
        %disp('yhteydessa')
    else
        NC = length(C);
        % Select only the cover sets the nearest to the component
        I = false(info(4),info(5),info(6));
        for j = 1:NC
            J = C(j);
            I(CC(J,1)-1:CC(J,1)+1,CC(J,2)-1:CC(J,2)+1,CC(J,3)-1:CC(J,3)+1) = true;
        end
        
        M = Fal;
        M(vertcat(Par{I})) = true;
        M(~Trunk) = false; % The nearest "Trunk" cover sets
        M = Ind(M);
        
        d = pdist2(P(Cen(C),:),P(Cen(M),:));
        %disp(size(d))
        if NC == 1 && length(M) == 1
            dt = d;
            I = C;
            J = M;
        elseif NC == 1
            [dt,J] = min(d);
            I = C;
            J = M(J);
        elseif length(M) == 1
            [dt,I] = min(d);
            J = M;
            I = C(I);
        else
            [d,I] = min(d);
            [dt,J] = min(d);
            I = I(J);
            I = C(I);
            J = M(J);
        end
        
        if (NC > 5) || ((NC <= 5) && (dt < 3*dmin))
            
            K(i) = true;
            Nei{I} = [Nei{I}; J];
            Nei{J} = [Nei{J}; I];
            Bale{I} = [I; Nei{I}];
            Bale{J} = [J; Nei{J}];
        end
        
    end
end



%% Determine the components not connected to the "Tree" and update the neighbor-relation

% Define the expanded trunk again
Trunk(Forb) = false;
Trunk(Base) = false;
Added = Trunk;
while any(Added)
    Added(vertcat(Bale{Added})) = true;
    Added(Trunk) = false;
    Added(Forb) = false;
    Added(Base) = false;
    Trunk(Added) = true;
end

% Define "Tree", the component starting from the base
Tree = Fal;
Tree(vertcat(Bale{Base})) = true;
Tree(Forb) = false;
Tree(Base) = false;
Added = Tree;
while any(Added)
    Added(vertcat(Bale{Added})) = true;
    Added(Tree) = false;
    Added(Forb) = false;
    Added(Base) = false;
    Tree(Added) = true;
end
%plot_segments(P,Bal,2,Ind(Forb),Ind(Tree),Ind(Trunk))
%plot_segments(P,Bal,3,Ind(Forb),Ind(Tree),Ind(Trunk),Ind)
%pause


Other = ~Fal;
Tree(Base) = false;
Other(Tree) = false;
Other(Forb) = false;
k0 = min(10,ceil(0.4/dmin));
k = k0;
Cmin = ceil(0.05/dmin);
%Cmin = 0;
dmax = 0;
W = 1000;
while any(Other)
    npre = nnz(Other);
    again = true;
    
    % Partition the centers of the cover sets into cubes
    Fast = true;
    [Par,CC,~,info] = partition_of_point_cloud(P(Cen,:),k*dmin,Fast);
    while any(Other) && again
        
        % Determine the components of "Other"
        [Comps,CompSize] = connected_components(Nei,Other,1,Fal);
        
        nc = length(CompSize);
        
        % Label the cover sets by their components
        CL = zeros(nb,1);
        CL(Trunk) = nc+1;
        CL(~Trunk&Tree) = nc+2;
        CL(Forb) = nc+3;
        for i = 1:nc
            CL(Comps{i}) = i;
        end
        
        % Check each component: part of "Tree", other components, or "Forb"
        for i = 1:nc
            C = Comps{i};
            NC = length(C);
            
            % Select only the cover sets the nearest to the component
            I = false(info(4),info(5),info(6));
            for j = 1:NC
                I(CC(C(j),1)-1:CC(C(j),1)+1,CC(C(j),2)-1:CC(C(j),2)+1,CC(C(j),3)-1:CC(C(j),3)+1) = true;
            end
            N = Fal;
            N(vertcat(Par{I})) = true;
            N(C) = false;  % The nearest cover sets
            N = Ind(N);

            L = CL(N);  % The component labels of the nearest cover sets
            Tru = L == nc+1; % "Trunk" sets
            Tre = L == nc+2;  % "Tree" sets
            F = L == nc+3;  % "Forb" sets
            
            O = N(~(Tre|F|Tru));
            Tru = N(Tru);
            Tre = N(Tre);
            F = N(F);
            
            % Determine the closest sets for "Other", "Trunk", "Tree" and "Forb"
            if ~isempty(O)
                if NC > W
                    m = ceil(NC/W);
                    dmin = 100;
                    for j = 1:m
                        if j < m
                            d = pdist2(P(Cen(C((j-1)*W+1:j*W)),:),P(Cen(O),:));
                            D = min(min(d));
                            if D < dmin
                                dmin = D;
                                dis = d;
                            end
                        else
                            d = pdist2(P(Cen(C((j-1)*W+1:end)),:),P(Cen(O),:));
                            D = min(min(d));
                            if D < dmin
                                dmin = D;
                                dis = d;
                            end
                        end
                    end
                    d = dis;
                else
                    d = pdist2(P(Cen(C),:),P(Cen(O),:));
                end
                if NC == 1 && length(O) == 1
                    do = d;
                    IO = 1;
                    JO = 1;
                elseif NC == 1
                    [do,JO] = min(d);
                    IO = 1;
                elseif length(O) == 1
                    [do,IO] = min(d);
                    JO = 1;
                else
                    [d,IO] = min(d);
                    [do,JO] = min(d);
                    IO = IO(JO);
                end
            else
                do = 9;
            end
            if ~isempty(Tru) 
                if NC > W
                    m = ceil(NC/W);
                    dmin = 100;
                    for j = 1:m
                        if j < m
                            d = pdist2(P(Cen(C((j-1)*W+1:j*W)),:),P(Cen(Tru),:));
                            D = min(min(d));
                            if D < dmin
                                dmin = D;
                                dis = d;
                                CTu = C((j-1)*W+1:j*W);
                            end
                        else
                            d = pdist2(P(Cen(C((j-1)*W+1:end)),:),P(Cen(Tru),:));
                            D = min(min(d));
                            if D < dmin
                                dmin = D;
                                dis = d;
                                CTu = C((j-1)*W+1:end);
                            end
                        end
                    end
                    d = dis;
                else
                    d = pdist2(P(Cen(C),:),P(Cen(Tru),:));
                end
                if NC == 1 && length(Tru) == 1
                    du = d;
                    IU = 1;
                    JU = 1;
                elseif NC == 1
                    [du,JU] = min(d);
                    IU = 1;
                elseif length(Tru) == 1
                    [du,IU] = min(d);
                    JU = 1;
                else
                    [d,IU] = min(d);
                    [du,JU] = min(d);
                    IU = IU(JU);
                end
            else
                du = 7;
            end
            if ~isempty(Tre)
                if NC > W
                    m = ceil(NC/W);
                    dmin = 100;
                    for j = 1:m
                        if j < m
                            d = pdist2(P(Cen(C((j-1)*W+1:j*W)),:),P(Cen(Tre),:));
                            D = min(min(d));
                            if D < dmin
                                dmin = D;
                                dis = d;
                                CTe = C((j-1)*W+1:j*W);
                            end
                        else
                            d = pdist2(P(Cen(C((j-1)*W+1:end)),:),P(Cen(Tre),:));
                            D = min(min(d));
                            if D < dmin
                                dmin = D;
                                dis = d;
                                CTe = C((j-1)*W+1:end);
                            end
                        end
                    end
                    d = dis;
                else
                    d = pdist2(P(Cen(C),:),P(Cen(Tre),:));
                end
                if NC == 1 && length(Tre) == 1
                    dt = d;
                    IT = 1;
                    JT = 1;
                elseif NC == 1
                    [dt,JT] = min(d);
                    IT = 1;
                elseif length(Tre) == 1
                    [dt,IT] = min(d);
                    JT = 1;
                else
                    [d,IT] = min(d);
                    [dt,JT] = min(d);
                    IT = IT(JT);
                end
            else
                dt = 8;
            end
            if ~isempty(F)
                if NC > W
                    m = ceil(NC/W);
                    dmin = 100;
                    for j = 1:m
                        if j < m
                            d = pdist2(P(Cen(C((j-1)*W+1:j*W)),:),P(Cen(F),:));
                            D = min(min(d));
                            if D < dmin
                                dmin = D;
                                dis = d;
                            end
                        else
                            d = pdist2(P(Cen(C((j-1)*W+1:end)),:),P(Cen(F),:));
                            D = min(min(d));
                            if D < dmin
                                dmin = D;
                                dis = d;
                            end
                        end
                    end
                    d = dis;
                else
                    d = pdist2(P(Cen(C),:),P(Cen(F),:));
                end
                df = min(d);
                if length(df) > 1
                    df = min(df);
                end
            else
                df = 10;
            end
            
            d = min([du do dt]);
            
            % Determine what to do with the component
            if NC < Cmin && d > 0.3
                % Remove small isolated component
                Forb(C) = true;
                Other(C) = false;
                CL(C) = nc+3;
            elseif df <= do && df <= dt && df <= du
                % Join the component to "Forb"
                Forb(C) = true;
                Other(C) = false;
                CL(C) = nc+3;
                
            elseif df == 10 && dt == 8 && do == 9 && du == 7
                % Isolated component, do nothing
            else
                if (du <= do && du < dt) || (du <= 2*dmin) 
                    % Join to "Trunk"
                    if NC > W
                        C = CTu;
                    end
                    I = C(IU);
                    J = Tru(JU);
                    if Tree(J)
                        C = Comps{i};
                        Other(C) = false;
                        Tree(C) = true;
                        CL(C) = nc+2;
                        if du > dmax
                            dmax = du;
                        end
                    end
                elseif (dt <= do) || (dt <= 2*dmin)
                    % Join to "Tree"
                    if NC > W
                        C = CTe;
                    end
                    I = C(IT);
                    J = Tre(JT);
                    C = Comps{i};
                    Other(C) = false;
                    Tree(C) = true;
                    CL(C) = nc+2;
                    if dt > dmax
                        dmax = dt;
                    end
                else
                    % Join to other components
                    I = C(IO);
                    J = O(JO);
                    if do > dmax
                        dmax = do;
                    end
                end
                Nei{I} = [Nei{I}; J];
                Nei{J} = [Nei{J}; I];
                Bale{I} = [I; Nei{I}];
                Bale{J} = [J; Nei{J}];
            end
        end
        
        % If "Other" has decreased, do another check with same "distance"
        if nnz(Other) < npre
            again = true;
            npre = nnz(Other);
        else
            again = false;
        end
    end
    k = k+k0;
    Cmin = 3*Cmin;
end
Forb(Base) = false;
%plot_segments(P,Bal,1,Ind(Forb),Base,Ind(Trunk),Ind)