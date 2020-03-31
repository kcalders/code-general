function [Segs, SPar, SChi, CiS, Sta, Axe, Rad, Len, CPar, CExt, CChi, SoC, SDir, Added, ChangeOfSegs]...
    = fill_gaps(Segs, SPar, SChi, CiS, Sta, Axe, Rad, Len, CPar, CExt, CChi, SoC, SDir, dmin, ChangeOfSegs)

% ---------------------------------------------------------------------
% FILL_GAPS.M    	Adds cylinders to fill gaps between cylinders
%
% Version 1.0
% Author        Pasi Raumonen
% Created       14 June 2013
%
% This software cannot be distributed without the permission of P. Raumonen
% and it is only for non-commercial use. These restrictions holds also for
% the modules or subprograms of the software:
% PARTITION_OF_POINT_CLOUD.M, DISTANCES_TO_LINE.M, NEARBY_POINTS
% ---------------------------------------------------------------------


%%
% Interpolates between cylinders or fills gaps between cylinders.
% Finds cylinders without recorded extension or parent cylinder and then
% tries to find one by checking the nearby possible candidates.
% If possible extension or parent cylinder is found, using only the
% existing cylinder data, a cylinder is fitted to this gap.
%
% Inputs:
% Segs      Cover set indexes of the tree segments, cell array
% SPar      Parent segments, vector
% SChi      Child segments, cell array
% CiS       Indexes of the cylinders in segments, cell array
% Sta       Starting points of the cylinders, matrix
% Axe       Axes of the cylinders, matrix
% Rad       Radii of the cylinders, vector
% Len       Lengths of the cylinders, vector
% CPar      Parent cylinders, vector
% CExt      Extension cylinders, vector
% CChi      Child cylinders, cell array
% SoC       Segment index of each sylinder, vector
% dmin      Minimum diameter of the cover sets
%
% Outputs are the same as the inputs but updated by the gap filling process
% Added     Logical vector of the cylinders added to fill gaps


%% Some parameters for finding possible nearby cylinders
%h = 0.5*max(Len);     % parameter defining the size of the nearby region
edgelength = max(Len)/10;  % the edge length of the cubical partition of the starting points
angle = 60;             % the maximum angle (in degrees) for extension
dmax = 0.3;             % the maximum distance between existing cylinders
cosi = cos(angle/180*pi);
%ChangeOfSegs = false;

%% Partition of starting points
% Partition the starting points of the existing cylinders in order to
% quickly find the nearby cylinders
[partition,~,~,info] = partition_of_point_cloud(Sta,edgelength,1);


%% Make the outputs longer
n = length(Rad);
n0 = n;
a = 2*n;
Sta(n+1:a,:) = zeros(n,3);
Axe(n+1:a,:) = zeros(n,3);
Rad(n+1:a) = zeros(n,1);
Len(n+1:a) = zeros(n,1);
CPar(n+1:a) = zeros(n,1);
CExt(n+1:a) = zeros(n,1);
CChi(n+1:a) = cell(n,1);
SoC(n+1:a) = zeros(n,1);


if ~ChangeOfSegs
    %% Check the trunk
    C = CiS{1};
    I = Sta(:,3) ~= 0;
    H = max(Sta(:,3))-min(Sta(I,3));
    h = max(Sta(C,3))-min(Sta(C,3));
    cont = true;
    while cont && h < 0.9*H
        C = CiS{1};
        c = vertcat(CChi{C(end)});
        I = Axe(c,:)*Axe(C(end),:)';
        [cosine,J] = max(I);
        if cosine > 0.8
            c = c(J);
            C = C(end);
            si = SoC(c);
            ChangeOfSegs = true;
            Segs{1} = [Segs{1}; Segs{si}];
            Segs{si} = cell(0,1);
            I = SPar == si;
            SPar(I) = 1;
            SChi{1} = setdiff(SChi{1},si);
            SChi{1} = [SChi{1}; SChi{si}];
            CChi{C} = setdiff(CChi{C},c);
            CExt(C) = c;
            I = CiS{si};
            SoC(I) = 1;
            CiS{1} = [CiS{1}; I];
            CiS{si} = zeros(0,1);
            
            C = CiS{1};
            h = max(Sta(C,3))-min(Sta(C,3));
        else
            cont = false;
        end
        
    end
    %S = Segs{1};
    %plot_segments(P,Bal,2,vertcat(S{:}))
    %pause
end

t = 0;
%% Then check cylinders and their extensions
for i = 1:n
    if CExt(i) == 0 % No extension exists, check if a suitable can be defined
        % The starting point of the possible new cylinder
        S = Sta(i,:)+Len(i)*Axe(i,:);
        C = CChi{i};
        if isempty(C)
            % No child cylinders, search suitable nearby cylinders
            Search = true;
        else
            % Child cylinders exists, check them first
            Cosines = Axe(C,:)*Axe(i,:)';
            I = Cosines > 0.8;
            if any(I)
                C = C(I);
                A = mat_vec_subtraction(Sta(C,:),Sta(i,:));
                h = A*Axe(i,:)';
                E = Sta(C,:)+[Len(C).*Axe(C,1) Len(C).*Axe(C,2) Len(C).*Axe(C,3)];
                A = mat_vec_subtraction(E,Sta(i,:));
                he = A*Axe(i,:)';
                I = h > 0.75*Len(i);
                J = he > Len(i);
                I = I&J;
                if nnz(I) == 1
                    C = C(I);
                    E = E(I,:);
                    V = E-S;
                    Len(C) = norm(V);
                    Axe(C,:) = V/Len(C);
                    Sta(C,:) = S;
                    CExt(i) = C;
                    CChi{i} = setdiff(CChi{i},C);
                    p = SoC(i);
                    e = SoC(C);
                    SoC(C) = p;
                    SoC(CiS{e}) = p;
                    Segs{p} = [Segs{p}; Segs{e}];
                    Segs{e} = cell(0);
                    SPar(SChi{e}) = p;
                    S = SChi{p};
                    J = S ~= e;
                    SChi{p} = S(J);
                    SChi{p} = [SChi{p}; SChi{e}];
                    SChi{e} = zeros(0,1);
                    CiS{p} = [CiS{p}; CiS{e}];
                    CiS{e} = zeros(0,1);
                    Search = false;
                    ChangeOfSegs = true;
                    
                elseif nnz(I) > 1
                    Search = false;
                else
                    Search = true;
                end
            else
                Search = true;
            end
                
        end
        
        if Search
            % Search the nearby startingpoints, given by the index vector N
            P = S+4*Len(i)*Axe(i,:); % P and S define the neighborhood used to search nearby cylinders
            Q = [S; P];
            px = floor((Q(:,1)-info(1))/edgelength)+2;
            py = floor((Q(:,2)-info(2))/edgelength)+2;
            pz = floor((Q(:,3)-info(3))/edgelength)+2;
            N = nearby_points(partition,px,py,pz,info,5);
            N = setdiff(N,i);
            I  = (CPar(N) == 0); % only those nearby cylinders without parent
            N = N(I);
            if ~isempty(N)  % nearby cylinders without parent exists
                D = mat_vec_subtraction(Sta(N,:),S); % displacement vectors
                d = sqrt(sum(D.*D,2));  % norms of the disp. vectors
                cosines = 1./d.*(D*Axe(i,:)');  % angles between cylinder and disp.vec.
                I = (cosines > cosi);  % only those disp.vec. with small angle
                J = d < min([5*Rad(i) dmax]);
                I = I&J;
                if nnz(I) > 0  % suitable cylinders for extensions exists
                    % Search the extension cylinder with the smallest angle
                    % with its new possible parent
                    N = N(I);
                    cosines = 1./d(I).*(sum(D(I,:).*Axe(N,:),2));
                    [m,ind] = max(cosines);
                    if m > cosi  % only add a cylinder with small angle with its extension
                        % Define the index of the selected disp. vector
                        K = 1:1:length(I);
                        N = N(ind);
                        K = K(I);
                        ind = K(ind);
                        
                        % add new cylinder and update relevant variables
                        t = t+1;
                        Sta(n+t,:) = S;
                        Axe(n+t,:) = D(ind,:)/norm(D(ind,:));
                        Rad(n+t) = 0.5*(Rad(i)+Rad(N));
                        Len(n+t) = norm(D(ind,:));
                        CExt(i) = n+t;
                        CExt(n+t) = N;
                        CPar(N) = n+t;
                        CPar(n+t) = i;
                        p = SoC(i);
                        e = SoC(N);
                        SoC(n+t) = p;
                        SoC(CiS{e}) = p;
                        Segs{p} = [Segs{p}; Segs{e}];
                        Segs{e} = cell(0);
                        SPar(SChi{e}) = p;
                        Se = SChi{p};
                        J = Se ~= e;
                        SChi{p} = Se(J);
                        SChi{p} = [SChi{p}; SChi{e}];
                        SChi{e} = zeros(0,1);
                        CiS{p} = [CiS{p}; n+t; CiS{e}];
                        CiS{e} = zeros(0,1);
                        ChangeOfSegs = true;
                        
                        % add the starting point of the new cylinder in the partition
                        px = floor((S(1)-info(1))/edgelength)+2;
                        py = floor((S(2)-info(2))/edgelength)+2;
                        pz = floor((S(3)-info(3))/edgelength)+2;
                        partition{px,py,pz} = [partition{px,py,pz}; n+t];
                        
                    end
                end
            end
        end
    end
end


%% Check cylinders and their parentcyls
% Add cylinders between a cylinder and its parent cylinder if there is a
% large enough gap between them, the cylinder is not the extension of 
% its parent
D = min(dmin,0.001);  % if the gap is longer than D, then the cylinder is added
for i = 2:n
    if CPar(i) > 0 && CExt(CPar(i)) ~= i % Parent exists
        PC = CPar(i);  % parent cylinder
        [d,V,h,B] = distances_to_line(Sta(i,:),Axe(PC,:),Sta(PC,:)); 
        if d > Rad(PC)
            % Determine the parameters of the cylinder to be (possibly) added 
            L = norm(V);
            A = V/L;
            if h < 0 % Cyl "i" under its parent
                S = Sta(PC,:)+Rad(PC)*A;
                V = Sta(i,:)-S;
                L = norm(V);
                A = V/L;
            elseif h > Len(PC) % Cyl "i" above its parent
                S = Sta(PC,:)+Len(PC)*Axe(PC,:)+Rad(PC)*A;
                V = Sta(i,:)-S;
                L = norm(V);
                A = V/L;
            else
                L = L-Rad(PC);
                S = Sta(PC,:)+B+Rad(PC)*A;
            end
            R = Rad(i);
            
            % Add cylinder if the gap is long enough, cylinder is elongated
            % enough, or its lenght (volume) is large enough compared to
            % its extension
            if L/R > 1 || L > D || L/Len(i) > 0.333
                t = t+1;
                Rad(n+t) = R;
                Len(n+t) = L;
                Axe(n+t,:) = A;
                Sta(n+t,:) = S;
                CExt(n+t) = i;
                CPar(i) = n+t;
                CPar(n+t) = PC;
                s = SoC(i);
                SoC(n+t) = s;
                CiS{s} = [n+t; CiS{s}];
                % add the starting point of the new cylinder in the partition
                px = floor((S(1)-info(1))/edgelength)+2;
                py = floor((S(2)-info(2))/edgelength)+2;
                pz = floor((S(3)-info(3))/edgelength)+2;
                partition{px,py,pz} = [partition{px,py,pz}; n+t];
                
            end
        end
    end
end
n = n+t;


%% Display results
str = ['    ',num2str(n-n0),' cylinders added'];
disp(str)
str = ['    Total of ',num2str(n),' cylinders in the model'];
disp(str)


%% Finalize the outputs
Sta = Sta(1:n,:);
Axe = Axe(1:n,:);
Rad = Rad(1:n);
Len = Len(1:n);
CExt = CExt(1:n);
CPar = CPar(1:n);
SoC = SoC(1:n);
CChi = CChi(1:n);
Added = false(n,1);
Added(n0+1:1:n) = true;


%% Plot results
% Plot the trunk cylinders
% C = vertcat(CiS{1});
% plot_cylinder_model(Rad(C),Len(C),Axe(C,:),Sta(C,:),1,10,1)

% Plot the added cylinders
%if n > n0
%    plot_cylinder_model(Rad(n0+1:end),Len(n0+1:end),Axe(n0+1:end,:),Sta(n0+1:end,:),2,10,1)
%end

% Plot all the cylinders
% plot_cylinder_model(Rad(1:n0),Len(1:n0),Axe(1:n0,:),Sta(1:n0,:),3,6,1)
% if n > n0
%     hold on
%     plot_cylinder_model(Rad(n0+1:end),Len(n0+1:end),Axe(n0+1:end,:),Sta(n0+1:end,:),2,36,1)
%     hold off
% end
% pause

end % End of main function


function NearbyPoints = nearby_points(partition,px,py,pz,info,n)

% For given cover sets "sets", defines the nearby cover sets "nbsets"
% which are in a cube containing the "sets" and extending "n"
% cubes in each direction from the cubes of the "sets".

s = info(4:6);
p = sort(px);
minx = p(1)-n;
if minx < 1
    minx = 1;
end
maxx = p(2)+n;
if maxx > s(1)
    maxx = s(1);
end
p = sort(py);
miny = p(1)-n;
if miny < 1
    miny = 1;
end
maxy = p(2)+n;
if maxy > s(2)
    maxy = s(2);
end
p = sort(pz);
minz = p(1)-n;
if minz < 1
    minz = 1;
end
maxz = p(2)+n;
if maxz > s(3)
    maxz = s(3);
end
NearbyPoints = vertcat(partition{minx:maxx,miny:maxy,minz:maxz});
end
