function [Rad,Len,Axe,Sta,CPar,CExt,CChi,SoC,CiS,Segs,SPar,SChi,ChaOfSegs] = ...
    cylinders(P,Bal,Cen,Segs,SPar,SDir,SChi,dmin,lcyl)

% ---------------------------------------------------------------------
% CYLINDERS.M       Fits cylinders to the segments
%
% Version 1.2
% Author        Pasi Raumonen
% Created       2 December 2013
%
% This software cannot be distributed without the permission of P. Raumonen
% and it is only for non-commercial use. These restrictions holds also for
% the modules or subprograms of the software:
% FIX_SEGMENTS, REGIONS, COVER_EXTENT, REFINEMENT_FITTING, PARENT_CYLINDER,
% ADJUSTMENTS, 
% ---------------------------------------------------------------------

% Analyses segments of a tree/roots by approximating them with cylinders.
% Subdivides each segment to smaller regions for which cylinders are
% fitted. Returns the cylinder information and the
% child/parent/extension-relation of the cylinders.
%
% Inputs:
% P         Point cloud, matrix
% Bal       Cover sets (balls), cell array
% Cen       Center points of the cover sets, vector
% Segs      Segments of the tree, cell array of cell arrays
% SPar      Parent segments (indexes), vector
% SDir      Directions of the segments at their bases, matrix
% SChi      Child segments (cell array)
% dmin      Minimum diameter of the cover sets
% lcyl      Cylinder length/radius ratio
%
% Outputs:
% Rad       Radii of the cylinders, vector
% Len       Lengths of the cylinders, vector
% Axe       Axes of the cylinders, matrix
% Sta       Starting points of the cylinders, matrix
% CPar      Parents of the cylinders, vector
% CExt      Extensions of the cylinders, vector
% CChi      Children of the cylinders (extension not included), cell array
% SoC       Segments the cylinders belong, vector
% CiS       Cylinders in the segments, cell array
% Segs, SPar, SChi      Updated segment structures
% ChaOfSegs     Logical truth value weather segments have changed


%% Initializations
% Initialize the outputs
NumOfSeg = max(size(Segs));   % number of segments
Rad = zeros(40*NumOfSeg,1);
Len = zeros(40*NumOfSeg,1);
Axe = zeros(40*NumOfSeg,3);
Sta = zeros(40*NumOfSeg,3);
CPar = zeros(40*NumOfSeg,1);
CExt = zeros(40*NumOfSeg,1);
CChi = cell(40*NumOfSeg,1);
CiS = cell(NumOfSeg,1);
SoC = zeros(40*NumOfSeg,1);
ChaOfSegs = false;

% Trunk position and direction
S = Segs{1};
S = vertcat(S{:});
X = cov(P(Cen(S),:));
[U,~,~] = svd(X);
Dt = U(:,1);
Mt = mean(P(Cen(S),:));
Bot = min(P(Cen(S),3));
Top = max(P(Cen,3));
Height = Top-Bot;
d = mat_vec_subtraction(P(Cen,1:2),Mt(1:2));
Width = max(sqrt(sum(d.*d,2)));

c = 1;  % number of cylinders determined

% Determine suitable order of segments (from trunk to "youngest" child)
S = 1;
t = 1;
Ind = zeros(NumOfSeg,1);
Ind(1) = 1;
while ~isempty(S)
    S = vertcat(SChi{S});
    n = length(S);
    Ind(t+1:t+n) = S;
    t = t+n;
end

% Fit cylinders individually for each segment
for k = 1:NumOfSeg
    si = Ind(k);
    if si > 0 && ~isempty(Segs{si})
        %% Some initialization about the segment
        Seg = Segs{si};      % the current segment under analysis
        Segment = vertcat(Seg{:});
        direc = SDir(si,:);      % direction of the current segment at its base
        nl = max(size(Seg));    % number of layers in the segment
        Base = Seg{1};          % the base of the segment
        if si == 1
            Trunk = true;
        else
            Trunk = false;
        end
        
        % Check if the base is at the end of the segment
        if nl > 4
            [Seg,Segment,Base,nl] = fix_segment(P,Cen,Seg,Segment,Base,nl);
        end
        
        nc = length(Segment);   % number of cover sets in the current segment
        Q = vertcat(Bal{Segment}); % the points in the segments
        np = length(Q);         % number of points in the segment
        I = vertcat(Bal{Base}); % the points of the base
        nb = length(I);         % number of points in the base

        if np > nb && (nc > 2) && (np > 20) && (~isempty(Base)) % analyse only large enough segments
            
            %% Determine the regions for cylinder fitting
            [Regs,Starts,Lengs,Axes,Rads,Mads,t,NL] = regions(P,Bal,Seg,direc,nl,Trunk,Cen,dmin,lcyl);
            
            %% Fit cylinders to the regions
            FQ = zeros(t,4);   % the fitting quality
            for j = 1:t
                if (length(Regs{j}) > 10) && (norm(Axes(j,:)) > 0) % fit cylinders to large enough subsegs
                    
                    Points = P(Regs{j},:);  % the coordinate points used for fitting
                    
                    % Initial estimates
                    CA0 = Axes(j,:)';     % cylinder axis
                    AP0 = Starts(j,:)'; % mean(Points)';   % point in the cylinder axis
                    R0 = Rads(j,1);
                    mad0 = Mads(j);
                    FQ(j,1) = mad0;
                    
                    % First fitting
                    [AP,CA,R,d,~,conv,rel] = lscylinder(Points,AP0,CA0,R0,0.1,0.1);
                    
                    % Conditions for second fitting and accepting the
                    % results
                    I1 = conv == 1 & rel == 1; % fitting converged and is reliable
                    I2 = ~(isnan(R)|any(isnan(AP))|any(isnan(CA))); % results are numbers
                    mad = mean(abs(d));  % mean distance to the fitted cylinder
                    md = max(d);  % maximum distance
                    I3 = mad < R0 & abs(CA0'*CA) > 0.8; % distances and the angle acceptable
                    I4 = R < 3*R0 & R > 0; % radius is acceptable
                    % second fitting if large enough absolute and relative "errors"
                    Second = mad > 0.02 & mad/R > 0.1 & md/R > 0.333 & R > 1.5*R0;
                    Accept = I1&I2&I3&I4; % accept the first fitting
                    Second = Accept & Second; % second cylinder fitting
                    
                    % Possible second fitting
                    if Second
                        
                        I = (d < 0.5*md);
                        if (nnz(I) < 10) || (nnz(I) < 0.25*length(I))
                            I = d < 0.25*md;
                            if (nnz(I) < 10) || (nnz(I) < 0.25*length(I))
                                I = true(length(d),1);
                            end
                        end
                        Points = Points(I,:);
                        
                        [AP,CA,R,d,~,conv,rel] = lscylinder(Points,AP,CA,R,0.1,0.1);
                        
                        I1 = conv == 1 & rel == 1; % fitting converged and is reliable
                        I2 = ~(isnan(R)|any(isnan(AP))|any(isnan(CA))); % results are numbers
                        mad = mean(abs(d));  % mean distance to the fitted cylinder
                        I3 = mad < R0 & abs(CA0'*CA) > 0.8; % distances and the angle acceptable
                        I4 = R < 3*R0 & R > 0; % radius is acceptable
                        Accept = I1&I2&I3&I4; % accept the first fitting
                        
                        H = Points*CA;
                        hmin = min(H);
                        L = abs(max(H)-hmin);
                        hpoint = CA'*AP;
                        SP = AP'-(hpoint-hmin)*CA';
                        CovExt = cover_extent(Points,CA,SP,R,L);
                    elseif Accept
                        H = Points*CA;
                        hmin = min(H);
                        L = abs(max(H)-hmin);
                        hpoint = CA'*AP;
                        SP = AP'-(hpoint-hmin)*CA';
                        CovExt = cover_extent(Points,CA,SP,R,L);
                    else
                        CovExt = cover_extent(Points,CA0,AP0,R0,Lengs(j));
                    end
                    
                    if Accept
                        Rads(j) = R;
                        Axes(j,:) = CA';
                        Lengs(j) = L;
                        Starts(j,:) = SP;
                        FQ(j,2:4) = [mad CovExt];
                    else
                        % do not accept least square fittings, use initial
                        % estimates
                        FQ(j,2:4) = [mad0 CovExt];
                    end
                end
            end
            
            %% "Refinement" fitting for "large" segments
            % Each region is devided into three parts where the
            % cylinders are fitted
            [Rads,Lengs,Starts,Axes,FQ,t] = ...
                refinement_fitting(P,Regs,Rads,Lengs,Starts,Axes,FQ,t,lcyl);
            
            
            %% Search possible parent cylinder
            if t > 0
                [PC,Starts,Lengs,Axes,Segs,SPar,SChi,CExt,t,ChaOfSegs] = ...
                    parent_cylinder(P,Regs,Segs,SPar,SChi,CiS,CExt,...
                    Rad,Axe,Sta,Len,Starts,Axes,Rads,Lengs,c,si,t,ChaOfSegs);
            end
            
            
            %% Adjust cylinders
            if t > 0
                [Rads,Starts,Axes,Lengs,t] = adjustments(P,Rad,Sta,CiS,Rads,...
                    Starts,Axes,Lengs,Regs,FQ,Width,Height,Dt,Mt,PC,si,t,SoC,SChi,SPar);
            end
             
            
            %% Adjust the radii for very small branches
            if max(Rads) < 0.005
                if t > 2
                    r = sort(Rads);
                    r = r(2:end-1);
                    a = 2*mean(r);
                    if a > max(r)
                        a = min(0.01,max(r));
                    end
                    b = min(0.5*min(Rads),0.001);
                    Rads = linspace(a,b,t)';
                else
                    if max(Rads) > 0.01
                        r = 0.01;
                    else
                        r = max(Rads);
                    end
                    if t == 1
                        Rads = r;
                    else
                        Rads = [r 0.5*r];
                    end
                end
            end
           
            %% Save the cylinders
            % if at least one acceptable cylinder, then save them
            I = sum(Axes.*Axes,2);
            J = sum(Starts.*Starts,2);
            if (t > 0) && (min(Rads(1:t)) > 0) && ~any(I == 0) && ~any(J == 0)
                % If the parent cylinder exists, set the parent-child relations
                if ~isempty(PC)
                    CPar(c) = PC;
                    if CExt(PC) == c
                        I = SoC(PC);
                        SoC(c:c+t-1) = I;
                        CiS{I} = [CiS{I}; linspace(c,c+t-1,t)'];
                    else
                        CChi{PC} = [CChi{PC}; c];
                        SoC(c:c+t-1) = si;
                        CiS{si} = linspace(c,c+t-1,t)';
                    end
                else
                    SoC(c:c+t-1) = si;
                    CiS{si} = linspace(c,c+t-1,t)';
                end
                                
                Rad(c:c+t-1) = Rads(1:t);
                Len(c:c+t-1) = Lengs(1:t);
                Axe(c:c+t-1,:) = Axes(1:t,:);
                Sta(c:c+t-1,:) = Starts(1:t,:);
                CPar(c+1:c+t-1) = linspace(c,c+t-2,t-1);
                CExt(c:c+t-2) = linspace(c+1,c+t-1,t-1);
                
                c = c+t; % number of cylinders so far (plus one)
                
            end
        end
    end
end
c = c-1; % number of cylinders 

str = ['    ',num2str(c),' cylinders fitted'];
disp(str)

%% Finalize outputs
Rad = Rad(1:c);
Len = Len(1:c);
Axe = Axe(1:c,:);
Sta = Sta(1:c,:);
SoC = SoC(1:c);
CPar = CPar(1:c);
CExt = CExt(1:c);
CChi = CChi(1:c,:);
for si = 1:NumOfSeg
    if size(CiS{si},2) > 1
        CiS{si} = CiS{si}';
    end
end
for si = 1:c
    if size(CChi{si},2) > 1
        CChi{si} = CChi{si}';
    end
end

%plot_cylinder_model(Rad,Len,Axe,Sta,3,20,0.3)


end % End of main function


function [Seg,Segment,Base,nl] = fix_segment(P,Cen,Seg,Segment,Base,nl)

% Project the points into the principal component and check if the base
% is at the other end. If not and segment is short, remove the segment
% from the cylinder fitting process. For larger such segments, try to
% fix the problem by first expanding the base. If that does not work.
% then remove the other end of the segment.

if nl > 4
    S = vertcat(Seg{1:5});
else
    S = vertcat(Seg{1:4});
end
X = cov(P(Cen(S),:));
[U,~,~] = svd(X);   % principal components
dps = P(Cen(S),:)*U(:,1); % project the segment
ds1 = min(dps);     % the minimum projection value of the segment
ds2 = max(dps);     % the maximum projection value of the segment
dpb = P(Cen(Base),:)*U(:,1); % project the base
db1 = min(dpb);     % the minimum projection value of the base
db2 = max(dpb);     % the maximum projection value of the base
if (ds1 ~= db1) && (ds2 ~= db2)
    % The base is not at the other, try to fix it by expanding the base
    Base = vertcat(Seg{1:4});   
    seg = cell(nl-3,1);         % update the Seg using the cell "seg"
    seg{1} = Base;
    seg(2:end,1) = Seg(5:end,1);

    % The updated inputs
    Seg = seg;
    Segment = vertcat(Seg{:});
    nl = nl-3;
end
end % End of function


function [Regs,Starts,Lengs,Axes,Rads,Mads,t,NL] = regions...
    (P,Bal,Seg,direc,nl,Trunk,Cen,dmin,lcyl)

% Define the regions for cylinder fitting

if nl > 3
    % Determine number of set layers (nlayer) for initial regions
    if Trunk
        i = 3;   i0 = 3;   H = 6;
    else
        i = 3;   i0 = 1;   H = lcyl;
    end
    while (i < nl) && (length(vertcat(Bal{vertcat(Seg{1:i})})) < 5)
        i = i+1;
    end
    if i+1 > nl
        i = nl-1;
    end
    Test = vertcat(Seg{i0:i+1});
    Test = vertcat(Bal{Test});
    B = vertcat(Bal{vertcat(Seg{i0:i})});
    if length(B) > 1
        B = mean(P(B,:));
    else
        B = P(B,:);
    end
    T = vertcat(Bal{Seg{i+1}});
    T = mean(P(T,:));
    V = T-B;
    [d,~,h] = distances_to_line(P(Test,:),V,B);
    R = median(d);
    if R < dmin
        H = max(16,ceil(dmin/R*H));
    end
    L = max(h)-min(h);
    while i < nl-1 && L < H*R
        i = i+1;
        T = vertcat(Bal{Seg{i+1}});
        T = mean(P(T,:));
        V = T-B;
        Test = vertcat(Seg{i0:i+1});
        Test = vertcat(Bal{Test});
        [d,~,h] = distances_to_line(P(Test,:),V,B);
        R = median(d);
        L = max(h)-min(h);
    end
    NL = min(20,i-i0);
    NL = max(3,NL);
    if NL < 5 && L/R > 10
        NL = 2;
    elseif NL >= 6 && ~Trunk && L/R > 12
        NL = max(3,ceil(NL/2));
    end
        
    if Trunk
        i = 1;  % include the base for the trunk
    else
        i = 1; % do not include the base for brances
    end
    t = 0;
    Points = zeros(nl,3);
    RegsExt = cell(nl,1);
    while i <= nl
        j = i+NL-1;
        if j > nl
            j = nl;
        end
        S = vertcat(Seg{i:j});
        if length(S) > 50
            t = t+1;
            RegsExt{t} = S;
            Points(t,:) = mean(P(Cen(S),:));
        else
            p = vertcat(Bal{S});
            if ~isempty(p)
                t = t+1;
                RegsExt{t} = S;
                if length(p) > 1
                    Points(t,:) = mean(P(p,:));
                elseif length(p) == 1
                    Points(t,:) = P(p,:);
                end
            end
        end
        i = i+NL;
    end
    Points = Points(1:t,:);
    RegsExt = RegsExt(1:t);
    
    if t > 1
        t = t-1;
        Regs = cell(t,1);    % the points (indexes) of the final subsegs.
        Lengs = zeros(t,1);     % Length of the subsegs
        Starts = zeros(t,3);   % starting points of the subsegs
        Axes = zeros(t,3);
        Rads = zeros(t,1);
        Mads = zeros(t,1);  % mean absolute distances to the cyl surface
        for j = 1:t
            if j == 1 % the first subseg
                B = vertcat(RegsExt{1:2});
                if length(B) < 60
                    I = vertcat(Bal{B});
                else
                    I = Cen(B);
                end
                V = Points(j+1,:)-Points(j,:);
                Axes(j,:) = V/norm(V);
                [d,~,h] = distances_to_line(P(I,:),V,Points(j,:));
                J = h <= norm(V);
                d = d(J);
                h = h(J);
                I = I(J);
                K = d <= 3*mean(d);
                Rads(j) = mean(d(K));
                Mads(j) = mean(abs(d(K)-Rads(j)));
                I = I(K);
                h = h(K);
                Regs{j} = I;
                if nnz(K) > 0
                    minH =  min(h);
                    Lengs(j) = max(h)-minH;
                    H = P(I,:)*Axes(j,:)';
                    hpoint = Points(j,:)*Axes(j,:)';
                    Starts(j,:) = Points(j,:)-(hpoint-min(H))*Axes(j,:);
                else
                    Lengs(j) = 0.01;
                    Starts(j,:) = Points(j,:);
                end
            elseif j == t  % the last subseg
                B = vertcat(RegsExt{j:j+1});
                if length(B) < 60
                    I = vertcat(Bal{B});
                else
                    I = Cen(B);
                end
                V = Points(j+1,:)-Points(j,:);
                Axes(j,:) = V/norm(V);
                [d,~,h] = distances_to_line(P(I,:),V,Points(j,:));
                J = h >= 0;
                d = d(J);
                h = h(J);
                I = I(J);
                J = d <= 3*mean(d);
                Rads(j) = mean(d(J));
                Mads(j) = mean(abs(d(J)-Rads(j)));
                I = I(J);
                h = h(J);
                if any(J) % not empty
                    Regs{j} = I;
                    minH = min(h);
                    Lengs(j) = max(h)-minH;
                    Starts(j,:) = Points(j,:);
                else
                    Starts(j,:) = Points(j,:);
                    Lengs(j) = Lengs(j-1);
                end
            else % the middle subsegs
                B = vertcat(RegsExt{j:j+1});
                if length(B) < 60
                    I = vertcat(Bal{B});
                else
                    I = Cen(B);
                end
                V = Points(j+1,:)-Points(j,:);
                L = norm(V);
                Axes(j,:) = V/L;
                [d,~,h] = distances_to_line(P(I,:),V,Points(j,:));
                J = h >= 0;
                K = h <= L;
                J = J&K;
                d = d(J);
                I = I(J);
                J = d <= 3*mean(d);
                Rads(j) = mean(d(J));
                Mads(j) = mean(abs(d(J)-Rads(j)));
                Regs{j} = I(J);
                Lengs(j) = L;
                Starts(j,:) = Points(j,:);
            end
        end
    else
        Regs = cell(1,1);
        Regs{1} = vertcat(Bal{RegsExt{:}});
        [d,~,h] = distances_to_line(P(Regs{1},:),direc,Points);
        Lengs = max(h)-min(h);
        Rads = median(d);
        Mads = mean(abs(d-Rads));
        Axes = direc;
        H = P(Regs{:},:)*Axes';
        hpoint = Points*Axes';
        Starts = Points-(hpoint-min(H))*Axes;
        t = 1;
    end
elseif nl > 1
    Regs = cell(1,1);
    Regs{1} = vertcat(Bal{vertcat(Seg{2:end})});
    Starts = mean(P(Regs{:},:));
    if max(size(Starts)) == 3
        [d,~,h] = distances_to_line(P(Regs{1},:),direc,Starts);
        Lengs = max(h)-min(h);
        Rads = median(d);
        Mads = mean(abs(d-Rads));
        Axes = direc;
        H = P(Regs{:},:)*Axes';
        hpoint = Starts*Axes';
        Starts = Starts-(hpoint-min(H))*Axes;
        t = 1;
        NL = 1;
    else
        t = 0;
        Axes = 0;
        Rads = 0;
        Lengs = 0;
        Mads = 0;
        NL = 0;
    end
end

if (t > 1) && (length(Regs{t}) < 11)
    Regs{t-1} = [Regs{t-1}; Regs{t}];
    t = t-1;
    Regs = Regs(1:t);
    Rads = Rads(1:t);
    Lengs = Lengs(1:t);
    Axes = Axes(1:t,:);
    Starts = Starts(1:t,:);
end
end % End of function


function CovExt = cover_extent(Points,CylAxis,AxisPoint,R,L)

% Calculates CovExt that measures how broadly/well the arc of circle 
% (the surface of cylinder) is covered. Divides the points into three parts
% and each of the parts into 16 sectors. Gives the mean of sectors covered
% and an indication if at least one of the parts has opposite sectors
% covered.

% Calculate the distances of the points to the axis
[d,V,h] = distances_to_line(Points, CylAxis, AxisPoint);

% Select points close to the cylinder surface and divide the points into
% three parts in the direction of the axis
I = d > 0.9*R;
J = d < 1.1*R;
I = I&J;
V = V(I,:);
h = h(I);
L1 = h < 0.333*L;  % first part
L2 = h < 0.666*L;  
L3 = ~L2;          % third part
L2 = ~L1 & ~L3;    % second part

% Define vectors orthonormal to the axis to be used as a reference for
% angle determination
[U,W] = orthonormal_vectors(CylAxis);

A = V*[U W];
l = sqrt(sum(A.*A,2));
A = [A(:,1)./l A(:,2)./l];

angles = atan2(A(:,1),A(:,2))+pi; % the angles of the points
sectors = ceil(angles*8/pi);  % the sectors (16 of them) of the points determined by the angles

I1 = false(16,1);   % sectors for the first part
I2 = I1;            % sectors for the second part
I3 = I1;
I1(sectors(L1)) = true;  % sectors occupied by the points in the first part
I2(sectors(L2)) = true;
I3(sectors(L3)) = true;

% The first element of the CovExt, which is the mean of sectors covered by
% the parts (the relative value)
CovExt = zeros(1,2);
CovExt(1) = mean([sum(I1)/16 sum(I2)/16 sum(I3)/16]);

% The second element of the CovExt, equals to one if in one of the parts is
% such that nearly opposite sectors are covered
J1 = I1(1:8)&I1(9:16);
J1 = any(J1);
K1 = I1&[I1(8:16); I1(1:7)];
K1 = any(K1);
L1 = I1&[I1(10:16); I1(1:9)];
L1 = any(L1);
J2 = I2(1:8)&I2(9:16);
J2 = any(J2);
K2 = I2&[I2(8:16); I2(1:7)];
K2 = any(K2);
L2 = I2&[I2(10:16); I2(1:9)];
L2 = any(L2);
J3 = I3(1:8)&I3(9:16);
J3 = any(J3);
K3 = I3&[I3(8:16); I3(1:7)];
K3 = any(K3);
L3 = I3&[I3(10:16); I3(1:9)];
L3 = any(L3);
CovExt(2) = (J1|K1|L1) | (J2|K2|L2) | (J3|K3|L3);

end % End of function


function [Rads,Lengs,Starts,Axes,FQ,t] = refinement_fitting...
    (P,Regs,Rads,Lengs,Starts,Axes,FQ,t,lcyl)

% Makes a new fitting such that every region is partitioned into three
% parts using the corresponding cylinder and new cylinders are fitted to
% these subregions using the previous cylinders as initial estimates.
               
rads = zeros(30*t,1);
lengs = zeros(30*t,1);
starts = zeros(30*t,3);
axes = zeros(30*t,3);
fq = zeros(30*t,4);  % fit quality for refined cylinders
c = 0;
for i = 1:t
    Points = P(Regs{i},:);  % the coordinate points used for fitting
    
    Q = mat_vec_subtraction(Points,Starts(i,:));
    d = Q*Axes(i,:)';
    
    ratio = Lengs(i)/Rads(i)/lcyl;
    if i > 1 && Rads(i) > 1.25*Rads(i-1)
        ratio = Rads(i)/Rads(i-1)*ratio;
    end
    if ratio > 4
        ratio = 4;
    end
    if ratio >= 1.25
        n = ceil(ratio);
        m = length(d);
        I = false(m,n);
        S = Starts(i,:);  L = Lengs(i);  A = Axes(i,:);  R = Rads(i);
        AP = zeros(n,3);
        rad = zeros(n,1); len = rad; sta = zeros(n,3); axe = sta;
        j = 1;
        Sub = true(n,1); % the subsegments/regions
        while j <= n
            J = d < j*L/n;
            K = d > (j-1)*L/n;
            if nnz(J&K) > 10
                I(:,j) = J&K;
                if j == 1
                    AP(j,:) = S;
                else
                    AP(j,:) = S+(j-1)/n*L*A;
                end
            else
                J = d < (j+1)*L/n;
                I(:,j) = J&K;
                if j == 1
                    AP(j,:) = S;
                else
                    AP(j,:) = S+(j-1)/n*L*A;
                end
                j = j+1;
                Sub(j) = false;
            end
            j = j+1;
        end
        if ~all(Sub)
            n = nnz(Sub);
            AP = AP(Sub,:);
            rad = rad(Sub);
            len = len(Sub);
            sta = sta(Sub,:);
            axe = axe(Sub,:);
        end
        
        d0 = distances_to_line(Points,A,S)-Rads(i);
        
        for j = 1:n
            J = I(:,j);
            
            if nnz(J) > 10
                mad0 = mean(abs(d0(J)));
                c = c+1;
                fq(c,1) = mad0;
                
                [s,a,r,d,~,conv,rel] = lscylinder(Points(J,:),AP(j,:)',A',R,0.1,0.1);
                
                mad = mean(abs(d));
                H = Points(I(:,j),:)*a;
                hmin = min(H);
                l = abs(max(H)-hmin);
                hpoint = s'*a;
                s = s'-(hpoint-hmin)*a';
                CovExt = cover_extent(Points(J,:),a,s,r,l);
                
                fq(c,2:4) = [mad CovExt];
                
                % Check if the fitting can be accepted
                I1 = isnan(r)|isnan(a(1))|isnan(s(1))|(conv == 0)|(rel == 0);
                I2 = r > 1.2*R & mad > 0.75*mad0;
                I3 = r > 1.5*R & r < 0.5*R;
                I4 = A*a < 0.8;
                Reject = I1|I2|I3|I4;
                if Reject
                    % Reject the fitting
                    r = R;
                    a = A';
                    s = AP(j,:)';
                    l = L/n;
                end
            else
                c = c+1;
                r = R;
                a = A';
                s = AP(j,:)';
                l = L/n;
                mad0 = mean(abs(d0(J)));
                CovExt = cover_extent(Points(J,:),a,s,r,l);
                fq(c,:) = [mad0 FQ(i,2) CovExt];
            end
            rad(j) = r;
            len(j) = l;
            sta(j,:) = s';
            axe(j,:) = a';
        end
        I = axe*A' < 0.35;
        if any(I) || sum(len) > 1.2*L
            % reject the fittings
            rads(c-n+1:c) = R;
            lengs(c-n+1:c) = L/n;
            starts(c-n+1:c,:) = AP;
            axes(c-n+1:c,:) = repmat(A,[n 1]);
            fq(c-n+1:c,:) = repmat(FQ(i,:),[n,1]);
        else
            rads(c-n+1:c) = rad;
            lengs(c-n+1:c) = len;
            starts(c-n+1:c,:) = sta;
            axes(c-n+1:c,:) = axe;
        end
    else
        c = c+1;
        fq(c,:) = FQ(i,:);
        rads(c) = Rads(i);
        lengs(c) = Lengs(i);
        starts(c,:) = Starts(i,:);
        axes(c,:) = Axes(i,:);
    end
end
rads = rads(1:c);
lengs = lengs(1:c);
starts = starts(1:c,:);
axes = axes(1:c,:);
fq = fq(1:c,:);

I = true(c,1);
for i = 1:c
    j = 0;
    while i+j+1 <= c  &&  rads(i+j+1) == rads(i)
        j = j+1;
    end
    if j > 0
        I(i+1:i+j) = false;
        lengs(i) = sum(lengs(i:i+j));
    end
end
t = nnz(I);
rads = rads(I);
lengs = lengs(I);
starts = starts(I,:);
axes = axes(I,:);


i = 2;
while i <= t
    if fq(i,2) > 0.01
        j = i+1;
        while (j <= t) && (fq(j,2) > 0.01)
            j = j+1;
        end
        j = j-1;
        if j == t
            n = t-i+1;
            R = linspace(rads(i-1),0.5*rads(i-1),n+1);
            R = R(2:end);
            St = Starts(end,:)+Lengs(end)*Axes(end,:);
        else
            n = j-i+1;
            R = linspace(rads(i-1),rads(j+1),n+2);
            R = R(2:end-1);
            St = starts(j+1,:);
        end
        Sb = starts(i-1,:)+lengs(i-1)*axes(i-1,:);
        A = St-Sb;
        L = norm(A);
        A = A/L;
        S = zeros(n,3);
        for k = 1:n
            S(k,:) = Sb+(k-1)/n*L*A;
        end
        
        rads(i:j) = R;
        lengs(i:j) = L/n;
        axes(i:j,:) = repmat(A,[n,1]);
        starts(i:j,:) = S;
        
        i = j+1;
    else
        i = i+1;
    end
end

FQ = fq;
Rads = rads;
Lengs = lengs;
Axes = axes;
Starts = starts;

end % End of function 


function [PC,Starts,Lengs,Axes,Segs,SPar,SChi,CExt,t,ChaOfSegs] = parent_cylinder...
    (P,Regs,Segs,SPar,SChi,CiS,CExt,Rad,Axe,Sta,Len,Starts,Axes,Rads,Lengs,c,si,t,ChaOfSegs)

% Finds the parent cylinder from the possible parent segment.
% Does this by checking if the axis of the cylinder, if continued, will
% cross the nearby cylinders in the parent segment.
% Adjust the cylinder so that it starts from the surface of its parent.

% PC     Parent cylinder

if SPar(si) > 0 % parent segment exists, find the parent cylinder
    s = SPar(si);
    PC = CiS{s}; % the cylinders in the parent segment
    
    if ~isempty(PC) 
        % select the closest cylinders for closer examination
        if length(PC) > 3
            D = mat_vec_subtraction(-Sta(PC,:),-Starts(1,:));
            d = sum(D.*D,2);
            [~,I] = sort(d);
            I = I(1:4);
            pc = PC(I);
            [~,pcMax] = max(pc);
        else
            pc = PC;
            [~,pcMax] = max(pc);
        end
        
        
        %% Check if the cylinder is the extension of its parent cylinder
        ParentFound = false;  % true when the parent cylinder is found
        % Cylinder must be nearly parallel and the parent must not have an
        % extension
        if (pcMax == pc(1)) && (abs(Axe(pc(1),:)*Axes(1,:)') > 0.7) && (CExt(pc(1)) == 0)
            
            % cylinder must be situated as an extension compared to the parent 
            Q = [Starts(1,:); Starts(1,:)+Lengs(1)*Axes(1,:)];
            [d,~,h] = distances_to_line(Q,Axe(pc(1),:),Sta(pc(1),:));
            if (d(1) < Rad(pc(1))+Rads(1)) && (h(1) > 0.75*Len(pc(1))) && (h(2) > Len(pc(1)))
                
                % Adjusted cylinder, which continues/extends its parent,
                % must be nearly parallel
                PC = pc(1);
                E = Starts(1,:)+Lengs(1)*Axes(1,:); % the top of the cylinder
                S = Sta(PC,:)+Len(PC)*Axe(PC,:); % the top of the parent
                V = E-S;    % new axis for the cylinder
                L = norm(V);
                V = V/L;
                if abs(V*Axe(pc(1),:)') > 0.7
                    % cylinder extends its parent
                    Lengs(1) = L;
                    Starts(1,:) = S;
                    Axes(1,:) = V;
                    CExt(PC) = c;
                    S = SPar(si);
                    Segs{S} = [Segs{S}; Segs{si}];
                    Segs{si} = zeros(0);
                    SPar(si) = 0;
                    SChi{S} = setdiff(SChi{S},si);
                    SChi{S} = [SChi{S}; SChi{si}];
                    SPar(SChi{si}) = S;
                    SChi{si} = zeros(0);
                    ParentFound = true;
                    ChaOfSegs = true;
                end
            end
        end
        
        %% Check if a small segment and its cylinders are part of its parent segment
        % These segments and cylinders are deleted
        n = length(pc);   % number of cylinders to be checked
        if (SPar(si) > 0) && ~ParentFound && (t <= 3)
            S = Regs{1};
            j = 1;
            while (j <= n) && ~ParentFound
                [d,~,h] = distances_to_line(P(S,:),Axe(pc(j),:),Sta(pc(j),:));
                mh = mean(h);
                if (mean(d) < 1.1*Rad(pc(j))) && (mh < Len(pc(j))) && (mh > 0)
                    ParentFound = true;
                    Segs{si} = zeros(0); % remove the segment
                    S = SPar(si);
                    SPar(si) = 0;
                    SChi{S} = setdiff(SChi{S},si);
                    SChi{S} = [SChi{S}; SChi{si}];
                    SPar(SChi{si}) = S;
                    SChi{si} = zeros(0);
                    t = 0;  % remove the cylinders
                    PC = pc(j);
                end
                j = j+1;
            end
        end
        
        %% Calculate the relations of the cylinder's start and end points to the possible parents
        if ~ParentFound
            Q = [Starts(1,:); Starts(1,:)+Lengs(1)*Axes(1,:)]; % start and end points
            D = zeros(n,2);  % distances
            H = zeros(n,2);  % heights of the closest axis point
            W = zeros(n,6);  % vectors between Q and the closest axis points
            AP = zeros(n,6);  % axis points closest to Q
            for j = 1:n
                [d,V,h,ap] = distances_to_line(Q,Axe(pc(j),:),Sta(pc(j),:));
                D(j,:) = d';
                H(j,:) = h';
                W(j,1:3) = V(1,:);
                W(j,4:6) = V(2,:);
                AP(j,1:3) = Sta(pc(j),:)+ap(1,:);
                AP(j,4:6) = Sta(pc(j),:)+ap(2,:);
            end
        end
        
        %% Check if the cylinder starts inside its parent and is inside the extended parent
        j = 1;
        while (j <= n) && ~ParentFound
            R = Rad(pc(j));  % radius of the possible parent
            L = Len(pc(j));  % length of the possible parent 
            if (t <= 3) && (isempty(SChi{si})) && (D(j,1) < R) && (D(j,2) < R) && (H(j,1) > 0) && (H(j,1) < L) && (H(j,2) > 0) && (H(j,2) < L)
                % cylinder's axis is completely inside its parent
                % discard the segment's cylinders
                PC = zeros(0);
                Segs{si} = zeros(0);
                S = SPar(si);
                SPar(si) = 0;
                SChi{S} = setdiff(SChi{S},si);
                SChi{S} = [SChi{S}; SChi{si}];
                SPar(SChi{si}) = S;
                SChi{si} = zeros(0);
                t = 0;
                ParentFound = true;
            elseif (D(j,1) < R) && (D(j,2) < R) && (H(j,1) > 0) && (H(j,1) < L)
                % cylinder's axis starts inside its parent and the end 
                % point is inside the extended parent
                if (j == pcMax) && (H(j,2) > L) && (H(j,1) > 0.5*L) && (abs(Axe(pc(j),:)*Axes(1,:)') > 0.7) && (CExt(pc(j)) == 0)
                    % cylinder extends its parent
                    PC = pc(j);
                    E = Starts(1,:)+Lengs(1)*Axes(1,:);
                    S = Sta(PC,:)+Len(PC)*Axe(PC,:);
                    V = E-S;
                    Lengs(1) = norm(V);
                    Starts(1,:) = S;
                    Axes(1,:) = V/Lengs(1);
                    CExt(PC) = c;
                    S = SPar(si);
                    Segs{S} = [Segs{S}; Segs{si}];
                    Segs{si} = zeros(0);
                    SPar(si) = 0;
                    SChi{S} = setdiff(SChi{S},si);
                    SChi{S} = [SChi{S}; SChi{si}];
                    SPar(SChi{si}) = S;
                    SChi{si} = zeros(0);
                    ParentFound = true;
                    ChaOfSegs = true;
                elseif (t < 3) && (isempty(SChi{si}))
                    % cylinder is not the extension of its parent, 
                    % discard the segment's cylinders
                    PC = zeros(0);
                    Segs{si} = zeros(0);
                    S = SPar(si);
                    SPar(si) = 0;
                    SChi{S} = setdiff(SChi{S},si);
                    SChi{S} = [SChi{S}; SChi{si}];
                    SPar(SChi{si}) = S;
                    SChi{si} = zeros(0);
                    t = 0;
                    ParentFound = true;
                else
                    PC = pc(j);
                    ParentFound = true;
                end
            end
            j = j+1;
        end
        
        
        %% Check if the cylinder is nearby its parent but not "normal"
        j = 1;
        while (j <= n) && ~ParentFound
            R = Rad(pc(j));
            L = Len(pc(j));
            if (D(j,1) < R) && (D(j,2) < R+Rads(1)) && (H(j,1) > 0) && (H(j,1) < L) && (H(j,2) > 0) && (H(j,2) < L)  && (isempty(SChi{si}))
                % cylinder's axis starts inside its parent and the end point
                % is little outside its parent but cylinder is "between"
                % its parent's start and end points remove the cylinders 
                PC = zeros(0);
                Segs{si} = zeros(0);
                S = SPar(si);
                SPar(si) = 0;
                SChi{S} = setdiff(SChi{S},si);
                SChi{S} = [SChi{S}; SChi{si}];
                SPar(SChi{si}) = S;
                SChi{si} = zeros(0);
                t = 0;
                ParentFound = true;
            elseif (D(j,1) < R+Rads(1)) && (D(j,2) < R) && (H(j,1) > 0) && (H(j,1) < L)
                % cylinder's axis starts little outside its parent but the 
                % end point is inside its extended parent
                if (j ~= pcMax) && (t < 3) && (isempty(SChi{si}))
                    % cylinder is not the extension of its parent, remove
                    PC = zeros(0);
                    Segs{si} = zeros(0);
                    S = SPar(si);
                    SPar(si) = 0;
                    SChi{S} = setdiff(SChi{S},si);
                    SChi{S} = [SChi{S}; SChi{si}];
                    SPar(SChi{si}) = S;
                    SChi{si} = zeros(0);
                    t = 0;
                elseif (j == pcMax) && (H(j,2) > L) && (abs(Axe(pc(j),:)*Axes(1,:)') > 0.7) && (CExt(pc(j)) == 0)
                    % cylinder is the extension of its parent
                    PC = pc(j);
                    E = Starts(1,:)+Lengs(1)*Axes(1,:);
                    S = Sta(PC,:)+Len(PC)*Axe(PC,:);
                    V = E-S;
                    Lengs(1) = norm(V);
                    Starts(1,:) = S;
                    Axes(1,:) = V/Lengs(1);
                    CExt(PC) = c;
                    S = SPar(si);
                    Segs{S} = [Segs{S}; Segs{si}];
                    Segs{si} = zeros(0);
                    SPar(si) = 0;
                    SChi{S} = setdiff(SChi{S},si);
                    SChi{S} = [SChi{S}; SChi{si}];
                    SPar(SChi{si}) = S;
                    SChi{si} = zeros(0);
                    ChaOfSegs = true;
                else
                    PC = pc(j);
                end
                ParentFound = true;
            elseif (D(j,1) > R) && (D(j,2) < R+Rads(1)) && (H(j,1) > 0) && (H(j,1) < L) && (Rads(1) < R)
                % cylinder starts outside its parent and cylinder's end is
                % little outside the extended parent, the radius is also
                % smaller than its parent's radius
                if (j ~= pcMax) && (t < 3) && (isempty(SChi{si}))
                    % the segment is small and the cylinder is not the
                    % extension of its parent, remove
                    PC = zeros(0);
                    Segs{si} = zeros(0);
                    S = SPar(si);
                    SPar(si) = 0;
                    SChi{S} = setdiff(SChi{S},si);
                    SChi{S} = [SChi{S}; SChi{si}];
                    SPar(SChi{si}) = S;
                    SChi{si} = zeros(0);
                    t = 0;
                elseif (j == pcMax) && (H(j,2) > L) && (abs(Axe(pc(j),:)*Axes(1,:)') > 0.7) && (CExt(pc(j)) == 0)
                    % cylinder is the extension of its parent
                    PC = pc(j);
                    E = Starts(1,:)+Lengs(1)*Axes(1,:);
                    S = Sta(PC,:)+Len(PC)*Axe(PC,:);
                    V = E-S;
                    Lengs(1) = norm(V);
                    Starts(1,:) = S;
                    Axes(1,:) = V/Lengs(1);
                    CExt(PC) = c;
                    S = SPar(si);
                    Segs{S} = [Segs{S}; Segs{si}];
                    Segs{si} = zeros(0);
                    SPar(si) = 0;
                    SChi{S} = setdiff(SChi{S},si);
                    SChi{S} = [SChi{S}; SChi{si}];
                    SPar(SChi{si}) = S;
                    SChi{si} = zeros(0);
                    ChaOfSegs = true;
                else
                    PC = pc(j);
                end
                ParentFound = true;
            elseif (t == 1) && (abs(Axe(pc(j),:)*Axes(1,:)') > 0.7) && ((D(j,1) < R+Rads(1)) || (D(j,2) < R+Rads(1))) ...
                    && (((H(j,1) > 0) && (H(j,1) < L)) || ((H(j,2) > 0) && (H(j,2) < L))) && (isempty(SChi{si}))
                % one cylinder which is nearly parallel to its parent and
                % is nearby its parent, delete
                PC = zeros(0);
                Segs{si} = zeros(0);
                S = SPar(si);
                SPar(si) = 0;
                SChi{S} = setdiff(SChi{S},si);
                SChi{S} = [SChi{S}; SChi{si}];
                SPar(SChi{si}) = S;
                SChi{si} = zeros(0);
                t = 0;
                ParentFound = true;
            end
            j = j+1;
        end
        
        
        %% Check possible crossing points
        if ~ParentFound
            % Calculate the possible crossing points of the cylinder axis, when
            % extended, on the surfaces of the parent candidate cylinders
            x = zeros(n,2);  % how much the starting point has to move to cross
            h = zeros(n,2);  % the crossing point height in the parent
            for j = 1:n
                % Crossing points solved from a quadratic equation
                A = Axes(1,:)-(Axes(1,:)*Axe(pc(j),:)')*Axe(pc(j),:);
                B = Starts(1,:)-Sta(pc(j),:)-(Starts(1,:)*Axe(pc(j),:)')*Axe(pc(j),:)...
                    +(Sta(pc(j),:)*Axe(pc(j),:)')*Axe(pc(j),:);
                e = A*A';
                f = 2*A*B';
                g = B*B'-Rad(pc(j))^2;
                di = sqrt(f^2 - 4*e*g);  % the discriminant
                x(j,1) = (-f + di)/(2*e);  % how much must the starting point be moved to cross
                x(j,2) = (-f - di)/(2*e);
                if isreal(x(j,1)) %% cylinders can cross
                    % the heights of the crossing points
                    h(j,1) = Starts(1,:)*Axe(pc(j),:)'+x(j,1)*Axes(1,:)*Axe(pc(j),:)'-...
                        Sta(pc(j),:)*Axe(pc(j),:)';
                    h(j,2) = Starts(1,:)*Axe(pc(j),:)'+x(j,2)*Axes(1,:)*Axe(pc(j),:)'-...
                        Sta(pc(j),:)*Axe(pc(j),:)';
                end
            end
            % Check first the clear cut cases
            % sp = starting points
            j = 1;
            while j <= n && ~ParentFound
                if x(j,1) > 0 && x(j,2) < 0 && h(j,1) > 0 && h(j,1) < Len(pc(j))
                    % sp inside the parent and crosses its surface
                    PC = pc(j);
                    if t == 1 && Rads(1) < 0.05*Rad(PC) && (Lengs(1)-x(j,1))/Rads(1) < 0.5 && (isempty(SChi{si}))
                        % one short cylinder which has small radius, remove
                        t = 0;
                        Segs{si} = zeros(0);
                        S = SPar(si);
                        SPar(si) = 0;
                        SChi{S} = setdiff(SChi{S},si);
                        SChi{S} = [SChi{S}; SChi{si}];
                        SPar(SChi{si}) = S;
                        SChi{si} = zeros(0);
                        PC = zeros(0);
                    elseif  Lengs(1)-x(j,1) > 0
                        Starts(1,:) = Starts(1,:)+x(j,1)*Axes(1,:);
                        Lengs(1) = Lengs(1)-x(j,1);
                        ParentFound = true;
                    end
                elseif x(j,1) < 0  &&  x(j,2) > 0  &&  h(j,2) > 0  &&  h(j,2) < Len(pc(j))  &&  Lengs(1)-x(j,2) > 0
                    % sp inside the parent and crosses its surface
                    PC = pc(j);
                    Starts(1,:) = Starts(1,:)+x(j,2)*Axes(1,:);
                    Lengs(1) = Lengs(1)-x(j,2);
                    ParentFound = true;
                elseif x(j,1) < 0 && x(j,2) < 0 && x(j,2) < x(j,1) && h(j,1) > 0 && h(j,1) < Len(pc(j)) && Lengs(1)-x(j,1) > 0
                    % sp outside the parent and crosses its surface when extended
                    % backwards
                    PC = pc(j);
                    Starts(1,:) = Starts(1,:)+x(j,1)*Axes(1,:);
                    Lengs(1) = Lengs(1)-x(j,1);
                    ParentFound = true;
                elseif x(j,1) < 0 && x(j,2) < 0 && x(j,2) > x(j,1) && h(j,2) > 0 && h(j,2) < Len(pc(j)) && Lengs(1)-x(j,2) > 0
                    % sp outside the parent and crosses its surface when extended
                    % backwards
                    PC = pc(j);
                    Starts(1,:) = Starts(1,:)+x(j,2)*Axes(1,:);
                    Lengs(1) = Lengs(1)-x(j,2);
                    ParentFound = true;
                elseif x(j,1) > 0 && x(j,2) > 0 && x(j,2) < x(j,1) && h(j,1) > 0 && h(j,1) < Len(pc(j)) && Lengs(1)-x(j,1) > 0
                    % sp outside the parent but crosses its surface when extended
                    % forward
                    PC = pc(j);
                    Starts(1,:) = Starts(1,:)+x(j,1)*Axes(1,:);
                    Lengs(1) = Lengs(1)-x(j,1);
                    ParentFound = true;
                elseif x(j,1) > 0 && x(j,2) > 0 && x(j,2) > x(j,1) && h(j,2) > 0 && h(j,2) < Len(pc(j)) && Lengs(1)-x(j,2) > 0
                    % sp outside the parent and crosses its surface when extended
                    % forward
                    PC = pc(j);
                    Starts(1,:) = Starts(1,:)+x(j,2)*Axes(1,:);
                    Lengs(1) = Lengs(1)-x(j,2);
                    ParentFound = true;
                end
                j = j+1;
            end
        end
        
        
        %% Check if a parent can still be found 
        if ~ParentFound
            j = 1;
            while (j <= n) && ~ParentFound
                if (H(j,1) > 0) && (H(j,1) < Len(pc(j)))
                    % cylinder is "between" its parent's start and end
                    % points
                    PC = pc(j);
                    V = W(j,1:3)/norm(W(j,1:3));
                    Starts(1,:) = AP(j,1:3)+Rad(PC)*V;
                    V = Q(2,:)-Starts(1,:);
                    L = norm(V);
                    V = V/L;
                    Axes(1,:) = V;
                    Lengs(1) = L;
                    ParentFound = true;
                end
                j = j+1;
            end
            
            j = 1;
            while (j <= n) && ~ParentFound
                % the parent is the one for which the cylinder is like
                % extension
                if (j == pcMax) && (H(j,1) > 0.75*Len(pc(j)))
                    PC = pc(j);
                    ParentFound = true;
                end
                j = j+1;
            end
            
            % The parent is the closest cylinder in the parent segment
            if ~ParentFound
                PC = pc(1);
            end
        end
    end
    
else
    % no parent segment exists
    PC = zeros(0);
end

end % End of function


function [Rads,Starts,Axes,Lengs,t] = adjustments(P,Rad,Sta,CiS,Rads,...
    Starts,Axes,Lengs,Regs,FQ,Width,Height,Dt,Mt,PC,si,t,SoC,SChi,SPar)


%% Determine the maximum radius
if ~isempty(PC)
    MaxR = Rad(PC);
    MaxR = max(MaxR,0.001);
elseif si == 1
    % For the trunk use 3 times the mean as the indicative maximum value
    MaxR = 3*mean(Rads);
elseif t > 0
    % Determine the maximum radius by position
    C = CiS{1};
    n = length(C);
    if n > 0
        if n > 6
            n = n-4;
        end
        j = 1;
        h = Starts(1,3);
        while (j <= n) && (Sta(j,3) < h)
            j = j+1;
        end
        if j > n
            MaxR = (1-(h-Sta(n,3))/(Height-Sta(n,3)))*Rad(n);
        else
            MaxR = Rad(j);
        end
        
        d = distances_to_line(Starts(1,:),Dt,Mt);
        MaxR = (Width-d)/Width*MaxR/3;
        
        MaxR = max([0.005 MaxR]);
    else
        MaxR = 0.005;
    end
else
    MaxR = max(Rads);
    MaxR = max(MaxR,0.001);
end

%% Check the radii against the maximum radius
if SPar(si) == 1 && ~isempty(PC)  % branch whose parent cylinder is part of trunk
    rMod = Rads(1:t) > 0.8*MaxR;
    Rads(rMod) = 0.8*MaxR;
elseif (t == 1) && (Rads(1) >= 4*MaxR) % small segment
    if isempty(SChi(si))
        t = 0;
    else
        rMod = Rads(1) > MaxR;
        Rads(rMod) = MaxR;
    end
elseif t > 0
    if t == 1;
        rMod = Rads(1) > 0.8*MaxR;
        Rads(rMod) = 0.8*MaxR;
    elseif ~isempty(PC) && SoC(PC) == 1
        rMod = Rads(1:t) > 0.9*MaxR;
        Rads(rMod) = 0.9*MaxR;
    else
        rMod = Rads(1:t) > MaxR;
        Rads(rMod) = MaxR;
    end
end

%% Check minimum radius
% Set the minimum radius to 1mm
if t > 0
    I = Rads(1:t) < 0.001;
    Rads(I) = 0.001;
    rMod = rMod | I;
end

%% Check the radii change
% Check if big radii changes between succesive cylinders and modify
% accordingly
if t > 1
    j = 2;
    while (j <= t)
        I1 = Rads(j) > 1.1*Rads(j-1);
        I2 = Rads(j) < 0.75*Rads(j-1);
        I3 = FQ(j,3) < 0.3125;
        I4 = FQ(j,3) < FQ(j-1,3);
        I5 = FQ(j,4) == 0;
        I6 = FQ(j,2)/Rads(j) > 0.05;
        I7 = FQ(j,2) > 2*FQ(j-1,2);
        I8 = Rads(j) > 1.25*Rads(j-1);
        I9 = FQ(j,2)/Rads(j) > 0.01;
        if (j > 3) && (Rads(j-2) < 3*Rads(j-1)) && (Rads(j) > 3*Rads(j-1))
            k = j-1;
            while (j <= t) && (Rads(j) > 3*Rads(k))
                j = j+1;
            end
            if j == t+1
                S = Starts(k,:)+Lengs(k)*Axes(k,:);
                E = Starts(t,:)+Lengs(t)*Axes(t,:);
                V = E-S;
                n = j-k-1;
                L = norm(V);
                for i = 1:n
                    Starts(k+i,:) = S+(i-1)/n*V;
                    Axes(k+i,:) = V/L;
                    Lengs(k+i) = L/n;
                    Rads(k+i) = Rads(k);
                end
            else
                S = Starts(k,:)+Lengs(k)*Axes(k,:);
                E = Starts(j,:);
                V = E-S;
                n = j-k-1;
                R = linspace(Rads(k),Rads(j),n+2);
                L = norm(V);
                for i = 1:n
                    Starts(k+i,:) = S+(i-1)/n*V;
                    Axes(k+i,:) = V/L;
                    Lengs(k+i) = L/n;
                    Rads(k+i) = R(i+1);
                end
            end
        elseif (I1 || I2) && I3 && I4 && I5 && I6
            k = 1;
            while (j+k <= t) && (Rads(j+k) > 1.25*Rads(j-1) || ...
                    (FQ(j+k,3) < 0.3125 && FQ(j+k,2) > FQ(j-1,2) && FQ(j+k,3) < FQ(j-1,3) && FQ(j+k,4) == 0))
                k = k+1;
            end
            if (j+k <= t) && (FQ(j+k,2)/Rads(j+k) <= 0.05) && ((FQ(j+k,3) >= 0.25) || (FQ(j+k,4) == 1))
                R = linspace(Rads(j-1),Rads(j+k),k+1);
                Rads(j:j+k-1) = R(2:end);
                S = Starts(j-1,:)+Lengs(j-1)*Axes(j-1,:);
                V = Starts(j+k,:)-S;
                if k == 1
                    Starts(j,:) = S;
                    Axes(j,:) = V/norm(V);
                    Lengs(j) = norm(V);
                else
                    L = norm(V);
                    v = V/L;
                    for m = 0:k-1
                        Starts(j+m,:) = S+m/k*V;
                        Axes(j+m,:) = v;
                        Lengs(j+m) = L/k;
                    end
                end
                j = j+k;
            else
                if j+k <= t
                    R = linspace(Rads(j-1),Rads(j+k),k+1);
                    Rads(j:j+k-1) = R(2:end);
                    rMod(j:j+k-1) = true;
                    j = j+k;
                elseif j == t
                    Rads(j) = Rads(j-1);
                    rMod(j) = true;
                    j = t+1;
                else
                    Rads(j:end) = Rads(j-1);
                    rMod(j:end) = true;
                    j = j+k;
                end
            end
        elseif (I1 || I2) && I5 && I6 && I7
            Rads(j) = Rads(j-1);
            rMod(j) = true;
            j = j+1;
        elseif I8 && I9 && I4 && (FQ(j-1,4) == 1) && (FQ(j-1,3) >= 0.25)
            Rads(j) = Rads(j-1);
            rMod(j) = true;
            j = j+1;
        elseif (j == t) && (Rads(j) > Rads(j-1))
            % Change always the last radius, if it is larger than the preceding one 
            Rads(j) = Rads(j-1);
            rMod(j) = true;
            j = j+1;
        elseif (j < t) && (I1||I2) && (FQ(j,2)/Rads(j) > 0.2) && (FQ(j-1,2)/Rads(j-1) < 0.05) && (FQ(j+1,2)/Rads(j+1) < 0.05)
            Rads(j) = 0.5*(Rads(j-1)+Rads(j+1));
            rMod(j) = true;
            j = j+1;
        else
            j = j+1;
        end
    end
end

if t > 2
    for i = 2:t-1
        if (Rads(i) > 1.25*Rads(i-1)) && (Rads(i) > 1.25*Rads(i+1))
            Rads(i) = 0.5*(Rads(i-1)+Rads(i+1));
            rMod(i) = true;
        end
    end
end

if t > 2
    for i = 2:t-1
        if (i < t) && (Rads(i) < 0.9*Rads(i-1)) && (Rads(i) < 0.9*Rads(i+1)) && (Rads(i+1) < Rads(i-1))
            Rads(i) = mean([Rads(i-1) Rads(i+1)]);
            rMod(i) = true;
        end
    end
end

%% If radii changes, adjust starting points and possibly other param.
if (t > 1) && any(rMod)
    if rMod(1) && ~rMod(2)
        Starts(1,:) = Starts(2,:)-Lengs(1)*Axes(1,:);
    elseif (t > 2) && rMod(1) && rMod(2) && ~rMod(3)
        Starts(1,:) = Starts(3,:)-Lengs(2)*Axes(2,:)-Lengs(1)*Axes(1,:);
        rMod(1) = false;
    end
    for j = 2:t-1
        if rMod(j)
            if ~rMod(j-1) && ~rMod(j+1)
                E = Starts(j-1,:)+Lengs(j-1)*Axes(j-1,:);
                V = Starts(j+1,:)-E;
                L = norm(V);
                Starts(j,:) = E;
                Axes(j,:) = V/L;
                Lengs(j) = L;
            elseif ~rMod(j-1) && rMod(j+1)
                Starts(j,:) = Starts(j-1,:)+Lengs(j-1)*Axes(j-1,:);
                rMod(j) = false;
            elseif rMod(j-1) && ~rMod(j+1);
                Starts(j,:) = Starts(j+1,:)-Lengs(j)*Axes(j,:);
                rMod(j) = false;
            end
        end
    end
    if (t > 1) && rMod(t) && ~rMod(t-1)
        E = Starts(t,:)+Lengs(t)*Axes(t,:);
        S = Starts(t-1,:)+Lengs(t-1)*Axes(t-1,:);
        V = E-S;
        L = norm(V);
        Starts(t,:) = S;
        Lengs(t) = L;
        Axes(t,:) = V/L;
    end
elseif (t == 1) && any(rMod)
    p = P(Regs{1},:);
    q = mean(p);
    V = mat_vec_subtraction(p,q);
    H = V*Axes(1,:)';
    h = min(H);
    Lengs(1) = max(H)-h;
    Starts(1,:) = q+h*Axes(1,:);
end

%% Continuous branches
% Make cylinders properly "continuous" by moving the starting points
% First check, move the starting point to the plane defined by parent
% cylinder's top
if t > 1
    for j = 2:t
        U = Starts(j,:)-Starts(j-1,:)-Lengs(j-1)*Axes(j-1,:);
        if (norm(U) > 0.0001)
            % First define vector V and W which are orthogonal to the
            % cylinder axis N
            N = Axes(j,:)';
            if norm(N) > 0
                [V,W] = orthonormal_vectors(N);
                % Now define the new starting point
                %warning off
                x = [N V W]\U';
                %warning on
                Starts(j,:) = Starts(j,:)-x(1)*N';
                if x(1) > 0
                    Lengs(j) = Lengs(j)+x(1);
                elseif Lengs(j)+x(1) > 0
                    Lengs(j) = Lengs(j)+x(1);
                end
            end
        end
    end
end

% Make cylinders continuous by moving the starting points.
% Second check, "straighten" cylinders with large difference between the
% end point and the extension's start point
if t > 1
    for j = 1:t-1
        E = Starts(j,:)+Lengs(j)*Axes(j,:);
        S = Starts(j+1,:);
        if norm(E-S) > Rads(j)/2
            if j+2 <= t
                V = Starts(j+2,:)-S;
                W = Starts(j+2,:)-E;
                if norm(V) > norm(W)
                    Lengs(j+1) = norm(W);
                    Axes(j+1,:) = W/Lengs(j+1);
                    Starts(j+1,:) = E;
                else
                    V = S-Starts(j,:);
                    Lengs(j) = norm(V);
                    Axes(j,:) = V/Lengs(j);
                end
            else
                
            end
        end
    end
end

%% Check if NaN
t0 = t;
A = isnan(Axes);
if any(A)
    j = 1;
    while (j <= t) && ~any(A(j,:))
        j = j+1;
    end
    if j <= t
        t = j-1;
    end
else
    A = isnan(Starts);
    if any(A)
        j = 1;
        while (j <= t) && ~any(A(j,:))
            j = j+1;
        end
        if j <= t
            t = j-1;
        end
    else
        A = isnan(Rads);
        if any(A)
            j = 1;
            while (j <= t) && ~any(A(j,:))
                j = j+1;
            end
            if j <= t
                t = j-1;
            end
        else
            A = isnan(Lengs);
            if any(A)
                j = 1;
                while (j <= t) && ~any(A(j,:))
                    j = j+1;
                end
                if j <= t
                    t = j-1;
                end
            end
        end
    end
end

%% Adjust variable lenghts
if t0 > t
    Axes = Axes(1:t,:);
    Starts = Starts(1:t,:);
    Rads = Rads(1:t);
    Lengs = Lengs(1:t);
end

end % End of function

