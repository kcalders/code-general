function [Rad,Len,Axe,Sta,CPar,CExt,CChi,SoC,CiS,Segs,SPar,SChi,Added] = ...
    cylinders(P,Bal,Cen,Segs,SPar,SChi,lcyl)

% ---------------------------------------------------------------------
% CYLINDERS.M       Fits cylinders to the segments
%
% Version 1.4
% Author        Pasi Raumonen
% Created       19 March 2014
%
% This software cannot be distributed without the permission of P. Raumonen
% and it is only for non-commercial use. These restrictions holds also for
% the modules or subprograms of the software:
% FIX_SEGMENTS, REGIONS, PARENT_CYLINDER, ADJUSTMENTS
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


%% Initializations
% Initialize the outputs
NumOfSeg = max(size(Segs));   % number of segments
n = 40*NumOfSeg;
Rad = zeros(n,1);
Len = zeros(n,1);
Axe = zeros(n,3);
Sta = zeros(n,3);
CPar = zeros(n,1);
CExt = zeros(n,1);
CChi = cell(n,1);
CiS = cell(NumOfSeg,1);
SoC = zeros(n,1);
Added = false(n,1);

c = 1;  % number of cylinders determined
Fal = false(size(P,1),1);

% Determine suitable order of segments (from trunk to the "youngest" child)
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
    Pass = true;
    if k > 1
        if isempty(CiS{SPar(si)}) || isempty(Segs{si}) || si == 0
            Pass = false;
        end
    end
    if Pass
        %% Some initialization about the segment
        Seg = Segs{si};      % the current segment under analysis
        Segment = vertcat(Seg{:});
        nl = max(size(Seg));    % number of layers in the segment
        Base = Seg{1};          % the base of the segment
        
        % Check if the base is at the end of the segment
        if nl > 4
            [Seg,Segment,Base,nl] = fix_segment(P,Cen,Seg,Segment,Base,nl);
        end
        
        nc = length(Segment);   % number of cover sets in the current segment
        Q = vertcat(Bal{Segment}); % the points in the segments
        np = length(Q);         % number of points in the segment
        nb = length(vertcat(Bal{Base})); % % number of points in the base

        % Reconstruct only large enough segments
        if nl > 1 && np > nb && nc > 2 && np > 20 && ~isempty(Base) 
            
            %% Determine the regions for cylinder fitting
            [Regs,Starts,Lengs,Axes,Rads,t] = regions(P,Q,Fal,Bal,Seg,nl,Cen,lcyl,si);
            Starts0 = Starts;
            Lengs0 = Lengs;
            Axes0 = Axes;
            
                
            %% Fit cylinders to the regions
            if t > 0
                [Rads,Lengs,Axes,Starts] = ...
                    cylinder_fitting(P,Regs,Rads,Lengs,Axes,Starts);
            end      
            
            
            %% Search possible parent cylinder
            if t > 0 && si > 1
                
                [PC,Rads,Starts,Lengs,Axes,t,added] = parent_cylinder...
                    (SPar,SChi,CiS,Rad,Axe,Sta,Len,Rads,Starts,Axes,Lengs,si,t);

                if added
                    Lengs0 = [Lengs(1); Lengs0];
                    Axes0 = [Axes(1,:); Axes0];
                    Starts0 = [Starts(1,:); Starts0];
                end
            elseif si == 1
                PC = zeros(0,1);
                added = false;
            else
                added = false;
            end
            
            
            %% Adjust cylinders
            if t > 0
                Rads = Rads(1:t);
                Lengs = Lengs(1:t,:);
                Axes = Axes(1:t,:);
                Starts = Starts(1:t,:);
                [Rads,Starts,Axes,Lengs] = adjustments(Rad,Rads,...
                    Starts,Axes,Lengs,Starts0,Axes0,Lengs0,PC,si);
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
                if added
                    Added(c) = true;
                end          
                c = c+t; % number of cylinders so far (plus one)
                
            end
        end
    end
end
c = c-1; % number of cylinders 

%% Finalize outputs
Rad = Rad(1:c);
Len = Len(1:c);
Axe = Axe(1:c,:);
Sta = Sta(1:c,:);
SoC = SoC(1:c);
CPar = CPar(1:c);
CExt = CExt(1:c);
CChi = CChi(1:c,:);
Added = Added(1:c);
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

% Display the number of fitted cylinders
str = ['    ',num2str(c),' cylinders fitted'];
disp(str)
if max(Rad) > 0.05
    I = Rad > 0.01;
    m = median(Len(I)./Rad(I));
else
    m = median(Len./Rad);
end
m = round(100*m)/100;
str = ['    Median relative cylinder length: ',num2str(m)];
disp(str)

%plot_cylinder_model(Rad,Len,Axe,Sta,2,20,0.3)
%pause


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


function [Regs,Starts,Lengs,Axes,Rads,t,NL] = regions...
    (P,Q,Fal,Bal,Seg,nl,Cen,lcyl,si)

% Define the regions for cylinder fitting

if nl > 3
    
    %% Define each region with approximate lenght of lcyl
    Fal(Q) = true;
    H = lcyl;
    Test = vertcat(Seg{1:4});
    Test = vertcat(Bal{Test});
    B = vertcat(Bal{vertcat(Seg{1:2})});
    if length(B) > 1
        B = mean(P(B,:));
    else
        B = P(B,:);
    end
    T = vertcat(Bal{vertcat(Seg{3:4})});
    T = mean(P(T,:));
    V = T-B;
    [d,~,h] = distances_to_line(P(Test,:),V,B);
    [~,~,ht] = distances_to_line(T,V,B);
    I = h <= ht;
    R = median(d(I));
    if si > 1
        J = d < 3*R;
    else
        J = d < 1.5*R;
    end
    I = I&J;
    R = median(d(I));
    L = max(h(I))-min(h);
    i = 4;   i0 = 1;   
    while (i <= nl-1 && L < H*R) || (i <= nl-1 && nnz(I) < 30)
        T = vertcat(Bal{Seg{i}});
        if length(T) > 1
            T = mean(P(T,:));
        else
            T = P(T,:);
        end
        V = T-B;
        i = i+1;
        Test = vertcat(Seg{i0:i});
        if length(Test) < 60
            Test = vertcat(Bal{Test});
        else
            Test = Cen(Test);
        end
        [d,~,h] = distances_to_line(P(Test,:),V,B);
        [~,~,ht] = distances_to_line(T,V,B);
        I = h <= ht;
        R = median(d(I));
        J = d < 3*R;
        I = I&J;
        R = median(d(I));
        L = max(h(I))-min(h);
    end
    Test = Test(I);
    NL = max(i-1,3);
    n = ceil(3*nl/NL);
    Regs = cell(n,1);
    Regs{1} = Test;
    Fal(Test) = false;
    Axes = zeros(n,3);
    Starts = zeros(n,3);
    Rads = zeros(n,1);
    Lengs = zeros(n,1);
    Axes(1,:) = V'/norm(V);
    Starts(1,:) = B+min(h)*Axes(1,:);
    Rads(1) = R;
    Lengs(1) = L;
    t = 1;
    
    i0 = NL-1;
    i = NL+ceil(NL/2);
    while i <= nl-1
        k = ceil(NL/2);
        B = T;
        T = vertcat(Bal{Seg{i}});
        if length(T) > 1
            T = mean(P(T,:));
        else
            T = P(T,:);
        end
        V = T-B;
        i = i+1;
        Test = vertcat(Seg{i0:i});
        if length(Test) < 60
            Test = vertcat(Bal{Test});
        else
            Test = Cen(Test);
        end
        [d,~,h] = distances_to_line(P(Test,:),V,B);
        [~,~,ht] = distances_to_line(T,V,B);
        I = h <= ht;
        J = h >= 0;
        I = I&J;
        if nnz(I) < 3
            I = h >= 0;
        end
        R = median(d(I));
        if R == 0
            R = mean(d(I));
            if R == 0
                R = max(d(I));
                if R == 0
                    R = mad(d);
                end
            end
        end
        J = d < 3*R;
        I = I&J;
        R = median(d(I));
        L = norm(V); %max(h(I));
        k = k+1;
        while (i <= nl-1 && L < H*R && k <= NL) || (i <= nl-1 && nnz(I) < 20)
            T = vertcat(Bal{Seg{i}});
            if length(T) > 1
                T = mean(P(T,:));
            else
                T = P(T,:);
            end
            V = T-B;
            i = i+1;
            Test = vertcat(Seg{i0:i});
            if length(Test) < 60
                Test = vertcat(Bal{Test});
            else
                Test = Cen(Test);
            end
            [d,~,h] = distances_to_line(P(Test,:),V,B);
            [~,~,ht] = distances_to_line(T,V,B);
            I = h <= ht;
            J = h >= 0;
            I = I&J;
            if nnz(I) < 3
                I = h >= 0;
            end
            R = median(d(I));
            if R == 0
                R = mean(d(I));
                if R == 0
                    R = max(d(I));
                    if R == 0
                        R = mad(d);
                    end
                end
            end
            J = d < 3*R;
            I = I&J;
            R = median(d(I));
            L = norm(V);
            k = k+1;
        end
        
        if i >= nl-1
            Test = vertcat(Seg{i0:nl});
            if length(Test) < 60
                Test = vertcat(Bal{Test});
            else
                Test = Cen(Test);
            end
            [d,~,h] = distances_to_line(P(Test,:),V,B);
            I = h >= 0;
            R = median(d(I));
            if R == 0
                R = mean(d(I));
                if R == 0
                    R = max(d(I));
                    if R == 0
                        R = mad(d);
                    end
                end
            end
            J = d < 3*R;
            I = I&J;
            R = median(d(I));
            L = max(h(I));
            Test = Test(I);
            I = Fal(Test);
            Test = Test(I);
            if length(Test) >= 20
                t = t+1;
                Regs{t} = Test;
                Axes(t,:) = V'/norm(V);
                Starts(t,:) = B;
                Rads(t) = R;
                Lengs(t) = L;
            else
                Regs{t} = [Regs{t}; Test];
            end
        else
            Test = Test(I);
            I = Fal(Test);
            Test = Test(I);
            if length(Test) >= 20
                t = t+1;
                Regs{t} = Test;
                Axes(t,:) = V'/norm(V);
                Starts(t,:) = B;
                Rads(t) = R;
                Lengs(t) = L;
            else
                Regs{t} = [Regs{t}; Test];
            end
        end
        i0 = i-1;
        i = i0+ceil(NL/2);
    end
    Axes = Axes(1:t,:);
    V = Starts(2:t,:)-Starts(1:t-1,:);
    L = sqrt(sum(V.*V,2));
    Axes(1:t-1,:) = [V(:,1)./L V(:,2)./L V(:,3)./L];
    Starts = Starts(1:t,:);
    Rads = Rads(1:t);
    Lengs = Lengs(1:t);
    Lengs(1:t-1) = L;
    Regs = Regs(1:t);
    
elseif nl > 1
    %% Define a region for small segments
    % Define the direction
    B = Seg{1};
    if length(B) > 1
        Bot = mean(P(Cen(B),:));
    else
        Bot = P(Cen(B),:);
    end
    T = Seg{nl};
    if length(T) > 1
        Top = mean(P(Cen(T),:));
    else
        Top = P(Cen(T),:);
    end
    Axes = Top-Bot;
    Axes = Axes/norm(Axes);
        
    % Define other outputs
    Regs = cell(1,1);
    Regs{1} = vertcat(Bal{vertcat(Seg{2:end})});
    Starts = mean(P(Regs{:},:));
    if max(size(Starts)) == 3
        [d,~,h] = distances_to_line(P(Regs{1},:),Axes,Starts);
        Lengs = max(h)-min(h);
        Rads = median(d);
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


function [Rads,Lengs,Axes,Starts] = cylinder_fitting...
    (P,Regs,Rads,Lengs,Axes,Starts)
% Fit cylinders to the regions
t = size(Rads,1);
for j = 1:t
    if (length(Regs{j}) > 10) && (norm(Axes(j,:)) > 0) % fit cylinders to large enough subsegs
        
        Points = P(Regs{j},:);  % the coordinate points used for fitting
        
        % Initial estimates
        CA0 = Axes(j,:)';     % cylinder axis
        AP0 = Starts(j,:)'; % mean(Points)';   % point in the cylinder axis
        R0 = Rads(j,1);
        
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
        elseif Accept
            H = Points*CA;
            hmin = min(H);
            L = abs(max(H)-hmin);
            hpoint = CA'*AP;
            SP = AP'-(hpoint-hmin)*CA';
        end
        
        if Accept
            Rads(j) = R;
            Axes(j,:) = CA';
            Lengs(j) = L;
            Starts(j,:) = SP;
        else
            % do not accept least square fittings, use initial
            % estimates
        end
    end
end
end % End of function


function [PC,Rads,Starts,Lengs,Axes,t,added] = parent_cylinder...
    (SPar,SChi,CiS,Rad,Axe,Sta,Len,Rads,Starts,Axes,Lengs,si,t)

% Finds the parent cylinder from the possible parent segment.
% Does this by checking if the axis of the cylinder, if continued, will
% cross the nearby cylinders in the parent segment.
% Adjust the cylinder so that it starts from the surface of its parent.

% PC     Parent cylinder

added = false;
if SPar(si) > 0 % parent segment exists, find the parent cylinder
    s = SPar(si);
    PC = CiS{s}; % the cylinders in the parent segment
    
    % select the closest cylinders for closer examination
    if length(PC) > 1
        D = mat_vec_subtraction(-Sta(PC,:),-Starts(1,:));
        d = sum(D.*D,2);
        [~,I] = sort(d);
        if length(PC) > 3
            I = I(1:4);
        end
        pc = PC(I);
        ParentFound = false;
    elseif length(PC) == 1
        ParentFound = true;
    else
        PC = zeros(0,1);
        ParentFound = true;
        %disp('No cylinders in the parent')
    end
    
    
    %% Check possible crossing points
    if ~ParentFound
        pc0 = pc;
        n = length(pc);
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
            x(j,1) = (-f + di)/(2*e);  % how much the starting point must be moved to cross
            x(j,2) = (-f - di)/(2*e);
            if isreal(x(j,1)) %% cylinders can cross
                % the heights of the crossing points
                h(j,1) = Starts(1,:)*Axe(pc(j),:)'+x(j,1)*Axes(1,:)*Axe(pc(j),:)'-...
                    Sta(pc(j),:)*Axe(pc(j),:)';
                h(j,2) = Starts(1,:)*Axe(pc(j),:)'+x(j,2)*Axes(1,:)*Axe(pc(j),:)'-...
                    Sta(pc(j),:)*Axe(pc(j),:)';
            end
        end
        
        %% Extend to crossing point in the extended parent
        j = 1;
        I = isreal(x(:,1));
        pc = pc0(I);
        x = x(I,:);
        h = h(I,:);
        n = nnz(I);
        X = zeros(n,3);
        while j <= n && ~ParentFound
            if x(j,1) > 0 && x(j,2) < 0
                % sp inside the parent and crosses its surface
                if h(j,1) >= 0 && h(j,1) <= Len(pc(j)) && Lengs(1)-x(j,1) > 0
                    PC = pc(j);
                    Starts(1,:) = Starts(1,:)+x(j,1)*Axes(1,:);
                    Lengs(1) = Lengs(1)-x(j,1);
                    ParentFound = true;
                elseif Lengs(1)-x(j,1) > 0
                    if h(j,1) < 0
                        X(j,:) = [x(j,1) abs(h(j,1)) 0];
                    else
                        X(j,:) = [x(j,1) h(j,1)-Len(pc(j)) 0];
                    end
                else
                    X(j,:) = [x(j,1) h(j,1) 1];
                end
            elseif x(j,1) < 0 && x(j,2) > 0 && Lengs(1)-x(j,2) > 0
                % sp inside the parent and crosses its surface
                if h(j,2) >= 0 && h(j,2) <= Len(pc(j)) && Lengs(1)-x(j,2) > 0
                    PC = pc(j);
                    Starts(1,:) = Starts(1,:)+x(j,2)*Axes(1,:);
                    Lengs(1) = Lengs(1)-x(j,2);
                    ParentFound = true;
                elseif Lengs(1)-x(j,2) > 0
                    if h(j,2) < 0
                        X(j,:) = [x(j,2) abs(h(j,2)) 0];
                    else
                        X(j,:) = [x(j,2) h(j,2)-Len(pc(j)) 0];
                    end
                else
                    X(j,:) = [x(j,2) h(j,2) 1];
                end
            elseif x(j,1) < 0 && x(j,2) < 0 && x(j,2) < x(j,1) && Lengs(1)-x(j,1) > 0
                % sp outside the parent and crosses its surface when extended
                % backwards
                if h(j,1) >= 0 && h(j,1) <= Len(pc(j)) && Lengs(1)-x(j,1) > 0
                    PC = pc(j);
                    Starts(1,:) = Starts(1,:)+x(j,1)*Axes(1,:);
                    Lengs(1) = Lengs(1)-x(j,1);
                    ParentFound = true;
                elseif Lengs(1)-x(j,1) > 0
                    if h(j,1) < 0
                        X(j,:) = [x(j,1) abs(h(j,1)) 0];
                    else
                        X(j,:) = [x(j,1) h(j,1)-Len(pc(j)) 0];
                    end
                else
                    X(j,:) = [x(j,1) h(j,1) 1];
                end
            elseif x(j,1) < 0 && x(j,2) < 0 && x(j,2) > x(j,1) && Lengs(1)-x(j,2) > 0
                % sp outside the parent and crosses its surface when extended
                % backwards
                if h(j,2) >= 0 && h(j,2) <= Len(pc(j)) && Lengs(1)-x(j,2) > 0
                    PC = pc(j);
                    Starts(1,:) = Starts(1,:)+x(j,2)*Axes(1,:);
                    Lengs(1) = Lengs(1)-x(j,2);
                    ParentFound = true;
                elseif Lengs(1)-x(j,2) > 0
                    if h(j,2) < 0
                        X(j,:) = [x(j,2) abs(h(j,2)) 0];
                    else
                        X(j,:) = [x(j,2) h(j,2)-Len(pc(j)) 0];
                    end
                else
                    X(j,:) = [x(j,2) h(j,2) 1];
                end
            elseif x(j,1) > 0 && x(j,2) > 0 && x(j,2) < x(j,1) && Lengs(1)-x(j,1) > 0
                % sp outside the parent but crosses its surface when extended
                % forward
                if h(j,1) >= 0 && h(j,1) <= Len(pc(j)) && Lengs(1)-x(j,1) > 0
                    PC = pc(j);
                    Starts(1,:) = Starts(1,:)+x(j,1)*Axes(1,:);
                    Lengs(1) = Lengs(1)-x(j,1);
                    ParentFound = true;
                elseif Lengs(1)-x(j,1) > 0
                    if h(j,1) < 0
                        X(j,:) = [x(j,1) abs(h(j,1)) 0];
                    else
                        X(j,:) = [x(j,1) h(j,1)-Len(pc(j)) 0];
                    end
                else
                    X(j,:) = [x(j,1) h(j,1) 1];
                end
            elseif x(j,1) > 0 && x(j,2) > 0 && x(j,2) > x(j,1) && Lengs(1)-x(j,2) > 0
                % sp outside the parent and crosses its surface when extended
                % forward
                if h(j,2) >= 0 && h(j,2) <= Len(pc(j)) && Lengs(1)-x(j,2) > 0
                    PC = pc(j);
                    Starts(1,:) = Starts(1,:)+x(j,2)*Axes(1,:);
                    Lengs(1) = Lengs(1)-x(j,2);
                    ParentFound = true;
                elseif Lengs(1)-x(j,2) > 0
                    if h(j,1) < 0
                        X(j,:) = [x(j,2) abs(h(j,2)) 0];
                    else
                        X(j,:) = [x(j,2) h(j,2)-Len(pc(j)) 0];
                    end
                else
                    X(j,:) = [x(j,2) h(j,2) 1];
                end
            end
            j = j+1;
        end
        
        if ~ParentFound && n > 0
            [~,I] = min(X(:,2));
            X = X(I,:);
            if X(3) == 0
                PC = pc(I);
                Starts(1,:) = Starts(1,:)+X(1)*Axes(1,:);
                Lengs(1) = Lengs(1)-X(1);
                ParentFound = true;
            else
                PC = pc(I);
                if t > 1 && X(1) <= Rads(1) && abs(X(2)) <= 1.5*Len(PC)
                    % Remove the first cylinder and adjust the second
                    S = Starts(1,:)+X(1)*Axes(1,:);
                    V = Starts(2)+Lengs(2)*Axes(2,:)-S;
                    Lengs(2) = norm(V);
                    Axes(2,:) = V/norm(V);
                    Starts(2,:) = S;
                    Starts = Starts(2:t,:);
                    Lengs = Lengs(2:t);
                    Axes = Axes(2:t,:);
                    Rads = Rads(2:t);
                    t = t-1;
                elseif t > 1
                    % Remove the first cylinder
                    Starts = Starts(2:t,:);
                    Lengs = Lengs(2:t);
                    Axes = Axes(2:t,:);
                    Rads = Rads(2:t);
                    t = t-1;
                elseif isempty(SChi{si})
                    % Remove the cylinder
                    t = 0;
                    PC = zeros(0,1);
                    %disp('Only one cylinder, no children')
                elseif X(1) <= Rads(1) && abs(X(2)) <= 1.5*Len(PC)
                    % Adjust the cylinder
                    Starts(1,:) = Starts(1,:)+X(1)*Axes(1,:);
                    Lengs(1) = abs(X(1));
                end
                ParentFound = true;
            end
        end
        
        if ~ParentFound
            % The parent is the closest cylinder in the parent segment
            % Add new cylinder
            pc = pc0;
            V = -mat_vec_subtraction(Sta(pc,:),Starts(1,:));
            L = sqrt(sum(V.*V,2));
            [L,I] = min(L);
            PC = pc(I);
            V = V(I,:);
            V = V/L;
            a = acos(V*Axe(PC,:)');
            h = sin(a)*L;
            S = Sta(PC,:)+Rad(PC)/h*L*V;
            L = (h-Rad(PC))/h*L;
            if L > 0.005
                t = t+1;
                Starts = [S; Starts];
                Rads = [Rads(1); Rads];
                Axes = [V; Axes];
                Lengs = [L; Lengs];
                added = true;
            end
        end
    end
else
    % no parent segment exists
    PC = zeros(0,1);
    display('No parent segment')
end

end % End of function


function [Rads,Starts,Axes,Lengs] = adjustments(Rad,Rads,...
    Starts,Axes,Lengs,Starts0,Axes0,Lengs0,PC,si)

t = size(Rads,1);
Rads0 = Rads;
%% Determine the maximum radius
if ~isempty(PC)
    MaxR = 0.95*Rad(PC);
    MaxR = max(MaxR,0.001);
elseif si == 1
    % For the trunk use 3 times the mean as the indicative maximum value
    MaxR = 3*mean(Rads);
else
    MaxR = 0.005;
end
MinR = 0.0025;

% plot_segs(P,Regs,2)
% hold on
% plot_cylinder_model(Rads,Lengs,Axes,Starts,2,20,0.3)
% hold off
% axis equal

%% Check maximum and minimum radii
I = Rads > MaxR;
Rads(I) = MaxR;
I = Rads < MinR;
Rads(I) = MinR;

if max(Rads) < 0.005
    
%     L = [0; cumsum(Lengs)];
%     R = [Rads(1); Rads];
%     figure(1)
%     plot(L,R,'-ob')
    
    %% Adjust radii of thin branches to be linearly decreasing
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
        r = max(Rads);
        if t == 1
            Rads = r;
        else
            Rads = [r; 0.5*r];
        end
    end
    
%     figure(1)
%     R = [Rads(1); Rads];
%     hold on
%     plot(L,R,'-oc')
%     hold off
    
elseif t > 4
    %% Parabola adjustment of maximum and minimum
    % Define parabola taper shape as maximum radii
    % "a" is the number first radii used to determine base radius 
    a = max(ceil(0.1*t),2);
    r = 1.05*mean(Rads(1:a)); % branch base radius
    r0 = MinR;
    l = sum(Lengs(1:t)); % branch length
    L = zeros(t,1); % middle points of cylinder as cumulative length from base
    for i = 1:t
        if i > 1
            L(i) = Lengs(i)/2+sum(Lengs(1:i-1));
        else
            L(i) = Lengs(i)/2;
        end
    end
    
%     figure(112)
%     plot(L,100*Rads,'-ob')
    
    % Determine data "S" for parabola fitting
    b = round(t/4);
    if b >= 3
        %  use 3 first 1/4-length parts to define data points as mean
        %  radii of those cylinders
        S = zeros(5,2);
        S(1,2) = r;
        I0 = 1;
        for i = 1:3
            [~,I] = min(abs(L-i/4*l));
            S(i+1,1) = L(I);
            S(i+1,2) = 1.05*mean(Rads(I0:I));
            I0 = I+1;
        end
    else
        j = 1;
        S = zeros(5,2);
        S(1,2) = r;
        for i = 1:3
            S(i+1,1) = L(round(j+b/2));
            S(i+1,2) = 1.05*mean(Rads(j:j+b-1));
            j = j+b;
        end
    end
    S(5,:) = [l r0];
    
    % Least square fitting of parabola to "S"
    A = [sum(S(:,1).^4) sum(S(:,1).^2); sum(S(:,1).^2) 5];
    y = [sum(S(:,2).*S(:,1).^2); sum(S(:,2))];
    x = A\y;
    R = x(1)*L.^2+x(2); % parabola
    I = Rads > R;
    Rads(I) = R(I);  % change values larger than parabola-values
    Q = 0.75*R;
    I = Q < r0;
    Q(I) = r0;
    I = Rads < Q;
    Rads(I) = Q(I);
    
%     figure(112)
%     hold on
%     plot(S(:,1),100*S(:,2),'*r','Markersize',20)
%     plot(L,100*R,'-or')
%     plot(L,100*Q,'-og')
%     plot(L,100*Rads,'-oc')
%     grid on
%     hold off
%     xlabel('Distance from base (m)')
%     ylabel('Radius (cm)')
%     pause
%     
else
    
%     L = [0; cumsum(Lengs)];
%     R = [Rads(1); Rads];
%     figure(112)
%     plot(L,R,'-ob')

    %% Adjust radii of short branches to be linearly decreasing
    if t > 2
        a = 2*mean(Rads(1:2));
        if a > max(Rads)
            a = max(Rads);
        end
        b = MinR;
        Rads = linspace(a,b,t)';
    else
        r = max(Rads);
        if t == 1
            Rads = r;
        else
            Rads = [r; 0.5*r];
        end
    end
    
%     figure(1)
%     R = [Rads(1); Rads];
%     hold on
%     plot(L,R,'-oc')
%     hold off
end

%% Check big adjustments of starting points
for i = 1:t
    d = abs(Rads(i)-Rads0(i));
    if d > 0.01 || d > 0.5*Rads0(i)
        S0 = Starts0(i,:);
        S = Starts(i,:);
        A0 = Axes0(i,:);
        d = distances_to_line(S,A0,S0);
        if d > 0.01 || d > 0.5*Rads0(i)
            %disp([1 i d Rads0(i)])
            Starts(i,:) = S0;
            Axes(i,:) = A0;
            Lengs(i) = Lengs0(i);
        end
    end
end
% disp('---')
% plot_segs(P,Regs,3)
% hold on
% plot_cylinder_model(Rads,Lengs,Axes,Starts,3,20,0.3)
% hold off
% axis equal

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

% plot_segs(P,Regs,4)
% hold on
% plot_cylinder_model(Rads,Lengs,Axes,Starts,4,20,0.3)
% hold off
% axis equal
% pause


end % End of function

