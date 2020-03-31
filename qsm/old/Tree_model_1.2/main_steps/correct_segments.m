function [Segs,SPar,SChi,SDir] = correct_segments(P,Cen,Segs,SPar,SChi,SDir,Bal,dmin)

% ---------------------------------------------------------------------
% CORRECT_SEGMENTS.M        Combines and removes segments.
%
% Version 1.2
% Author        Pasi Raumonen
% Created       2 December 2013
%
% This software cannot be distributed without the permission of P. Raumonen
% and it is only for non-commercial use. These restrictions holds also for
% the modules or subprograms of the software:
% DISTANCES_TO_LINE.M, MODIFY_SEGMENTS, SEGMENT_DIRECTION
% ---------------------------------------------------------------------

% Combine segments to form new segments. If a small segment is along a
% larger segment as a child segment, then remove the small one. If a
% child segment continues its parent segment, then combine theM to a single
% longer segment. Modify the segments so that the trunk segment is defined
% from the top to bottom.

% Inputs:
% P         Point cloud
% Cen       Indexes of center points
% Segs      Segments
% SPar      Parent segments
% SChi      Child segments
% SDir      Direction of segments at their base
% Bal       Cover sets
% dmin      Minimum cover set diameter

% Outputs:
% Segs      Segments
% SPar      Parent segments
% SChi      Child segments
% SDir      Direction of segments at their base


%% Remove small child segments
ns = size(Segs,1);
I = true(ns,1);
for i = 1:ns
    C = SChi{i};
    if ~isempty(C) % child segments
        n = length(C);
        S = Segs{i};
        s = size(S,1);
        for j = 1:n % check each child separately
            nl = SPar(C(j),2);  % the index of the layer in the parent the child begins
            if nl > 10
                a = nl-10;
            else
                a = 1;
            end
            if s-nl > 10
                b = nl+10;
            else
                b = s;
            end
            B = S{b};
            if length(B) > 1
                B = mean(P(Cen(B),:));
            else
                B = P(Cen(B),:);
            end
            A = S{a};
            if length(A) > 1
                A = mean(P(Cen(A),:));
            else
                A = P(Cen(A),:);
            end
            V = B-A;
            V = V/norm(V);
            Q = vertcat(S{a:b});
            R = Segs{C(j)};
            R = vertcat(R{:});
            R = [Q; R];
            dR = max(distances_to_line(P(Cen(R),:),V,A));
            dQ = max(distances_to_line(P(Cen(Q),:),V,A));
            if (dR-dQ < 0.02) || (dR/dQ < 1.2 && dR-dQ < 0.06)
                K = SChi{C(j)};
                c = length(K);
                if isempty(K)
                    % Remove, no child segments
                    I(C(j)) = false; 
                    Segs{C(j)} = zeros(0,1);
                    SPar(C(j),:) = zeros(1,2);
                    SChi{i} = setdiff(C,C(j));
                else
                    L = SChi(K);
                    L = vertcat(L{:}); % child child segments
                    if isempty(L)
                        J = false(c,1);
                        for k = 1:c
                            r = Segs{K(k)};
                            r = vertcat(r{:});
                            r = [r; R];
                            dr = max(distances_to_line(P(Cen(r),:),V,A));
                            if (dr-dQ < 0.02) || (dr/dQ < 1.2 && dr-dQ < 0.06)
                                J(k) = true;
                            end
                        end
                        if all(J)
                            % Remove
                            K = [K; C(j)];
                            c = length(K);
                            Segs(K) = cell(c,1);
                            I(K) = false; 
                            SPar(K,:) = zeros(c,2);
                            SChi{i} = setdiff(C,C(j));
                            SChi(K) = cell(c,1);
                        end
                    end
                end
            end
        end
    end
end
Segs = Segs(I);
n = nnz(I);
Ind = (1:1:ns)';
J = (1:1:n)';
Ind(I) = J;
Ind(~I) = 0;
SPar = SPar(I,:);
SDir = SDir(I,:);
J = SPar(:,1) > 0;
SPar(J,1) = Ind(SPar(J,1));
% Modify SChi
for i = 1:ns
    if I(i)
        C = SChi{i};
        if ~isempty(C)
            C = nonzeros(Ind(C));
            SChi{i} = C;
        end
    end
end
SChi = SChi(I);
ns = size(Segs,1);


%% Define trunk from top to bottom
nb = size(Bal,1);
SoB = zeros(nb,1);
TopS = [1 min(P(Cen,3))];
for i = 1:ns
    S = Segs{i};
    S = vertcat(S{:});
    SoB(S) = i;
    if max(P(Cen(S),3)) > TopS(2)
        TopS = [i max(P(Cen(S),3))];
    end
end
S = TopS(1);
if S ~= 1
    T = zeros(100,1);
    T(1) = S;
    t = 1;
    while S ~= 1
        S = SPar(S);
        t = t+1;
        T(t) = S;
    end
    T = T(1:t);
    
    % define trunk
    RemSeg = cell(t,1);
    for i = 1:t-2
        I = T(t-i); % segment to be combined to the first segment
        J = T(t-i-1); % above segment's child to be combined next
        SP = SPar(I,2);  % layer index of the child in the parent
        SegP = Segs{1};
        SegC = Segs{I};
        N = size(SegP,1);
        sp = SPar(J,2);  % layer index of the child's child in the child
        if SP >= N-5 % Use the whole parent
            Segs{1} = [SegP; SegC(1:sp)];
            if sp < size(SegC,1) % use only part of the child segment
                Segs{I} = SegC(sp+1:end);
                SPar(I,2) = N+sp;
                
                C = SChi{I};
                K = SPar(C,2) <= sp;
                c = C(~K);
                SChi{I} = c;
                SPar(c,2) = SPar(c,2)-sp+1;
                C = C(K);
                SChi{1} = [SChi{1}; C];
                SPar(C,1) = 1;
                SPar(C,2) = N+SPar(C,2);
                
            else % use the whole child segment
                Segs{I} = cell(0,1);
                RemSeg{i} = I;
                SPar(I,1) = 0;
                
                C = SChi{I};
                SChi{I} = cell(0,1);
                c = setdiff(SChi{1},I);
                if size(c,2) > size(c,1)
                    c = c';
                end
                SChi{1} = [c; C];
                SPar(C,1) = 1;
                SPar(C,2) = N+SPar(C,2);
            end
            
            T(t-i) = 1;
        else % divide the parent segment into two parts
            ns = ns+1;
            Segs{ns,1} = SegP(SP+1:end); % the top part of the parent forms a new segment
            SPar(ns,1) = 1;
            SPar(ns,2) = SP;
            SDir(ns,:) = segment_direction(P,Cen,zeros(1,3),Segs{ns},1);
            Segs{1} = [SegP(1:SP); SegC(1:sp)];
            
            C = SChi{1};
            K = SPar(C,2) > SP;
            SChi{1} = C(~K);
            C = C(K);
            SChi{ns} = C;
            SPar(C,1) = ns;
            SPar(C,2) = SPar(C,2)-SP;
            SChi{1} = [SChi{1}; ns];
            if sp < size(SegC,1) % use only part of the child segment
                Segs{I} = SegC(sp+1:end);
                SPar(I,2) = SP+sp;
                    
                C = SChi{I};
                K = SPar(C,2) <= sp;
                SChi{I} = C(~K);
                C = C(K);
                SChi{1} = [SChi{1}; C];
                SPar(C,1) = 1;
                SPar(C,2) = SP+SPar(C,2);
                
            else % use the whole child segment
                Segs{I} = cell(0);
                RemSeg{i} = I;
                SPar(I,1) = 0;
                
                C = SChi{I};
                c = setdiff(SChi{1},I);
                if size(c,2) > size(c,1)
                    c = c';
                end
                SChi{1} = [c; C];
                SPar(C,1) = 1;
                SPar(C,2) = SP+SPar(C,2);
                
            end
            T(t-i) = 1;
        end
        
%         % plot the progress
%         I = T;
%         for j = 1:t
%             I(j) = T(t-j+1);
%         end
%         S = Segs(I);
%         s = cell(t,1);
%         for j = 1:t
%             A = S{j};
%             if j > 1 && I(j) == I(j-1)
%                 
%             else
%                 s{j} = vertcat(A{:});
%             end
%         end
%         plot_segs(P,s,2,Bal)
%         pause
    end
    
    % Combine the last segment to the trunk
    I = T(1);
    SP = SPar(I,2);
    SegP = Segs{1};
    SegC = Segs{I};
    N = size(SegP,1);
    if SP >= N-5 % Use the whole parent
        Segs{1} = [SegP; SegC];
        Segs{I} = cell(0);
        SPar(I,1) = 0;
        RemSeg{t} = I;
        
        C = SChi{I};
        c = setdiff(SChi{1},I);
        if size(c,2) > size(c,1)
            c = c';
        end
        SChi{1} = [c; C];
        SPar(C,1) = 1;
        SPar(C,2) = N+SPar(C,2);
        
    else % divide the parent segment into two parts
        ns = ns+1;
        Segs{ns,1} = SegP(SP+1:end);
        SPar(ns,1) = 1;
        SDir(ns,:) = segment_direction(P,Cen,zeros(1,3),Segs{ns},1);
        Segs{1} = [SegP(1:SP); SegC];
        Segs{I} = cell(0);
        SPar(I,1) = 0;
        RemSeg{t} = I;
        
        C = SChi{1};
        K = SPar(C,2) > SP;
        SChi{1} = C(~K);
        C = C(K);
        SChi{ns} = C;
        SPar(C,1) = ns;
        SPar(C,2) = SPar(C,2)-SP;
        
        C = SChi{I};
        c = setdiff(SChi{1},I);
        if size(c,2) > size(c,1)
            c = c';
        end
        SChi{1} = [c; C];
        SPar(C,1) = 1;
        SPar(C,2) = SP+SPar(C,2);
        
    end
    
    % Modify indexes by removing empty segments
    RemSeg = vertcat(RemSeg{:});
    E = true(ns,1);
    E(RemSeg) = false;
    Segs = Segs(E);
    Ind = (1:1:ns)';
    n = nnz(E);
    I = (1:1:n)';
    Ind(E) = I;
    SPar = SPar(E,:);
    SDir = SDir(E,:);
    J = SPar(:,1) > 0;
    SPar(J,1) = Ind(SPar(J,1));
    for i = 1:ns
        if E(i)
            C = SChi{i};
            if ~isempty(C)
                C = Ind(C);
                SChi{i} = C;
            end
        end
    end
    SChi = SChi(E);
    ns = n;
    
    % Modify SChi
    for i = 1:ns
        C = SChi{i};
        if size(C,2) > size(C,1) && size(C,1) > 0
            SChi{i} = C';
        elseif size(C,1) == 0 || size(C,2) == 0
            SChi{i} = zeros(0,1);
        end
    end
end


%% Combine segments
ns0 = ns+1;
while ns < ns0
    E = true(ns,1);
    for i = 2:ns
        C = SChi{i};
        S = Segs{i};
        n = size(S,1);
        s = SPar(C,2);
        I = s >= n-1;
        C = C(I);
        if ~isempty(C) && ~isempty(S)
            if n > 5
                S0 = S{n-5};
            else
                S0 = S{1};
            end
            if length(S0) > 1
                Q0 = mean(P(Cen(S0),:));
            else
                Q0 = P(Cen(S0),:);
            end
            S = S{end};
            if length(S) > 1
                Q = mean(P(Cen(S),:));
            else
                Q = P(Cen(S),:);
            end
            V = Q-Q0;
            V = V/norm(V);
            nc = length(C);
            J = false(nc,1);
            q = zeros(nc,3);
            for j = 1:nc
                s = Segs{C(j)};
                if ~isempty(s)
                    if length(s) > 5
                        S0 = s{6};
                    else
                        S0 = s{end};
                    end
                    if length(S0) > 1
                        Q = mean(P(Cen(S0),:));
                    else
                        Q = P(Cen(S0),:);
                    end
                    s = s{1};
                    J(j) = true;
                    if length(s) > 1
                        Q0 = mean(P(Cen(s),:));
                    else
                        Q0 = P(Cen(s),:);
                    end
                    W = Q-Q0;
                    W = W/norm(W);
                    q(j,:) = W;
                end
            end
            if any(J)
                I = q*V';
                [m,k] = max(I);
                C = C(k);
                
                S = Segs{i};
                s = Segs{C};
                Segs{i} = [S; s];
                Segs{C} = zeros(0,1);
                c = SChi{i};
                s = SChi{C};
                SPar(s,:) = [i*ones(length(s),1) size(S,1)+SPar(s,2)];
                c = union(c,s);
                c = setdiff(c,C);
                SChi{i} = c;
                E(C) = false;
            end
        end
    end
    
    % Modify indexes after removing the empty segments
    Segs = Segs(E);
    Ind = (1:1:ns)';
    n = nnz(E);
    I = (1:1:n)';
    Ind(E) = I;
    SPar = SPar(E,:);
    SDir = SDir(E,:);
    J = SPar(:,1) > 0;
    SPar(J,1) = Ind(SPar(J,1));
    for i = 1:ns
        if E(i)
            C = SChi{i};
            if ~isempty(C)
                C = Ind(C);
                SChi{i} = C;
            end
        end
    end
    SChi = SChi(E);
    ns0 = ns;
    ns = nnz(E);
end


%% Modify the base of the segments 
for i = 2:ns
    SegC = Segs{i};
    SP = SPar(i,1);
    SegP = Segs{SP};
    [SegP,Base] = modify_segment(P,Bal,Cen,SDir(SP,:),SegP,SegC,SPar(i,2),SDir(i,:),dmin);
    Segs{SP} = SegP;
end
SPar = SPar(:,1);

% Modify SChi
for i = 1:ns
    C = SChi{i};
    if size(C,2) > size(C,1) && size(C,1) > 0
        SChi{i} = C';
    elseif size(C,1) == 0 || size(C,2) == 0
        SChi{i} = zeros(0,1);
    end
end

str = ['    ',num2str(ns),' segments after correction'];
disp(str)

end % End of main function


function [SegP,Base] = modify_segment(P,Bal,Cen,SD,SegP,SegC,nl,DComp,dmin)

% Expands the base of the branch backwards into its parent segment and
% then removes the expansion from the parent segment.

% Determine the center of Base
Base = SegC{1};
if length(Base) > 1
    B = mean(P(Cen(Base),:));
    db = distances_to_line(P(Cen(Base),:), DComp', B); % distances of the sets in the base to the axis of the branch
    DiamBase = 2*max(db);  % diameter of the base
elseif length(Bal{Base}) > 1
    B = mean(P(Bal{Base},:));
    db = distances_to_line(P(Bal{Base},:), DComp', B);
    DiamBase = 2*max(db);
else
    B = P(Cen(Base),:);
    DiamBase = 0;
end

% Define the segment direction
if nl > 0
    DS = segment_direction(P,Cen,SD,SegP,nl);
else
    DS = zeros(3,1);
end

% Determine the number of cover set layers "n" to be checked
A = abs(DComp*DS);  % abs of cosine of the angle between component and segment directions
n = max([3,ceil(A*2*DiamBase/dmin/2)]);
if n > nl  % can go only to the bottom of the segment
    n = nl;
end

i = 0;
base = cell(n+1,1);
base{1} = Base;
while i < n
    S = SegP{nl-i};
    if length(S) > 1
        q = mean(P(Cen(S),:));
    else
        q = P(Cen(S),:);
    end
    dSets = distances_to_line(P(Cen(S),:), DComp', B); % % distances of the sets in the segment to the axis of the branch
    VBase = mat_vec_subtraction(P(Cen(S),:),B);  % vectors from base's center to sets in the segment
    VSeg = mat_vec_subtraction(P(Cen(S),:),q);  % vectors from segments's center to sets in the segment
    dBase = sqrt(sum(VBase.*VBase,2)); % lengths of VBase
    dSeg = sqrt(sum(VSeg.*VSeg,2)); % lengths of VSeg
    if A < 0.9
        K = dBase < 1.1/(1-0.5*A^2)*dSeg;     % sets closer to the base's center than segment's center
        J = dSets < 1.25*DiamBase;   % sets close enough to the axis of the branch
        I = K&J;
    else % branch almost parallel to parent
        I = dSets < 1.25*DiamBase; % only the distance to the branch axis counts
    end
    
    if all(I) || ~any(I) % stop the proces if all the segment's or no segment's sets
        i = n;
    else
        SegP{nl-i} = S(not(I));
        base{i+1} = S(I);
        i = i+1;
    end
end
Base = vertcat(base{:});
end % End of function


function D = segment_direction(P,Cen,SD,Seg,nl)

% Defines the direction and center of the segment under the study region.

% Direction
if nl-5 > 0
    b = nl-5;
else
    b = 1;
end
if nl+5 <= size(Seg,1)
    t = nl+5;
else
    t = size(Seg,1);
end

if t > b
    B = Seg{b};
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
    V = Top-Bot;
    D = V'/norm(V);
elseif t == b
    D = SD';
else
    D = zeros(3,1);
end

end % End of function
