function [Segs,SChi,SPar] = correct_segments(P,Cen,Segs,SPar,SChi,Bal,dmin,Mod,Add)

% ---------------------------------------------------------------------
% CORRECT_SEGMENTS.M        Combines and removes segments.
%
% Version 2.0
% Author        Pasi Raumonen
% Created       15 March 2014
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

if nargin == 7
    Mod = false;
elseif nargin == 8
    Add = false;
end

% Determine maximum branch order
C = SChi{1};
i = 0;
while ~isempty(C)
    i = i+1;
    C = vertcat(SChi{C});
end
str = ['    Maximum branch order before correction: ',num2str(i)];
disp(str)


%% Remove small child segments
[Segs,SPar,SChi] = remove_small(P,Cen,Segs,SPar,SChi);


%% Make stem and branches as long as possible
[Segs,SPar,SChi] = modify_branches(P,Cen,Bal,Segs,SPar,SChi);


% Determine the maximum branching order
C = SChi{1};
i = 0;
while ~isempty(C)
    i = i+1;
    C = vertcat(SChi{C});
end
str = ['    Maximum branch order after correction:  ',num2str(i)];
disp(str)

if Mod
    %% Modify the base of the segments
    ns = size(Segs,1);
    if Add
        % Add the expanded base to the child and remove it from parent
        for i = 2:ns
            SegC = Segs{i};
            SegP = Segs{SPar(i,1)};
            [SegP,Base] = modify_segment(P,Bal,Cen,SegP,SegC,SPar(i,2),dmin);
            Segs{SPar(i,1)} = SegP;
            SegC{1} = Base;
            Segs{i} = SegC;
        end
    else
        % Only remove the expanded base from parent
        for i = 2:ns
            SegC = Segs{i};
            SegP = Segs{SPar(i,1)};
            SegP = modify_segment(P,Bal,Cen,SegP,SegC,SPar(i,2),dmin);
            Segs{SPar(i,1)} = SegP;
        end
    end
end

SPar = SPar(:,1);

% Modify SChi
for i = 1:ns
    C = SChi{i};
    if size(C,2) > size(C,1) && size(C,1) > 0
        SChi{i} = int32(C');
    elseif size(C,1) == 0 || size(C,2) == 0
        SChi{i} = zeros(0,1);
    else
        SChi{i} = int32(C);
    end
    S = Segs{i};
    for j = 1:size(S,1)
        S{j} = int32(S{j});
    end
    Segs{i} = S;
end

% Display the final number of segments
str = ['    ',num2str(ns),' segments after correction'];
disp(str)

end % End of main function


function [Segs,SPar,SChi] = remove_small(P,Cen,Segs,SPar,SChi)
% Remove small child segments
ns = size(Segs,1);
I = true(ns,1);
for i = 1:ns
    C = SChi{i};  % child segments
    if ~isempty(C) % child segments exists
        n = length(C); % number of children
        S = Segs{i};  % current or parent segment
        s = size(S,1);  % layers in the parent
        for j = 1:n % check each child separately
            nl = SPar(C(j),2);  % the index of the layer in the parent the child begins
            if nl > 10
                a = nl-10; % starting layer index in parent
            else
                a = 1;
            end
            if s-nl > 10
                b = nl+10;  % ending layer index in parent
            else
                b = s;
            end
            B = S{b};  % Ending layer
            if length(B) > 1
                B = mean(P(Cen(B),:)); % Center of ending layer
            else
                B = P(Cen(B),:);
            end
            A = S{a}; % Staring layer
            if length(A) > 1
                A = mean(P(Cen(A),:));  % Center of starting layer
            else
                A = P(Cen(A),:);
            end
            V = B-A;  % Vector between starting and ending centers
            V = V/norm(V);  % normalize
            Q = vertcat(S{a:b});  % cover sets of the parent between centers
            R = Segs{C(j)};  % cover sets of the child
            R = vertcat(R{:});
            R = [Q; R];  % cover sets of child and parent
            dR = max(distances_to_line(P(Cen(R),:),V,A));  % maximum distance in child and parent
            dQ = max(distances_to_line(P(Cen(Q),:),V,A));  % maximum distance in parent alone
            if (dR-dQ < 0.02) || (dR/dQ < 1.2 && dR-dQ < 0.06)
                K = SChi{C(j)};
                c = length(K);
                if isempty(K)
                    % Remove, no child segments
                    I(C(j)) = false; 
                    Segs{C(j)} = zeros(0,1);
                    SPar(C(j),:) = zeros(1,2);
                    c = setdiff(C,C(j));
                    if size(c,1) < size(c,2)
                        c = c';
                    end
                    SChi{i} = c;
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
                            d = setdiff(C,C(j));
                            if size(d,1) < size(d,2)
                                d = d';
                            end
                            SChi{i} = d;
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
J = SPar(:,1) > 0;
SPar(J,1) = Ind(SPar(J,1));
% Modify SChi
for i = 1:ns
    if I(i)
        C = SChi{i};
        if ~isempty(C)
            C = nonzeros(Ind(C));
            if size(C,1) < size(C,2)
                SChi{i} = C';
            else
                SChi{i} = C;
            end
        else
            SChi{i} = zeros(0,1);
        end
    end
end
SChi = SChi(I);
end % End of function


function [Segs,SPar,SChi] = modify_branches(P,Cen,Bal,Segs,SPar,SChi)

% Make stem and branches as long as possible
ns = size(Segs,1);
nc = ceil(ns/5);
for j = 1:ns
    seg = Segs{j};
    C = SChi{j};
    if size(C,1) < size(C,2)
        C = C';
        SChi{j} = C;
    end
    
    if ~isempty(seg) && ~isempty(C)
        % Find sub-segments
        if j > 1 % Branches
            S = zeros(nc,1); % sub-segments
            S(1) = j;
            t = 2;
            while ~isempty(C)
                n = length(C);
                S(t:t+n-1) = C;
                for i = 1:n
                    c = SChi{C(i)};
                    if isempty(c)
                        SChi{C(i)} = zeros(0,1);
                    end
                end
                C = vertcat(SChi{C});
                t = t+n;
            end
            if t > 2
                t = t-n;
            end
            S = S(1:t);
            
%             R = cell(t,1);
%             for i = 1:t
%                 r = Segs{S(i)};
%                 R{i} = vertcat(r{:}); 
%             end
%             plot_segs(P,R,5,Bal)
%             pause
            
            % Find tip-points
            T = zeros(t,3);
            for i = 1:t
                s = Segs{S(i)};
                p = vertcat(Bal{s{end}});
                if length(p) > 1
                    T(i,:) = mean(P(p,:));
                else
                    T(i,:) = P(p,:);
                end
            end
            
            % Define bottom of the branch
            p = vertcat(Bal{seg{1}});
            if length(p) > 1
                B = mean(P(p,:));
            else
                B = P(p,:);
            end
            
            % End segment is the segment whose tip has greatest distance to
            % the bottom of the branch
            V = mat_vec_subtraction(T,B);
            d = sum(V.*V,2);
            [~,I] = max(d);
            S = S(I(1));
            
        else % Stem
            % For stem the end segment is the highest segment
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
        end
        
        if S ~= j
            % Tip point was not in the current segment, modify segments
            T = zeros(100,1);
            T(1) = S;
            t = 1;
            while S ~= j
                S = SPar(S,1);
                t = t+1;
                T(t) = S;
            end
            T = T(1:t);
            
%             R = cell(t,1);
%             for i = 1:t
%                 r = Segs{T(i)};
%                 R{i} = vertcat(r{:}); 
%             end
%             plot_segs(P,R,3,Bal)
%             pause
            
            % refine branch
            for i = 1:t-2
                I = T(t-i); % segment to be combined to the first segment
                J = T(t-i-1); % above segment's child to be combined next
                SP = SPar(I,2);  % layer index of the child in the parent
                SegP = Segs{j};
                SegC = Segs{I};
                N = size(SegP,1);
                sp = SPar(J,2);  % layer index of the child's child in the child
                if SP >= N-5 % Use the whole parent
                    Segs{j} = [SegP; SegC(1:sp)];
                    if sp < size(SegC,1) % use only part of the child segment
                        Segs{I} = SegC(sp+1:end);
                        SPar(I,2) = N+sp;
                        
                        C = SChi{I};
                        if size(C,1) < size(C,2)
                            C = C';
                        end
                        K = SPar(C,2) <= sp;
                        c = C(~K);
                        SChi{I} = c;
                        SPar(c,2) = SPar(c,2)-sp;
                        C = C(K);
                        SChi{j} = [SChi{j}; C];
                        SPar(C,1) = j;
                        SPar(C,2) = N+SPar(C,2);
                        
                    else % use the whole child segment
                        Segs{I} = cell(0,1);
                        SPar(I,1) = 0;
                        
                        C = SChi{I};
                        if size(C,1) < size(C,2)
                            C = C';
                        end
                        SChi{I} = zeros(0,1);
                        c = setdiff(SChi{j},I);
                        if size(c,2) > size(c,1)
                            c = c';
                        end
                        SChi{j} = [c; C];
                        SPar(C,1) = j;
                        SPar(C,2) = N+SPar(C,2);
                        
                    end
                    
                    T(t-i) = j;
                else % divide the parent segment into two parts
                    ns = ns+1;
                    Segs{ns} = SegP(SP+1:end); % the top part of the parent forms a new segment
                    SPar(ns,1) = j;
                    SPar(ns,2) = SP;
                    
                    Segs{j} = [SegP(1:SP); SegC(1:sp)];
                    
                    C = SChi{j};
                    if size(C,1) < size(C,2)
                        C = C';
                    end
                    K = SPar(C,2) > SP;
                    SChi{j} = C(~K);
                    C = C(K);
                    SChi{ns} = C;
                    SPar(C,1) = ns;
                    SPar(C,2) = SPar(C,2)-SP;
                    SChi{j} = [SChi{j}; ns];
                    if sp < size(SegC,1) % use only part of the child segment
                        Segs{I} = SegC(sp+1:end);
                        SPar(I,2) = SP+sp;
                        
                        C = SChi{I};
                        if size(C,1) < size(C,2)
                            C = C';
                        end
                        K = SPar(C,2) <= sp;
                        SChi{I} = C(~K);
                        SPar(C(~K),2) = SPar(C(~K),2)-sp;
                        C = C(K);
                        SChi{j} = [SChi{j}; C];
                        SPar(C,1) = j;
                        SPar(C,2) = SP+SPar(C,2);
                        
                    else % use the whole child segment
                        Segs{I} = cell(0);
                        SPar(I,1) = 0;
                        
                        C = SChi{I};
                        if size(C,1) < size(C,2)
                            C = C';
                        end
                        c = setdiff(SChi{j},I);
                        if size(c,2) > size(c,1)
                            c = c';
                        end
                        SChi{j} = [c; C];
                        SPar(C,1) = j;
                        SPar(C,2) = SP+SPar(C,2);
                        
                    end
                    T(t-i) = j;
                end
                
            end
            
            % Combine the last segment to the branch
            I = T(1);
            SP = SPar(I,2);
            SegP = Segs{j};
            SegC = Segs{I};
            N = size(SegP,1);
            if SP >= N-5 % Use the whole parent
                Segs{j} = [SegP; SegC];
                Segs{I} = cell(0);
                SPar(I,1) = 0;
                
                C = SChi{I};
                if size(C,1) < size(C,2)
                    C = C';
                end
                c = setdiff(SChi{j},I);
                if size(c,2) > size(c,1)
                    c = c';
                end
                SChi{j} = [c; C];
                SPar(C,1) = j;
                SPar(C,2) = N+SPar(C,2);
                
            else % divide the parent segment into two parts
                ns = ns+1;
                Segs{ns} = SegP(SP+1:end);
                SPar(ns,:) = [j SP];
                Segs{j} = [SegP(1:SP); SegC];
                Segs{I} = cell(0);
                SPar(I,1) = 0;
                
                C = SChi{j};
                if size(C,1) < size(C,2)
                    C = C';
                end
                K = SPar(C,2) > SP;
                SChi{j} = [C(~K); ns];
                C = C(K);
                SChi{ns} = C;
                SPar(C,1) = ns;
                SPar(C,2) = SPar(C,2)-SP;
                
                C = SChi{I};
                c = setdiff(SChi{j},I);
                if size(c,2) > size(c,1)
                    c = c';
                end
                SChi{j} = [c; C];
                SPar(C,1) = j;
                SPar(C,2) = SP+SPar(C,2);
                
            end
            
%             T
%             T = unique(T);
%             t = length(T);
%             R = cell(t,1);
%             I = true(t,1);
%             for i = 1:t
%                 r = Segs{T(i)};
%                 R{i} = vertcat(r{:});
%                 if isempty(r)
%                     I(t) = false;
%                 end
%             end
%             I
%             plot_segs(P,R(I),4,Bal)
%             pause
            
        end
    end
    
end

% Modify indexes by removing empty segments
E = true(ns,1);
for i = 1:ns
    if isempty(Segs{i})
        E(i) = false;
    end
end
Segs = Segs(E);
Ind = (1:1:ns)';
n = nnz(E);
I = (1:1:n)';
Ind(E) = I;
SPar = SPar(E,:);
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
end % End of function


function [SegP,Base] = modify_segment(P,Bal,Cen,SegP,SegC,nl,dmin)

% Expands the base of the branch backwards into its parent segment and
% then removes the expansion from the parent segment.

% Define the directions of the segments
DComp = segment_direction(P,Cen,SegC,1);
DS = segment_direction(P,Cen,SegP,nl);

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


% Determine the number of cover set layers "n" to be checked
A = abs(DComp'*DS);  % abs of cosine of the angle between component and segment directions
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
        base{i+2} = S(I);
        i = i+1;
    end
end
Base = vertcat(base{:});
end % End of function


function D = segment_direction(P,Cen,Seg,nl)

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
    T = Seg{t};
    if length(T) > 1
        Top = mean(P(Cen(T),:));
    else
        Top = P(Cen(T),:);
    end
    V = Top-Bot;
    D = V'/norm(V);
else
    D = zeros(3,1);
end

if size(D,2) == 3
    D = D';
end
 
end % End of function
