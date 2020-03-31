function [Segs,SPar,SChi,SDir] = segments(P,Bal,Nei,Cen,Base,Forb,dmin)

% ---------------------------------------------------------------------
% SEGMENTS.M        Segments the covered point cloud into branches.
%
% Version 1.2
% Author        Pasi Raumonen
% Created       2 December 2013
%
% This software cannot be distributed without the permission of P. Raumonen
% and it is only for non-commercial use. These restrictions holds also for
% the modules or subprograms of the software:
% DEFINE_CUT, CUT_COMPONENTS, DETERMINE_REGION_SIZES, DEFINE_STUDY, 
% CHECK_STUDY_COMPONENT, UPDATE_STUDY, SEGMENT_DIRECTION, 
% COMPONENT_CLASSIFICATION, COMPONENT_DIRECTION, DIRECTION_AT_BASE
% ---------------------------------------------------------------------

% Segments the tree into branches and records their parent-child-relations. 
% Bifurcations are recognized by studying connectivity of a "study"
% region moving along the tree. In case of multiple connected components 
% in "study", the components are checked in multistep process if they are 
% real branches.
%
% Inputs:
% P             Point cloud
% Bal           Cover sets
% Nei           Neigboring cover sets
% Cen           Center points of the cover sets
% Dir           Direction lines of the cover sets
% Prin          Principal components of the cover sets
% Base          Base of the tree
% Forb          Cover sets not part of the tree
% dmin          Minimum diameter of the cover sets
%
% Outputs:
% Segs          Segments found, (n_seg x 1)-cell, each cell contains a cell array the cover sets
% SPar          Parent segment of each segment, (n_seg x 1)-vector,
%                   equals to zero if no parent segment
% SChi          Children segments of each segment, (n_seg x 1)-cell
% SDir          Direction of each segment at its base


nb = size(Bal,1);       % The number of cover sets
a = nb;                 % Estimate for maximum number of segments
SBas = cell(a,1);       % The segment bases found
Segs = cell(a,1);       % The segments found
SPar = zeros(a,2);      % The parent segment of each segment
SChi = cell(a,1);       % The children segments of each segment
SDir = zeros(a,3);      % The direction of each segment at its base

C = zeros(200,1);
for i = 1:a
    SChi{i} = C;
end
NChi = zeros(a,1);

Fal = false(nb,1);      % Logical false-vector for cover sets

s = 1;                  % The index of the segment under expansion
b = s;                  % The index of the latest found base

SBas{s} = Base;
Seg = cell(200000,1);    % The cover set layers in the current segment
Seg{1} = Base;
CutPre = Base;       % First cut to determine the study region

ForbAll = Fal;       % The forbidden sets
ForbAll(Forb) = true;
ForbAll(Base) = true;
Forb = ForbAll;      % The forbidden sets for the segment under expansion

Continue = true; % True as long as the component can be segmented further 
NewSeg = true;   % True if the first Cut for the current segment
Trunk = true;    % True for trunk segment
nl = 1;          % The number of cover set layers currently in the segment

% Segmenting stops when there are no more segments to be found
% or when the maximum number of segments is found
while Continue && (b < nb)
    
    % Update the forbidden sets
    Forb(CutPre) = true;
    
    % Define the study
    if NewSeg
        Cut = define_cut(Nei,CutPre,Forb);
        CutSize = length(Cut);
        if CutSize > 0
            
            [CutComps,CompSize] = cut_components(Nei,Cut,CutSize,Fal);
            
            if length(CompSize) == 1
                nc = 1;
            else
                if NewSeg
                    NewSeg = false;
                    ns = determine_region_sizes(P,Cen,Cut,dmin,Trunk);
                end
                
                [StudyComps,Cont] = define_study(Nei,CutComps,Forb,ns,Fal);
                nc = length(Cont);

            end
        else
            nc = 0;
        end
    else
        
        [StudyComps,Cont,Cut] = update_study(Nei,StudyComps,Cont,Forb,Fal,ns);
        nc = length(Cont);
        CutSize = length(Cut);
        
    end
    
    % Classify study components
    if nc == 1
        CutPre = Cut;
        nl = nl+1;
        if size(Cut,2) > 1
            Seg{nl} = Cut';
        else
            Seg{nl} = Cut;
        end
    elseif nc > 1
        % Determine the direction of the segment below "Study" region.
        Thick = CutSize > 10;
        [DSeg,CSeg,RSeg] = segment_direction(P,Bal,Cen,Nei,Segs,Seg,...
            2*ns,nl,Trunk,Thick,SPar(s,1));
        
        % Classify the components of the Study region
        [Class,DComp] = component_classification(P,Bal,Nei,Cen,...
            StudyComps,Cont,DSeg,CSeg,RSeg,Seg,nl,nc,CutSize);
        
        I = true(nc,1);  
        for i = 1:nc
            if Class(i) == 1
                I(i) = false;
                C = StudyComps{i};
                Base = C{1};
                Comp = vertcat(C{:});
                
                ForbAll(Base) = true;
                Forb(Comp) = true;
                Forb(Base) = true;
                
                J = Forb(Cut);
                Cut = Cut(~J);
                
                b = b+1;
                SBas{b} = Base;
                SPar(b,:) = [s nl];
                NChi(s) = NChi(s)+1;
                SChi{s}(NChi(s)) = b;
                
            elseif Class(i) == 2
                % Remove small number of points
                I(i) = false;
                Comp = StudyComps{i};
                Comp = vertcat(Comp{:});
                
                ForbAll(Comp) = true;
                Forb(Comp) = true;
                
                J = Forb(Cut);
                Cut = Cut(~J);
            end
        end
        
        StudyComps = StudyComps(I);
        Cont = Cont(I);
        
        % Define a new cut.
        % If cut is empty, determine the new base
        if isempty(Cut)
            Segs{s} = Seg(1:nl);
            SDir(s,:) = direction_at_base(P,Bal,Cen,Seg,nl,ns);
            S = vertcat(Seg{1:nl});
            ForbAll(S) = true;

            if s < b;
                
                s = s+1;
                CutPre = SBas{s};
                Seg = cell(10000,1);
                Seg{1} = CutPre;
                Forb = ForbAll;
                NewSeg = true;
                Trunk = false;
                nl = 1;
            else
                Continue = false;
            end
        else
            if size(Cut,2) > 1
                Cut = Cut';
            end
            CutPre = Cut;
            nl = nl+1;
            Seg{nl} = Cut;
        end
    
    else
        % If the study region has zero size, then the current segment is
        % complete and determine a new base
        Segs{s} = Seg(1:nl);
        SDir(s,:) = direction_at_base(P,Bal,Cen,Seg,nl,ns);
        S = vertcat(Seg{1:nl});
        ForbAll(S) = true;
        
        if s < b;
            s = s+1;
            CutPre = SBas{s};
            Seg = cell(1000,1);
            Seg{1} = CutPre;
            Forb = ForbAll;
            NewSeg = true;
            Trunk = false;
            nl = 1;
        else
            Continue = false;
        end
    end
end
Segs = Segs(1:b);
SPar = SPar(1:b,:);
schi = SChi(1:b);
SDir = SDir(1:b,:);

SChi = cell(b,1);
for i = 1:b
    if NChi(i) > 0
        SChi{i} = schi{i}(1:NChi(i));
    end
end

% Display the number of segments found
str = ['    ',num2str(b),' segments found'];
disp(str)


end % End of the main function


% Define subfunctions

function Cut = define_cut(Nei,CutPre,Forb)

% Defines the "Cut" region
Cut = unique(vertcat(Nei{CutPre}));
I = Forb(Cut);
Cut = Cut(~I);
end % End of function 


function [Components,CompSize] = cut_components(Nei,Cut,CutSize,Fal)

% Define the connected components of the Cut
if CutSize == 1
    % Cut is connected and therefore Study is also
    CompSize = 1;
    Components = Cut;
elseif CutSize == 2
    I = Nei{Cut(1)} == Cut(2);
    if any(I)
        Components = Cut;
        CompSize = 1;
    else
        Components = cell(2,1);
        Components{1} = Cut(1);
        Components{2} = Cut(2);
        CompSize = [1 1];
    end
elseif CutSize == 3
    I = Nei{Cut(1)} == Cut(2);
    J = Nei{Cut(1)} == Cut(3);
    K = Nei{Cut(2)} == Cut(3);
    if any(I)+any(J)+any(K) >= 2
        CompSize = 1;
        Components = Cut;
    elseif any(I)
        Components = cell(2,1);
        Components{1} = Cut(1:2);
        Components{2} = Cut(3);
        CompSize = [2 1];
    elseif any(J)
        Components = cell(2,1);
        Components{1} = Cut([1 3]');
        Components{2} = Cut(2);
        CompSize = [2 1];
    elseif any(K)
        Components = cell(2,1);
        Components{1} = Cut(2:3);
        Components{2} = Cut(1);
        CompSize = [2 1];
    else
        CompSize = [1 1 1];
        Components = cell(3,1);
        Components{1} = Cut(1);
        Components{2} = Cut(2);
        Components{3} = Cut(3);
    end
else
    Components = cell(CutSize,1);
    CompSize = zeros(CutSize,1);
    Fal(Cut) = true;
    nc = 0;      % number of components found
    m = Cut(1);
    i = 0;
    while i < CutSize
        Added = Nei{m};
        I = Fal(Added);
        Added = Added(I);
        a = length(Added);
        Comp = zeros(CutSize,1);
        Comp(1) = m;
        Fal(m) = false;
        t = 1;
        while a > 0
            Comp(t+1:t+a) = Added;
            Fal(Added) = false;
            t = t+a;
            Ext = unique(vertcat(Nei{Added}));
            I = Fal(Ext);
            Added = Ext(I);
            a = length(Added);
        end
        i = i+t;
        nc = nc+1;
        Components{nc} = Comp(1:t);
        CompSize(nc) = t;
        if i < CutSize
            I = Fal(Cut);
            m = Cut(I);
            m = m(1);
        end
    end
    Components = Components(1:nc);
    CompSize = CompSize(1:nc);
end

end % End of function


function ns = determine_region_sizes(P,Cen,Cut,dmin,Trunk)

% Determines the size of "Study" region as the number of cover set layers. 
% The size is such that "Study" region is about as long as its diameter 
% but at least two cover set layers
if length(Cut) > 1
    D = dimensions(P(Cen(Cut),:)); % approximate diameter of the segment
else
    D = 2*dmin;
end
if Trunk
    ns = max(ceil(0.2/dmin),ceil(D(1)/dmin));
else
    ns = max([4 ceil(0.015/dmin) 2*ceil(D(1)/dmin)]);
    ns = min(ns,ceil(0.2/dmin));
end

end % End of function


function [Comps,Cont] = define_study(Nei,CutComps,Forb,ns,Fal)

% Define the study components from given cut components

% Initialize the study components
n = size(CutComps,1);
Comps = cell(n,1);
for i = 1:n
    C = cell(ns,1);
    C{1} = CutComps{i};
    Comps{i} = C;
end

% Expand the components one layer at a time anc check if the components
% overlap or are connected, in which case combine those components
i = 2; % layer under construction
Cont = true(n,1);
Bpre = CutComps;
Forb(vertcat(Bpre{:})) = true;
while i <= ns && any(Cont)
    
    % First expand components by one layer
    B = cell(n,1); % the layers of each component
    for j = 1:n
        if Cont(j)
            b = unique(vertcat(Nei{Bpre{j}}));
            I = Forb(b);
            b = b(~I);
            if isempty(b)
                Cont(j) = false;
                C = Comps{j};
                Comps{j} = C(1:i-1);
            elseif size(b,2) > 1
                B{j} = b';
            else
                B{j} = b;
            end
        end
    end
    
    % Check if the layers overlap or are connected (are neighbors), and
    % combine those overlaping or connected components
    I = true(n,1);
    for j = 1:n-1
        if Cont(j) && I(j)
            for k = j+1:n
                if Cont(k) && I(k)
                    Fal(B{j}) = true;
                    J = Fal(B{k});
                    K = Fal(vertcat(Nei{B{k}}));
                    Fal(B{j}) = false;
                    if any(J) || any(K)
                        b = union(B{j},B{k});
                        if size(b,2) > 1
                            B{j} = b';
                        else
                            B{j} = b;
                        end
                        I(k) = false;
                        C = Comps{j};
                        c = Comps{k};
                        for l = 1:i-1
                            L = union(C{l},c{l});
                            if size(L,2) > 1
                                C{l} = L';
                            else
                                C{l} = L;
                            end
                        end
                        Comps{j} = C;
                    end
                end
            end
        end
    end
    
    % Remove components combined to other components
    if ~all(I)
        B = B(I);
        Comps = Comps(I);
        Cont = Cont(I);
        n = nnz(I);
    end
    
    % Update the i:th layer of components
    for j = 1:n
        if Cont(j)
            C = Comps{j};
            C{i} = B{j};
            Comps{j} = C;
        end
    end
    Bpre = B;
    Forb(vertcat(Bpre{:})) = true;
    i = i+1;
end

if size(Cont,1) < size(Cont,2)
    Cont = Cont';
end

end % End of function


function [Comps,Cont] = check_study_component(Nei,StudyComp,CutComps,Cont0,Forb,ns,Fal)

% If the first layer of a component has multiple components (CutComps),
% check if the component itself has multiple components and define those
% components

% Initialize the possible new components
n = size(CutComps,1);
Comps = cell(n,1);
for i = 1:n
    C = cell(ns,1);
    C{1} = CutComps{i};
    Comps{i} = C;
end

% Expand the components one layer at a time anc check if the components
% overlap or are connected, in which case combine those components
i = 2; % layer under construction
Cont = true(n,1);
Bpre = CutComps;
Forb(vertcat(Bpre{:})) = true;
while i <= ns && any(Cont)
    
    % First expand components by one layer
    B = cell(n,1); % the layers of each component
    for j = 1:n
        if Cont(j) % expand only those components that continue
            b = unique(vertcat(Nei{Bpre{j}}));
            I = Forb(b);
            b = b(~I);
            if isempty(b)
                Cont(j) = false;
                C = Comps{j};
                Comps{j} = C(1:i-1);
            elseif size(b,2) > 1
                B{j} = b';
            else
                B{j} = b;
            end
        end
    end
    
    % Check if the layers overlap or are connected (are neighbors), and
    % combine those overlaping or connected components
    I = true(n,1);
    for j = 1:n-1
        if Cont(j) && I(j)
            for k = j+1:n
                if Cont(k) && I(k)
                    Fal(B{j}) = true;
                    J = Fal(B{k});
                    K = Fal(vertcat(Nei{B{k}}));
                    Fal(B{j}) = false;
                    if any(J) || any(K)
                        b = union(B{j},B{k});
                        if size(b,2) > 1
                            B{j} = b';
                        else
                            B{j} = b;
                        end
                        I(k) = false;
                        C = Comps{j};
                        c = Comps{k};
                        for l = 1:i
                            L = union(C{l},c{l});
                            if size(L,2) > 1
                                C{l} = L';
                            else
                                C{l} = L;
                            end
                        end
                        Comps{j} = C;
                    end
                end
            end
        end
    end
    
    % Remove components combined to other components
    if ~all(I)
        B = B(I);
        Comps = Comps(I);
        Cont = Cont(I);
        n = nnz(I);
    end
    
    % Update the i:th layer of components
    for j = 1:n
        if Cont(j)
            C = Comps{j};
            C{i} = B{j};
            Comps{j} = C;
        end
    end
    
    % If more than one component still left, continue expansion
    if n > 1
        i = i+1;
        Bpre = B;
        Forb(vertcat(Bpre{:})) = true;
    else
        i = ns+1;
        Comps = StudyComp;
        Cont = Cont0;
    end
end

if size(Cont,1) < size(Cont,2)
    Cont = Cont';
end

end % End of function


function [StudyComps,Cont,Cut] = update_study(Nei,StudyComps,Cont,Forb,Fal,ns)

% Update the components of the Study

forb = Forb;

% Update the Forb for study (components) expansion
n = length(Cont);
for i = 1:n
    S = StudyComps{i};
    Forb(vertcat(S{:})) = true;
end

% Update (expand) the study (components)
I = true(n,1);  % components of non-zero length
for i = 1:n
    C = StudyComps{i};
    if Cont(i)
        b = unique(vertcat(Nei{C{end}}));
        J = Forb(b);
        b = b(~J);
        if isempty(b)
            Cont(i) = false;
            C = C(2:end);
        else
            C(1:end-1) = C(2:end);
            if size(b,2) > 1
                C{end} = b';
            else
                C{end} = b;
            end
        end
        if isempty(C)
            I(i) = false;
        end
        StudyComps{i} = C;
    else
        if size(C,1) > 1
            StudyComps{i} = C(2:end);
        else
            I(i) = false;
        end
    end
end

% Remove zero-length components
if ~all(I)
    StudyComps = StudyComps(I);
    Cont = Cont(I);
    I = I(I);
    n = length(I);
end

% Check if the expansions overlap and combine overlapping components
if n > 1
    for i = 1:n-1
        if Cont(i) && I(i)
            C = StudyComps{i};
            B = C{end};
            for j = i+1:n
                if Cont(j) && I(j)
                    c = StudyComps{j};
                    b = c{end};
                    Fal(B) = true;
                    J = Fal(b);
                    Fal(B) = false;
                    if any(J)
                        B = union(B,b);
                        I(j) = false;
                        for k = 1:ns
                            b = union(C{k},c{k});
                            if size(b,2) > 1
                                C{k} = b';
                            else
                                C{k} = b;
                            end
                        end
                        StudyComps{i} = C;
                    end
                end
            end
        end
    end
end

% Remove components combined to other components
if ~all(I)
    StudyComps = StudyComps(I);
    Cont = Cont(I);
    I = I(I);
    n = length(I);
end

for i = 1:n
    C = StudyComps{i};
    cut = C{1};
    
    % Determine the connected components of cut (first layer of the current component)
    [CutComps,CompSize] = cut_components(Nei,cut,length(cut),Fal);
    
    % If cut has multiple component, check the current component if it has
    % many components, in which case create new components
    if length(CompSize) > 1
        I = true(n,1);
        I(i) = false;
        F = forb;
        for j = 1:n
            if I(j)
                S = StudyComps{j};
                F(vertcat(S{:})) = true;
            end
        end
        
        [Comps,cont] = check_study_component(Nei,C,CutComps,Cont(i),F,ns,Fal);
        m = length(cont);
        
        if m > 1
            StudyComps{i} = Comps{1};
            StudyComps(end+1:end+m-1) = Comps(2:m);
            Cont(i) = cont(1);
            Cont(end+1:end+m-1) = cont(2:m);
        end
    end
end

n = length(Cont);
Cut = cell(n,1);
for i = 1:n
    C = StudyComps{i};
    Cut{i} = C{1}; 
end
Cut = vertcat(Cut{:});

if size(Cont,1) < size(Cont,2)
    Cont = Cont';
end

end % End of function


function [D,C,R] = segment_direction(P,Bal,Cen,Nei,Segs,SL,ss,nl,Trunk,Thick,ParSeg)

% Defines the direction and center of the segment under the study region.

% Direction
if nl-ss > 0
    b = nl-ss+1;
else
    b = 1;
end

if nl > b
    B = SL{b};
    if length(B) > 1
        Bot = mean(P(Cen(B),:));
    else
        Bot = P(Cen(B),:);
    end
    T = SL{nl};
    if length(T) > 1
        Top = mean(P(Cen(T),:));
    else
        Top = P(Cen(T),:);
    end
    V = Top-Bot;
    D = V/norm(V);
elseif ParSeg > 0
    S = Segs{ParSeg};
    S = vertcat(S{:});
    B = SL{1};
    M = unique(vertcat(Nei{B}));
    N = intersect(M,S);
    if ~isempty(N)
        if length(N) > 1
            Bot = mean(P(Cen(N),:));
        else
            Bot = P(Cen(N),:);
        end
        if length(B) > 1
            Top = mean(P(Cen(B),:));
        else
            Top = P(Cen(B),:);
        end
        V = Top-Bot;
        D = V/norm(V);
    else
        D = zeros(1,3);
    end
else
    D = zeros(1,3);
end

% Center
if nl <= 4
    s = vertcat(SL{1:nl});
else
    s = vertcat(SL{nl-4:nl});
end
if length(s) > 1
    C = mean(P(Cen(s),:));
else
    C = P(Cen(s),:);
end

if (Trunk || Thick)
    % Radius
    B = vertcat(SL{b:nl});
    if length(B) > 10
        Points = P(Cen(B),:);
    else
        Points = P(vertcat(Bal{B}),:);
    end
    R0 = median(distances_to_line(Points,D,C));
    if size(Points,1) > 10
        [C0, D0, R] = lscylinder(Points, C', D', R0, 0.1, 0.1);
        %C = C0';
        %D = D0';
    else
        R = R0;
    end
else
    R = 0;
end

if size(D,1) < size(D,2)
    D = D';
end
if size(C,1) < size(C,2)
    C = C';
end

end % End of function


function [Class,DComp] = component_classification(P,Bal,Nei,Cen,...
    Comps,Cont,DSeg,CSeg,RSeg,Seg,nl,nc,CutSize)

% Classifies study region components:
% Class(i) == 0 continuations
% Class(i) == 1 branch
% Class(i) == 2 remove


CS = cell(nc,1);        % Layer size of the components
H = zeros(nc,2);        % Maximum heights for the components
CompSize = zeros(nc,1); 
for i = 1:nc
    Comp = Comps{i};
    NoL = size(Comp,1);  % Number of layers
    L = zeros(NoL,1);
    for j = 1:NoL
        L(j) = length(Comp{j});
    end
    CS{i} = L;
    CompSize(i) = sum(L);
    [h,I] = max(P(Cen(Comp{end}),3));
    H(i,:) = [h I];
end
StudySize = sum(CompSize);

Class = ones(nc,1);     % true if a component is a branch to be further segmented
DComp = zeros(nc,3);    % direction lines of the branch components
Conti = 0;               % index of the continuation component
% Simple initial classification
for i = 1:nc
    Comp = Comps{i};
    comp = vertcat(Comp{:});
    NoL = size(Comp,1);  % Number of layers
    L = CS{i};
    BaseSize = L(1);
    np = length(vertcat(Bal{comp})); % Number of points in the component
    if NoL == 1 && np < 10 && ~Cont(i)
        % component has no expansion, not a branch
        Class(i) = 0;
    elseif BaseSize == 1 && CompSize(i) <= 2 && np < 10 && ~Cont(i)
        % component has very small expansion, not a branch
        Class(i) = 0;
    elseif BaseSize/CutSize < 0.05 && (2*BaseSize >= CompSize(i) || np < 30) && ~Cont(i)
        % component has very small expansion or is very small, not a branch
        Class(i) = 0;
    elseif CompSize(i) <= 3 && np < 5 && ~Cont(i)
        % very small component, not a branch
        Class(i) = 0;
    elseif BaseSize/CutSize >= 0.7 || CompSize(i) >= 0.7*StudySize
        % continuation of the segment
        Class(i) = 0;
        Conti = i;
    else
        % Component is probably a branch
        DComp(i,:) = component_direction(P,Bal,Cen,Comp,CSeg,DSeg);
    end
end
Cosines = DComp*DSeg;  % the cosines of the angles between the
                        % components and the segment

% If continuation does not exists after initial classification, then try to
% define one
if Conti == 0
    % Define possible continuation
    I = Cosines > 0.85;          % nearly parallel components
    J = CompSize > 0.35*StudySize; % the largest component
    K = 1:1:nc;
    L = I&J&Cont;
    L = K(L);
    if length(L) == 1
        Class(L) = 0;
        Conti = K(L);
    elseif length(L) > 1
        [~,I] = max(CompSize(L));
        L = L(I);
        Class(L) = 0;
        Conti = L;
    end
end


I = Class == 1;
if Conti > 0 && any(I)
    if RSeg > 0
        Comp = Comps{Conti};
        a = 2;
        if a > size(Comp,1)
            a = size(Comp,1);
        end
       	Con = Comp(1:a);
        Con = vertcat(Con{:});
        if length(Con) > 5
            Points = P(Cen(Con),:);
        else
            Points = P(vertcat(Bal{Con}),:);
        end
        AP0 = mean(Points)';
        CA0 = component_direction(P,Bal,Cen,Comp,CSeg,DSeg);
        %d = distances_to_line(Points,CA0,AP0);
        %Ds = mean(d);
        for i = 1:nc
            if I(i)
                Comp = Comps{i};
                n = 3;
                if n > size(Comp,1)
                    n = size(Comp,1);
                end
                D = zeros(n,2);
                for j = 1:n
                    C = Comp{j};
                    if length(C) > 5
                        points = P(Cen(C),:);
                    else
                        points = P(vertcat(Bal{C}),:);
                    end
                    d = distances_to_line(points,CA0,AP0);
                    D(j,:) = [min(d) max(d)];
                    if j == 1
                        S = intersect(vertcat(Nei{C}),Seg{nl});
                        if length(S) > 5
                            points = P(Cen(S),:);
                        else
                            points = P(vertcat(Bal{S}),:);
                        end
                        d = distances_to_line(points,CA0,AP0);
                        D(1,1) = mean(d);
                    end
                end
                if any(D(2:n,2) < 1.1*D(1,1))
                    Class(i) = 0;
%                     Ds
%                     disp(D)
%                     C = vertcat(Comp{:});
%                     S = vertcat(Seg{nl-5:nl});
%                     plot_segments(P,Bal,1,Con,C,S)
%                     pause
                end
            end
        end
    end
end


% % Determine continuation for the trunk if it does not yet exist
% if Trunk && Conti == 0
%     D = zeros(nc,3);
%     H = zeros(nc,1);
%     for i = 1:nc
%         C = Comps{i};
%         Cend = C{end};
%         [d,~,h] = distances_to_line(P(Cen(Cend),:),DSeg,CSeg);
%         D(i,2) = min(d);
%         %H(i) = max(h);
%         H(i) = max(P(Cen(Cend),3));
%         C = vertcat(C{:});
%         d = distances_to_line(P(Cen(C),:),DSeg,CSeg);
%         %D(i,1) = max(d);
%         I = d < 1.2*RSeg;
%         D(i,3) = nnz(I)/length(I);
%     end
%     [~,I1] = max(H);
%     [m,I2] = max(D(:,3));
%     [~,I3] = min(D(:,2));
%     [~,I4] = max(Cosines);
%     if I1 == I4 && I1 == I3 && I1 == I2
%         Conti = I1;
%     elseif I1 == I4 && I1 == I2 
%         Conti = I1;
%     elseif m ~= 0  && I1 == I2
%         Conti = I1;
%     elseif m ~= 0  && I2 == I4
%         Conti = I1;
%     elseif I1 == I4
%         Conti = I1;
%     elseif m ~= 0 && Cosines(I2) ~= 0
%         Conti = I2;
%     elseif m == 0 && I1 == I4
%         Conti = I1;
%     elseif m == 0
%         Conti = I4;
%     else
%         Conti = I4;
%     end
%     Class(Conti) = 0;
%     
%     if nl > 5
%         Cosines
%         Conti
%         Class
%         D
%         H
%         plot_segs(P,Comps{1},2,Bal)
%         for i = 2:nc
%             hold on
%             plot_segs(P,Comps{i},2,Bal)
%         end
%         S = vertcat(Seg{nl-5:nl});
%         hold on
%         plot_segments(P,Bal,2,S)
%         hold on
%         quiver3(CSeg(1),CSeg(2),CSeg(3),DSeg(1),DSeg(2),DSeg(3),'.g','Linewidth',10)
%         for i = 1:nc
%             hold on
%             C = Comps{i};
%             C = C{1};
%             M = mean(P(vertcat(Bal{C}),:));
%             quiver3(M(1),M(2),M(3),DComp(i,1),DComp(i,2),DComp(i,3),'r','Linewidth',10)
%         end
%         hold off
%         pause
%     end
% end
   
% if Conti == 0
%     % Define possible continuation
%     I = Cosines > 0.85;          % nearly parallel components
%     J = CompSize > 0.35*StudySize; % the largest component
%     K = 1:1:nc;
%     L = I&J&Cont;
%     % Define possible continuation
%     I = Cosines > 0.9;          % nearly parallel components
%     J = CompSize > 0.5*StudySize; % the largest component
%     K = 1:1:nc;
%     if any(I&J&Cont)
%         % continuation is the component with the largest base
%         % (also quite parallel to the segment)
%         Class(I&J&Cont) = 0;
%         Conti = K(I&J&Cont);
%     elseif any(I&Cont)
%         % continuation is the largest component nearly parallel to the segment
%         [~,J] = max(CompSize(I));
%         K = K(I);
%         Class(K(J)) = 0;
%         Conti = K(J);
%     elseif any(Cont)
%         K = K(Cont);
%         [m,I] = max(Cosines(K));
%         Class(K(I(1))) = 0;
%         Conti = K(I(1));
%     else
%         % continuation is the component most parallel to the segment
%         [m,I] = max(Cosines);
%         if m > 0.8
%             Class(I(1)) = 0;
%             Conti = K(I(1));
% %         else
% %             I = Class == 1;
% %             if nnz(I) == 1
% %                 Conti = K(I);
% %                 Class(I) = 0;
% %             end
%         end
%     end
%end

% if continuation exists and there are branch components nearly parallel to
% the continuation, then define the center, direction and radius of the
% continuation
% if Conti > 0 && any(Class == 1)
%     Cosines = DComp*DComp(Conti,:)';
%     I = Class == 1;
%     J = abs(Cosines) > 0.95;
%     if (Trunk || RSeg > 0) && any(I&J)
%         Comp = Comps{Conti};
%         Comp = vertcat(Comp{:});
%         Points = P(vertcat(Bal{Comp}),:);
%         AP0 = mean(Points)';
%         CA0 = DComp(Conti,:)';
%         R0 = RSeg;
%         d = distances_to_line(Points,CA0,AP0);
%         I = d < 1.2*R0;
%         if nnz(I) > 10
%             Points = Points(I,:);
%         else
%             I = d < 1.5*R0;
%             if nnz(I) > 10
%                 Points = Points(I,:);
%             end
%         end
%         
%         if size(Points,1) > 10
%             [APCon, DCon, RCon] = lscylinder(Points, AP0, CA0, R0, 0.1, 0.1);
%         else
%             APCon = AP0;
%             DCon = CA0;
%             RCon = R0;
%         end
%         h = Points*DCon;
%         hmin = min(h);
%         hpoint = DCon'*APCon;
%         APCon = APCon-(hpoint-hmin)*DCon;
%     end
% end

% for i = 1:nc
%     % Check if a branch nearly parallel to the continuation is part of the
%     % continuation
%     if (Trunk || RSeg > 0) && (Conti > 0) && (Class(i) == 1) && (abs(Cosines(i)) > 0.95)
%         Comp = Comps{i};
%         Comp = vertcat(Comp{:});
%         Points = P(vertcat(Bal{Comp}),:);
%         AP0 = mean(Points)';
%         CA0 = DComp(i,:)';
%         R0 = RCon;
%         d = distances_to_line(Points,CA0,AP0);
%         I = d < 1.2*R0;
%         if nnz(I) > 10
%             Points = Points(I,:);
%         else
%             I = d < 1.5*R0;
%             if nnz(I) > 10
%                 Points = Points(I,:);
%             end
%         end
%         
%         if size(Points,1) > 10
%             [AP, CA, R] = lscylinder(Points, AP0, CA0, R0, 0.1, 0.1);
%         else
%             AP = AP0;
%             CA = CA0;
%             R = R0;
%         end
%         h = Points*CA;
%         hmin = min(h);
%         hpoint = CA'*AP;
%         AP = AP-(hpoint-hmin)*CA;
%         d = distances_to_line(AP',DCon,APCon);
%         if (R > 0.75*RCon) && (d < 0.5*RCon)
%             Class(i) = 0;
%         end
%     end
%         
%     % Check if a branch component is bifurcated at its base, in which case
%     % it is not classified as a branch but a "continuation"
%     if (Class(i) == 1) && ((abs(Cosines(i)) > 0.6) || (Trunk && abs(Cosines(i)) > 0.4))
%         L = CS{i};
%         n = length(L);
%         if (n >= 3) && (L(1)/CutSize > 0.2) && (L(1) > 2)
%             Comp = Comps{i};
%             N = ones(n-2,1);
%             j = 1;
%             next = true;
%             while (j <= n-1) && next
%                 comp = vertcat(Comp{j+1:end});
%                 [~,S] = connected_components(Nei,comp,2);
%                 N(j) = length(S);
%                 if j+2 < n
%                     comp = vertcat(Comp{1:j+2});
%                 else
%                     comp = vertcat(Comp{1:n});
%                 end
%                 A = P(vertcat(Bal{comp}),:);
%                 [D,dir] = dimensions(A);
%                 k = abs(dir(1:3)*DComp(i,:)');
%                 
%                 if (k < 0.8 || D(4) < 0.5) && (sum(N(1:j)) < 2*j+2)
%                     if (k < 0.8 || D(4) < 0.5) && (D(1)/3 < j*dmin) && (sum(N(1:j)) < j+2)
%                         Class(i) = 1;
%                         next = false;
%                     end
%                 else
%                     if sum(N(1:j)) > j+2
%                         Class(i) = 0;
%                     end
%                     next = false;
%                 end
%                 j = j+1;
%             end
%         end
%     end
% end
end % End of function


function D = component_direction(P,Bal,Cen,Comp,CSeg,DSeg)

% Defines the direction of the component

n = size(Comp,1);
if n > 3
    if length(Comp{1}) > 1
        Bot = mean(P(Cen(Comp{1}),:));
    else
        B = vertcat(Bal{Comp{1}});
        if length(B) > 1
            Bot = mean(P(B,:));
        elseif length(B) == 1
            Bot = P(B,:);
        else
            Bot = CSeg';
        end
    end
else
    Bot = CSeg';
end

Top = Comp{n};
if length(Top) > 1
    [d,~,h] = distances_to_line(P(Cen(Top),:),DSeg,CSeg);
    if min(d) < max(d)+0.02
        I = d < min(d)+0.02;
        Top = Top(I);
        if length(Top) > 1
            Top = mean(P(Cen(Top),:));
        elseif length(Top) == 1
            Top = P(Cen(Top),:);
        else
            Top = Comp{n};
            Top = mean(P(Cen(Top),:));
        end
    end
elseif length(Top) == 1
    Top = P(Cen(Top),:);
else
    Top = vertcat(Comp{1:n});
    if length(Top) > 1
        [d,~,h] = distances_to_line(P(Cen(Top),:),DSeg,CSeg);
        if min(d) < max(d)+0.02
            I = d < min(d)+0.02;
            Top = Top(I);
            if length(Top) > 1
                Top = mean(P(Cen(Top),:));
            elseif length(Top) == 1
                Top = P(Cen(Top),:);
            else
                Top = vertcat(Comp{1:n});
                Top = mean(P(Cen(Top),:));
            end
        end
    elseif length(Top) == 1
        Top = P(Cen(Top),:);
    end
end

I = false(8,1);
I(1) = isempty(Top);
I(2) = isempty(Bot);
I(3) = any(isnan(Top));
I(4) = any(isnan(Bot));
I(5) = size(Top,1) ~= 1;
I(6) = size(Top,2) ~= 3;
I(7) = size(Bot,1) ~= 1;
I(8) = size(Bot,2) ~= 3;

if any(I)
    D = zeros(1,3);
else
    V = Top-Bot;
    D = V/norm(V);
end

end % End of function


function D = direction_at_base(P,Bal,Cen,Seg,nl,ns)

% Defines the direction line of the segment at its base
if length(Seg{1}) > 1
    Bot = mean(P(Cen(Seg{1}),:));
else
    B = vertcat(Bal{Seg{1}});
    if length(B) > 1
        Bot = mean(P(B,:));
    else
        Bot = P(B,:);
    end
end
if nl >= 2*ns
    if length(Seg{2*ns}) > 1
        Top = mean(P(Cen(Seg{2*ns}),:));
    else
        B = vertcat(Bal{Seg{2*ns}});
        if length(B) > 1
            Top = mean(P(B,:));
        else
            Top = P(B,:);
        end
    end
elseif nl >= ns
    if length(Seg{ns}) > 1
        Top = mean(P(Cen(Seg{ns}),:));
    else
        B = vertcat(Bal{Seg{ns}});
        if length(B) > 1
            Top = mean(P(B,:));
        else
            Top = P(B,:);
        end
    end
elseif nl >= 2
    if length(Seg{nl}) > 1
        Top = mean(P(Cen(Seg{nl}),:));
    else
        B = vertcat(Bal{Seg{nl}});
        if length(B) > 1
            Top = mean(P(B,:));
        else
            Top = P(B,:);
        end
    end
else
    Top = Bot+randn(1,3);
end

if isempty(Top) || isempty(Bot)
    D = zeros(1,3);
else
    V = Top-Bot;
    D = V/norm(V);
end

end % End of function
