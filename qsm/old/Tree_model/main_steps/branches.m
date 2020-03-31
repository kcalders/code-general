function [BoC,BOrd,BPar,BChi,CiB,BVol,BLen,BAng,FCB,BSeg] = branches(Segs,SChi,CiS,CChi,CPar,Rad,Len,Axe,Added)

% ---------------------------------------------------------------------
% BRANCHES.M        Determines the branching structure.
%
% Version 1.0
% Author        Pasi Raumonen
% Created       14 June 2013
%
% This software cannot be distributed without the permission of P. Raumonen
% and it is only for non-commercial use. These restrictions holds also for 
% the modules or subprograms of the software.
% ---------------------------------------------------------------------

% Determines the branches (which cylinders define a branch), their order
% and topological parent-child-relation. Branch number one is the trunk and
% its order is zero. Notice that branch number does not tell its age in the 
% sense that branch number two would be the oldest branch and the number 
% three the second oldest. However, the branch numbers for first-order
% branches are smaller than the numbers for second-order branches, whose
% numbers are smaller than the numbers for third-order-branches, etc.
%
% Inputs:
% Segs      Segments, cell-array
% SChi      Child segments, cell-array
% CiS       Cylinders in the segments, cell-array
% CChi      Child cylinder, cell-raay
% CPar      Parent cylinder, vector
% Rad       Radii of the cylinders
% Len       Lengths of the cylinders
% Axe       Axes of the cylinders
% Added     Logical vector indicating if the cylinder is added to fill gaps    
%
% Outputs:
% BoC       Branch of the cylinder, vector
% BOrd      Branch order, vector
% BPar      Parent branch, vector
% BChi      Child branches, cell-array
% CiB       Cylinders in the branches, cell-array
% BVol      Volumes of the branches
% BLen      Lengths of the branches
% BAng      Branching angles of the branches
% FCB       First cylinders of the branches
% BSeg      Segments forming the branches


%% Branches
nc = size(CChi,1);  % number of cylinder
ns = size(Segs,1);  % number of segments
BOrd = zeros(ns,1);
CiB = cell(ns,1);
BoC = zeros(nc,1);
C = CiS{1};
CiB{1} = C;
BoC(C) = 1;
BVol = zeros(ns,1);
BLen = zeros(ns,1);
BAng = zeros(ns,1);
BVol(1) = pi*sum(Len(C).*Rad(C).^2);
BLen(1) = sum(Len(C));
FCB = zeros(ns,1);
FCB(1) = 1;
BSeg = zeros(ns,1);
BSeg(1) = 1;

S = SChi{1};    % segments under inspection
b = 1;          % branches determined so far
BO = 0;         % branch order under inspection
while ~isempty(S)
    BO = BO+1;
    n = length(S);
    for j = 1:n
        C = CiS{S(j)};
        if ~isempty(C)
            b = b+1;
            CiB{b} = C;
            BOrd(b) = BO;
            BoC(C) = b;
            BVol(b) = pi*sum(Len(C).*Rad(C).^2);
            BLen(b) = sum(Len(C));
            I = Added(C(1));
            if I
                FC = C(2);          % first cyl in the branch
                PC = CPar(C(I));    % parent cylinder of the branch
            else
                FC = C(1);
                PC = CPar(FC);
            end
            FCB(b) = FC;
            if PC > 0
                BAng(b) = 180/pi*acos(Axe(FC,:)*Axe(PC,:)');
            end
            BSeg(b) = S(j);
        end
    end
    S = vertcat(SChi{S});
end
CiB = CiB(1:b);
BOrd = BOrd(1:b);
BVol = 1000*BVol(1:b);  % volumes in liters
BLen = BLen(1:b);
BAng = BAng(1:b);
FCB = FCB(1:b);
BSeg = BSeg(1:b);


%% Branching structure (topology, parent-child-relation)
BChi = cell(b,1);
BPar = zeros(b,1);
for i = 1:b
    C = CiB{i};
    ChildCyls = unique(vertcat(CChi{C}));
    CB = unique(BoC(ChildCyls));  % Child branches
    BChi{i} = CB;
    BPar(CB) = i; 
end