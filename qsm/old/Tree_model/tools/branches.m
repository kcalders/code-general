function [Branches,BranchOrder,BoC] = branches(Segs,SChi,CiS,CPar,CChi,Rad,Len,Axe,Sta,Added)

% Determines the branches (the cylinders forming the branches) and lot of
% attributes of the branches such as the volume and length and their
% distributon by the branch order.


%% Branches, their order, volume, length and angle
nc = size(Rad,1);
ns = size(Segs,1);
Branches = cell(ns,1);
BranchOrder = zeros(ns,1);
BranchVol = zeros(ns,1);
BranchLen = zeros(ns,1);
BranchAng = zeros(ns,1);
BoC = zeros(nc,1);
C = CiS{1};
Branches{1} = C;
BranchVol(1) = pi*sum(Len(C).*Rad(C).^2);
BranchLen(1) = sum(Len(C));
BoC(C) = 1;
FirstCyl = zeros(ns,1);
FirstCyl(1) = 1;
S = SChi{1};  % segment under inspection
b = 1;      % branches determined
BO = 0;     % branch order
while ~isempty(S)
    BO = BO+1;
    n = length(S);
    for j = 1:n
        C = CiS{S(j)};
        if ~isempty(C)
            b = b+1;
            Branches{b} = C;
            BranchOrder(b) = BO;
            BranchVol(b) = pi*sum(Len(C).*Rad(C).^2);
            BranchLen(b) = sum(Len(C));
            BoC(C) = b;
            I = Added(C(1));
            if I
                FC = C(2);          % first cyl in the branch
                PC = CPar(C(I));    % parent cylinder of the branch
            else
                FC = C(1);
                PC = CPar(FC);
            end
            FirstCyl(b) = FC;
            if PC > 0
                BranchAng(b) = 180/pi*acos(Axe(FC,:)*Axe(PC,:)');
            end
        end
    end
    S = unique(vertcat(SChi{S}));
end
Branches = Branches(1:b);
BranchOrder = BranchOrder(1:b);
BranchVol = BranchVol(1:b);
BranchLen = BranchLen(1:b);
BranchAng = BranchAng(1:b);
FirstCyl = FirstCyl(1:b);

BO = max(BranchOrder);

%% Branching structure (topology, parent-child-relation)
ChildBranches = cell(b,1);
ParentBranch = zeros(b,1);
for i = 1:b
    C = Branches{i};
    ChildCyls = unique(vertcat(CChi{C}));
    CB = unique(BoC(ChildCyls));  % Child branches
    ChildBranches{i} = CB;
    ParentBranch(CB) = i; 
end


%% Order distributions
VolByOrder = zeros(BO+1,1);
LenByOrder = zeros(BO+1,1);
AngByOrder = zeros(BO+1,4);
FreByOrder = zeros(BO+1,1);
for i = 0:BO
    I = BranchOrder == i;
    VolByOrder(i+1) = sum(BranchVol(I));
    LenByOrder(i+1) = sum(BranchLen(I));
    if i > 0
        A = BranchAng(I);
        J = A > 0;
        A = A(J);
        if ~isempty(A)
            AngByOrder(i+1,:) = [min(A) mean(A) max(A) std(A)];
        end
    end
    FreByOrder(i+1) = nnz(I)/b;
end
% Normalize volume and length
VolByOrder = VolByOrder/sum(BranchVol);
LenByOrder = LenByOrder/sum(BranchLen);


%% Height distributions
% Branch volume, length, frequency as functions of branch base height
G = min(Sta(:,3));
H = ceil(max(Sta(:,3))-G);
if H >= 5
    VolByBranHeight = zeros(H,1);
    LenByBranHeight = zeros(H,1);
    FreByBranHeight = zeros(H,1);
    for i = 1:H
        I = Sta(FirstCyl,3)-G < i;
        J = Sta(FirstCyl,3)-G >= i-1;
        I = I&J;
        if i == 1
            I(i) = false;
        end
        VolByBranHeight(i) = sum(BranchVol(I));
        LenByBranHeight(i) = sum(BranchLen(I));
        FreByBranHeight(i) = nnz(I)/(b-1);
    end
    h1 = 1:1:H;
else
    VolByBranHeight = zeros(2*H,1);
    LenByBranHeight = zeros(2*H,1);
    FreByBranHeight = zeros(2*H,1);
    for i = 1:2*H
        I = Sta(FirstCyl,3)-G < 0.5*i;
        J = Sta(FirstCyl,3)-G >= 0.5*(i-1);
        I = I&J;
        if i == 1
            I(i) = false;
        end
        VolByBranHeight(i) = sum(BranchVol(I));
        LenByBranHeight(i) = sum(BranchLen(I));
        FreByBranHeight(i) = nnz(I)/(b-1);
    end
    h1 = 0.5:0.5:H;
end
% Normalize volume and length
VolByBranHeight = VolByBranHeight/sum(BranchVol);
LenByBranHeight = LenByBranHeight/sum(BranchLen);

% Volume, length, frequency of cylinders as functions of height
G = min(Sta(:,3));
H = ceil(max(Sta(:,3))-G);
if H >= 5
    VolByHeight = zeros(H,1);
    LenByHeight = zeros(H,1);
    FreByHeight = zeros(H,1);
    for i = 1:H
        I = Sta(:,3)-G < i;
        J = Sta(:,3)-G >= i-1;
        I = I&J;
        VolByHeight(i) = sum(Len(I).*Rad(I).^2);
        LenByHeight(i) = sum(Len(I));
        FreByHeight(i) = nnz(I)/length(I);
    end
    h2 = 1:1:H;
else
    VolByHeight = zeros(2*H,1);
    LenByHeight = zeros(2*H,1);
    FreByHeight = zeros(2*H,1);
    for i = 1:2*H
        I = Sta(FirstCyl,3)-G < 0.5*i;
        J = Sta(FirstCyl,3)-G >= 0.5*(i-1);
        I = I&J;
        if i == 1
            I(i) = false;
        end
        VolByHeight(i) = sum(Len(I).*Rad(I).^2);
        LenByHeight(i) = sum(Len(I));
        FreByHeight(i) = nnz(I)/length(I);
    end
    h2 = 0.5:0.5:H;
end
% Normalize volume and length
VolByHeight = VolByHeight/sum(Len.*Rad.^2);
LenByHeight = LenByHeight/sum(Len);


%% Plot the order and height distributions
figure(7)
bar(0:1:BO,[VolByOrder LenByOrder FreByOrder])
t = title('Branch distributions by order');
legend('Volume','Length','Frequency')
set(t, 'FontSize', 16);
set(t,'FontWeight','bold');
x = xlabel('Order');
set(x, 'FontSize', 16);
set(x,'FontWeight','bold');
y = ylabel('Relative value');
set(y, 'FontSize', 16);
set(y,'FontWeight','bold');
set(gca,'fontsize',16);
set(gca,'FontWeight','bold');

figure(8)
bar(1:1:BO,AngByOrder(2:end,:))
t = title('Branching angle distribution by order');
legend('min','mean','max','std')
set(t, 'FontSize', 16);
set(t,'FontWeight','bold');
x = xlabel('Order');
set(x, 'FontSize', 16);
set(x,'FontWeight','bold');
y = ylabel('Angle in degrees');
set(y, 'FontSize', 16);
set(y,'FontWeight','bold');
set(gca,'fontsize',16);
set(gca,'FontWeight','bold');

figure(9)
bar(h1,[VolByBranHeight LenByBranHeight FreByBranHeight])
t = title('Branch distributions by height');
legend('Volume','Length','Frequency')
set(t, 'FontSize', 16);
set(t,'FontWeight','bold');
x = xlabel('Height (m)');
set(x, 'FontSize', 16);
set(x,'FontWeight','bold');
y = ylabel('Relative value');
set(y, 'FontSize', 16);
set(y,'FontWeight','bold');
set(gca,'fontsize',16);
set(gca,'FontWeight','bold');

figure(10)
bar(h2,[VolByHeight LenByHeight FreByHeight])
t = title('Cylinder distributions by height');
legend('Volume','Length','Frequency')
set(t, 'FontSize', 16);
set(t,'FontWeight','bold');
x = xlabel('Height (m)');
set(x, 'FontSize', 16);
set(x,'FontWeight','bold');
y = ylabel('Relative value');
set(y, 'FontSize', 16);
set(y,'FontWeight','bold');
set(gca,'fontsize',16);
set(gca,'FontWeight','bold');



% I = true(b,1);
% I(1) = false;
% for i = 1:b
%     if BranchOrder(i) > 10
%         C = ChildBranches{i};
%         if isempty(C)
%             I(i) = false;
%         else
%             C = vertcat(ChildBranches{C});
%             if isempty(C)
%                 I(i) = false;
%             end
%         end
%     end
% end
I = BranchOrder <= 10;


disp([round(1000*sum(BranchVol)) round(1000*sum(BranchVol(I))) round(sum(BranchLen)) round(sum(BranchLen(I)))])
disp(max(BranchOrder(I)))
C = vertcat(Branches{I});
plot_cylinder_model(Rad(C),Len(C),Axe(C,:),Sta(C,:),2,30,0.2)
%plot_cylinder_model(Rad,Len,Axe,Sta,1,30,0.2)