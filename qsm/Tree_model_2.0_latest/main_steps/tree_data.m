function [TreeData,Taper,BODist,BSDist,CODist,CSDist,strX,strS,BO] = ...
    tree_data(Rad,Len,Sta,Axe,BOrd,CiB,BVol,BLen,trunk,string)

% ---------------------------------------------------------------------
% TREE_DATA.M       Calculates some tree attributes. 
%
% Version 1.3
% Author        Pasi Raumonen
% Created       16 March 2014
%
% This software cannot be distributed without the permission of P. Raumonen
% and it is only for non-commercial use. 
% ---------------------------------------------------------------------

% Inputs:
% Rad       Radii of the cylinders
% Len       Lengths of the cylinders
% Sta       Starting points of the cylinders
% Axe       Axes of the cylinders
% BOrd      Branch order data
% CiB       Cylinders in the branches
%
% Output:
% TreeData      Tree attributes


%% Tree attributes from cylinders
% Trunk cylinders
nc = length(Rad);
Trunk = false(nc,1);
Trunk(CiB{1}) = true;

% Volumes, lengths, and area
TotVol = 1000*pi*Rad.^2'*Len;
TrunkVol = 1000*pi*Rad(Trunk).^2'*Len(Trunk);
BranVol = 1000*pi*Rad(~Trunk).^2'*Len(~Trunk);
bottom = min(Sta(:,3));
[top,i] = max(Sta(:,3));
if Axe(i,3) > 0
    top = top+Len(i)*Axe(i,3);
end
TotHei = top-bottom;
TrunkLen = sum(Len(Trunk));
BranLen = sum(Len(~Trunk));
TotArea = 2*pi*sum(Rad.*Len);

% Diameter at breast height (dbh)
i = 1;
while sum(Len(1:i)) < 1.3
    i = i+1;
end
DBH = 200*Rad(i);

% Trunk taper and cumulative volume
n = nnz(Trunk);
Taper = zeros(n+1,4);
Taper(1,2) = Rad(1);
for i = 1:n
    Taper(i+1,:) = [sum(Len(1:i)) 200*Rad(i) 1000*pi*sum(Rad(1:i).^2.*Len(1:i)) 0];
end
Taper(:,4) = 100*Taper(:,3)/TrunkVol;

% Branch volume, length and frequency by order
b = length(BOrd)-1; % number of branches
BO = max(BOrd); % maximum branch order
BODist = zeros(BO,4);
for i = 1:BO
    I = BOrd == i;
    BODist(i,1) = sum(BVol(I)); % volumes
    BODist(i,2) = sum(BLen(I)); % lengths
    BODist(i,3) = nnz(I); % number of ith-order branches
    BODist(i,6) = 100*nnz(I)/(b-1); % relative value
end
V = sum(BODist(:,1));
L = sum(BODist(:,2));
BODist(:,4) = 100*BODist(:,1)/V; % relative values
BODist(:,5) = 100*BODist(:,2)/L;

% Branch size distribution
rad = Rad(~Trunk);
len = Len(~Trunk);
M = ceil(100*max(rad));
if M >= 3
    strX = 'Diameter (cm)';
    BSDist = zeros(M,4);
    for i = 1:M
        I = rad <= i/200;
        J = rad > (i-1)/200;
        I = I&J;
        BSDist(i,1) = 1000*pi*sum(rad(I).^2.*len(I)); % volumes
        BSDist(i,2) = sum(len(I)); % lengths
    end
else
    strX = 'Diameter (mm)';
    M = 10*M;
    BSDist = zeros(M,4);
    for i = 1:M
        I = rad <= i/2000;
        J = rad > (i-1)/2000;
        I = I&J;
        BSDist(i,1) = 1000*pi*sum(rad(I).^2.*len(I)); % volumes
        BSDist(i,2) = sum(len(I)); % lengths
    end
end
BSDist(:,3) = 100*BSDist(:,1)/BranVol; % relative values
BSDist(:,4) = 100*BSDist(:,2)/BranLen;

% Cumulative volume, length and frequency by order
CODist = zeros(BO+1,5);
CODist(:,1) = (0:1:BO)';
for i = 0:BO
    I = BOrd == i;
    CODist(i+1,2) = sum(BVol(I)); % volumes
    CODist(i+1,3) = sum(BLen(I)); % lengths
end
L = sum(CODist(:,3));
CODist(:,4) = 100*CODist(:,2)/TotVol; % relative values
CODist(:,5) = 100*CODist(:,3)/L;
CODist(:,2:5) = cumsum(CODist(:,2:5));

% Cumulative volume and length by size
M = ceil(100*max(Rad));
if M >= 3
    strS = 'Diameter (cm)';
    CSDist = zeros(M,4);
    for i = 1:M
        I = Rad <= i/200;
        J = Rad > (i-1)/200;
        I = I&J;
        CSDist(i,1) = 1000*pi*sum(Rad(I).^2.*Len(I)); % volumes
        CSDist(i,2) = sum(Len(I)); % lengths
    end
else
    strS = 'Diameter (mm)';
    M = 10*M;
    CSDist = zeros(M,4);
    for i = 1:M
        I = Rad <= i/2000;
        J = Rad > (i-1)/2000;
        I = I&J;
        CSDist(i,1) = 1000*pi*sum(Rad(I).^2.*Len(I)); % volumes
        CSDist(i,2) = sum(Len(I)); % lengths
    end
end
CSDist(:,3) = 100*CSDist(:,1)/TotVol; % relative values
CSDist(:,4) = 100*CSDist(:,2)/L;
CSDist = cumsum(CSDist);


%% Trunk volume and DBH from triangulation
% Determine suitable cylinders
n = nnz(Trunk);
i = 2;
while i < n && Rad(i) > 0.333*Rad(1) && Axe(i,:)*Axe(i-1,:)' > 0.97
    i = i+1;
end
i = i-1;
maxL = sum(Len(1:i));
maxL = round(100*maxL)/100;

% Set the parameters for triangulation
if maxL < 2;
    CL = 4*Rad(1);
    H = Rad(1)/2;
    NA = 18;
else
    if Rad(1) < 0.5
        CL = 1;
        H = 0.05;
        NA = 36;
    else
        CL = 2;
        H = 0.1;
        NA = 72;
    end
end

% Select the trunk point set used for triangulation
[~,~,h] = distances_to_line(trunk,Axe(1,:),Sta(1,:));
I = h < maxL+2*H;
trunk = trunk(I,:);

% Calculate the volumes
Vtcyl = 1000*pi*sum(Rad(1:i).^2.*Len(1:i));
[Vtrunk,Diam] = triangulated_cylinder_surface(trunk,NA,H,CL,maxL);

% Heights
Htri = round(Diam(end,4))/100; % Height of the triangulated surface
if Htri < maxL-0.2
    % if the triangulation was shortened, shorten the cylinder case also
    i = 2;
    while i < n && sum(Len(1:i)) < Htri 
        i = i+1;
    end
    i = i-1;
    maxL = sum(Len(1:i));
    maxL = round(100*maxL)/100;
    Vtcyl = 1000*pi*sum(Rad(1:i).^2.*Len(1:i));
end

% Dbh from triangulation
d = abs(Diam(:,end)-130);
[~,I] = min(d);
DBHtri = Diam(I,2);

%% Tree data
TreeData = zeros(33,1);         % Tree attributes
TreeData(1) = TotVol;           % Total volume of the tree
TreeData(2) = TrunkVol;         % Volume of the trunk
TreeData(3) = BranVol;          % Total volume of all the branches
TreeData(4) = TotHei;           % Total height of the tree
TreeData(5) = TrunkLen;         % Length of the trunk
TreeData(6) = BranLen;          % Total length of all the branches
TreeData(7) = b;                % Total number of branches
TreeData(8) = BO;               % Maximum branch order
TreeData(9) = TotArea;          % Total area of cylinders
TreeData(10) = DBH;             % Diameter at breast height (cylinder)
TreeData(11) = DBHtri;          % DBH from triangulation
TreeData(12) = Vtcyl;           % Trunk volume of over 33.3% diameter part (cylinders)
TreeData(13) = Vtrunk;          % Trunk volume of over 33.3% diameter part (triangulation)
TreeData(14) = maxL;            % Trunk length of over 33.3% diam part (cylinders)
TreeData(15) = Htri;            % Trunk length of over 33.3% diam part (triangulation)
for i = 1:min(6,BO)
    TreeData(15+i) =  BODist(i,3); % Number of ith-order branches
end
for i = 1:min(6,BO)
    TreeData(21+i) =  BODist(i,1); % Volume of ith-order branches
end
for i = 1:min(6,BO)
    TreeData(27+i) =  BODist(i,2); % Length of ith-order branches
end

% Round tree data
for i = 1:33
    D = TreeData(i);
    if D > 100
        D = round(D);
    elseif D > 10
        D = round(10*D)/10;
    elseif D > 1
        D = round(100*D)/100;
    else
        D = round(1000*D)/1000;
    end
    TreeData(i) = D;
end


%% Display tree data
disp('------------')
disp('Tree attributes:')
disp(['Total volume = ' num2str(TreeData(1)) ' liters'])
disp(['Trunk volume = ' num2str(TreeData(2)) ' liters'])
disp(['Branch volume = ' num2str(TreeData(3)) ' liters'])
disp(['Total height = ' num2str(TreeData(4)) ' meters'])
disp(['Trunk length = ' num2str(TreeData(5)) ' meters'])
disp(['Branch length = ' num2str(TreeData(6)) ' meters'])
disp(['Number of branches = ' num2str(TreeData(7))])
disp(['Maximum branch order = ' num2str(TreeData(8))])
disp(['Total cylinder area = ' num2str(TreeData(9)) ' square meters'])
disp(['Dbh (cylinder) = ' num2str(TreeData(10)) ' centimeters'])
disp(['Dbh (triangulation) = ' num2str(TreeData(11)) ' centimeters'])
disp(['Trunk volume (cylinders) = ' num2str(TreeData(12)) ' liters'])
disp(['Trunk volume (triangulation) = ' num2str(TreeData(13)) ' liters'])
disp(['Trunk length (cylinders) = ' num2str(TreeData(14)) ' meters'])
disp(['Trunk length (triangulation) = ' num2str(TreeData(15)) ' meters'])
disp('-----')
disp('Branch order data:')
S = ['1st'; '2nd'; '3rd'; '4th'; '5th'; '6th'];
for i = 1:min(6,BO)
    str = ['Number of ',S(i,:),'-order branches = ' num2str(TreeData(15+i))];
    disp(str)
end
for i = 1:min(6,BO)
    str = ['Volume of ',S(i,:),'-order branches = ' num2str(TreeData(21+i)),' liters'];
    disp(str)
end
for i = 1:min(6,BO)
    str = ['Length of ',S(i,:),'-order branches = ' num2str(TreeData(27+i)),' meters'];
    disp(str)
end
disp('------------')



% %% Plot and save distributions
% % Volumes
% h = plot_dist(BSDist(:,1),3,'Branch volume by size',strX,'Volume (L)');
% str = ['vol_size_abs_',string];
% saveas(h,str,'epsc')
% pause(0.1)
% h = plot_dist(BSDist(:,3),3,'Branch volume by size',strX,'Percent of total');
% str = ['vol_size_rel_',string];
% saveas(h,str,'epsc')
% pause(0.1)
% h = plot_dist(BODist(:,1),3,'Branch volume by order','Branch order','Volume (L)');
% str = ['vol_ord_abs_',string];
% saveas(h,str,'epsc')
% pause(0.1)
% h = plot_dist(BODist(:,4),3,'Branch volume by order','Branch order','Percent of total');
% str = ['vol_ord_rel_',string];
% saveas(h,str,'epsc')
% pause(0.1)
% % Lengths
% h = plot_dist(BSDist(:,2),3,'Branch length by size',strX,'Length (m)');
% str = ['len_size_abs_',string];
% saveas(h,str,'epsc')
% pause(0.1)
% h = plot_dist(BSDist(:,4),3,'Branch length by size',strX,'Percent of total');
% str = ['len_size_rel_',string];
% saveas(h,str,'epsc')
% pause(0.1)
% h = plot_dist(BODist(:,2),3,'Branch length by order','Branch order','Length (m)');
% str = ['len_ord_abs_',string];
% saveas(h,str,'epsc')
% pause(0.1)
% h = plot_dist(BODist(:,5),3,'Branch length by order','Branch order','Percent of total');
% str = ['len_ord_rel_',string];
% saveas(h,str,'epsc')
% pause(0.1)
% % Frequency
% h = plot_dist(BODist(:,3),3,'Number of branches by order','Branch order','Number');
% str = ['fre_ord_abs_',string];
% saveas(h,str,'epsc')
% pause(0.1)
% h = plot_dist(BODist(:,6),3,'Number of branches by order','Branch order','Percent of total');
% str = ['fre_ord_rel_',string];
% saveas(h,str,'epsc')
% pause(0.1)
% 
