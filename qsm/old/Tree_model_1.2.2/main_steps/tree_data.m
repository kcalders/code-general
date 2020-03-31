function [TreeData,Vert,Facets,fvd,trunk] = tree_data(Rad,Len,Sta,Axe,BOrd,CiB,trunk)

% ---------------------------------------------------------------------
% TREE_DATA.M       Calculates some tree attributes. 
%
% Version 1.2
% Author        Pasi Raumonen
% Created       2 Dec 2013
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


%% Define the trunk cylinders
TrunkCyls = CiB{1};
Trunk = false(length(Rad),1);
Trunk(TrunkCyls) = true;


%% Order frequency
b = length(BOrd);
BO = max(BOrd);
FreByOrder = zeros(BO+1,1);
for i = 0:BO
    I = BOrd == i;
    FreByOrder(i+1) = nnz(I)/b;
end


%% Compute the tree attributes "data"
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

if BO < 6
    FreByOrder(BO+2:7) = zeros(6-BO,1); 
end

TotArea = 2*pi*sum(Rad.*Len);

% Determine diameter at breast height (dbh) from cylinders and cones
i = 1;
while sum(Len(1:i)) < 1.3
    i = i+1;
end
DBH = 200*Rad(i);

% Calculate the trunk volume and DBH with triangulation
n = length(TrunkCyls);
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
[Vtrunk,Diam,~,Vert,Facets,fvd] = triangulated_cylinder_surface(trunk,NA,H,CL,maxL);
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
% Calculate the DBH from triangulation
d = abs(Diam(:,end)-130);
[~,I] = min(d);
DBHtri = Diam(I,2);
DBHtri = round(10*DBHtri)/10;


%% Round some attributes
if TotVol >= 100
    TotVol = round(TotVol);
elseif TotVol >= 10
    TotVol = round(10*TotVol)/10;
else
    TotVol = round(100*TotVol)/100;
end
if TrunkVol >= 100
    TrunkVol = round(TrunkVol);
elseif TrunkVol >= 10
    TrunkVol = round(10*TrunkVol)/10;
else
    TrunkVol = round(100*TrunkVol)/100;
end
if BranVol >= 100
    BranVol = round(BranVol);
elseif BranVol >= 10
    BranVol = round(10*BranVol)/10;
else
    BranVol = round(100*BranVol)/100;
end
if BranLen >= 100
    BranLen = round(BranLen);
elseif BranLen >= 10
    BranLen = round(10*BranLen)/10;
else
    BranLen = round(100*BranLen)/100;
end
TotHei = round(10*TotHei)/10;
TrunkLen = round(10*TrunkLen)/10;
TotArea = round(10*TotArea)/10;
DBH = round(10*DBH)/10;
if Vtcyl >= 100
    Vtcyl = round(Vtcyl);
elseif Vtcyl >= 10
    Vtcyl = round(10*Vtcyl)/10;
else
    Vtcyl = round(100*Vtcyl)/100;
end

%% Save the data
TreeData = zeros(21,1);         % Tree attributes
TreeData(1) = TotVol;           % Total volume of the tree
TreeData(2) = TrunkVol;         % Volume of the trunk
TreeData(3) = BranVol;          % Total volume of all the branches
TreeData(4) = TotHei;           % Total height of the tree
TreeData(5) = TrunkLen;         % Length of the trunk
TreeData(6) = BranLen;          % Total length of all the branches
TreeData(7) = b-1;              % Total number of branches
TreeData(8) = FreByOrder(2)*b;  % Number of 1st-order branches
TreeData(9) = FreByOrder(3)*b;  % Number of 2nd-order branches
TreeData(10) = FreByOrder(4)*b; % Number of 3rd-order branches
TreeData(11) = FreByOrder(5)*b; % Number of 4th-order branches
TreeData(12) = FreByOrder(6)*b; % Number of 5th-order branches
TreeData(13) = FreByOrder(7)*b; % Number of 6th-order branches
TreeData(14) = BO;              % Maximum branch order
TreeData(15) = TotArea;         % Total area of cylinders
TreeData(16) = DBH;             % Diameter at breast height (cylinder)
TreeData(17) = DBHtri;          % DBH from triangulation
TreeData(18) = Vtcyl;           % Trunk volume of over 33.3% diameter part (cylinders)
TreeData(19) = Vtrunk;          % Trunk volume of over 33.3% diameter part (triangulation)
TreeData(20) = maxL;            % Trunk length of over 33.3% diam part (cylinders)
TreeData(21) = Htri;            % Trunk length of over 33.3% diam part (triangulation)


%% Display the tree data
disp('------------')
disp(['Total volume = ' num2str(TotVol) ' liters'])
disp(['Trunk volume = ' num2str(TrunkVol) ' liters'])
disp(['Branch volume = ' num2str(BranVol) ' liters'])
disp(['Total height = ' num2str(TotHei) ' meters'])
disp(['Trunk length = ' num2str(TrunkLen) ' meters'])
disp(['Branch length = ' num2str(BranLen) ' meters'])
disp(['Number of branches = ' num2str(b-1)])
disp(['Number of 1st-order branches = ' num2str(FreByOrder(2)*b)])
disp(['Number of 2nd-order branches = ' num2str(FreByOrder(3)*b)])
disp(['Number of 3rd-order branches = ' num2str(FreByOrder(4)*b)])
disp(['Number of 4th-order branches = ' num2str(FreByOrder(5)*b)])
disp(['Number of 5th-order branches = ' num2str(FreByOrder(6)*b)])
disp(['Number of 6th-order branches = ' num2str(FreByOrder(7)*b)])
disp(['Maximum branch order = ' num2str(BO)])
disp(['Total cylinder area = ' num2str(TotArea) ' square meters'])
disp(['Dbh (cylinder) = ' num2str(DBH) ' centimeters'])
disp(['Dbh (triangulation) = ' num2str(DBHtri) ' centimeters'])
disp(['Trunk volume (cylinders) = ' num2str(Vtcyl) ' liters'])
disp(['Trunk volume (triangulation) = ' num2str(Vtrunk) ' liters'])
disp(['Trunk length (cylinders) = ' num2str(maxL) ' meters'])
disp(['Trunk length (triangulation) = ' num2str(Htri) ' meters'])
disp('------------')
