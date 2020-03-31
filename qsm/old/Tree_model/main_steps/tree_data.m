function TreeData = tree_data(Rad,Len,Sta,Axe,BOrd,CiB)

% ---------------------------------------------------------------------
% TREE_DATA.M       Calculates some tree attributes. 
%
% Version 1.0
% Author        Pasi Raumonen
% Created       14 June 2013
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


%% Compute the tree attributes "data" and display them
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

TreeData = zeros(14,1);          % Tree attributes
TreeData(1) = TotVol;           % Total volume of the tree
TreeData(2) = TrunkVol;         % Volume of the trunk
TreeData(3) = BranVol;          % Total volume of all the branches
TreeData(4) = TotHei;           % Total height of the tree
TreeData(5) = TrunkLen;         % Length of the trunk
TreeData(6) = BranLen;          % Total length of all the branches
TreeData(7) = FreByOrder(2)*b;  % Number of 1st-order branches
TreeData(8) = FreByOrder(3)*b;  % Number of 2nd-order branches
TreeData(9) = FreByOrder(4)*b;  % Number of 3rd-order branches
TreeData(10) = FreByOrder(5)*b;  % Number of 4th-order branches
TreeData(11) = FreByOrder(6)*b;  % Number of 5th-order branches
TreeData(12) = FreByOrder(7)*b;  % Number of 6th-order branches
TreeData(13) = b-1;              % Total number of branches
TreeData(14) = BO;               % Maximum branch order

% Round some attributes for the display
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


% display the tree data
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
disp('------------')