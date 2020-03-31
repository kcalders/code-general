function [TreeData, Sta, Axe, Rad, Len, CPar, CExt, BoC, BOrd, BPar, BVol, BLen, BAng, BSeg, FCB,...
    BChi, CiB, CChi, CiS, Added, P, Bal, Cen, Bale, Nei, Segs, SPar,...
    SChi, SDir, SoC, Dim, Dir, Prin, Base, Forb]...
    = cylinder_tree_model(P, dmin, rcov, nmin, string, rfil1, nfil1, rfil2, nfil2)

% ---------------------------------------------------------------------
% CYLINDER_TREE_MODEL.M     Creates cylinder tree models from point clouds 
%                           scanned from trees.
%
% Version 1.0
% Author        Pasi Raumonen
% Created       14 June 2013
%
% This software cannot be distributed without the permission of P. Raumonen
% and it is only for non-commercial use. These restrictions holds also for 
% the modules or subprograms of the software:
% FILTERING.M, COVER_SETS.M, TREE_SETS.M, SEGMENTS.M, CORRECT_SEGMENTS.M,
% CYLINDERS.M, FILL_GAPS.M, BRANCHES.M, TREE_DATA.M
% ---------------------------------------------------------------------

% Produces a cylindrical tree model from point cloud, which is a sample of
% the tree surface. The point cloud is covered with small sets, which are
% along the surface, such that each point belongs at most one cover set; 
% i.e. the cover is a partition of the point cloud. 
% 
% The cover is generated such that at first the point cloud is covered 
% with balls with radius RCOV. This first cover is such that 
% 1) the minimum distance between the centers is DMIN, and 
% 2) the maximum distance from any point to nearest center is also DMIN.
% Then the first cover of RCOV-balls is used to define a second cover:
% each RCOV-ball A defines corresponding cover set B in the second cover
% such that B contains those points of A that are nearer to the center of
% A than any other center of RCOV-balls. The RCOV-balls also define 
% the neighbors for the second cover: Let CA and CB denote cover sets in 
% the and second cover, and BA and BB their RCOV-balls. Then CB is 
% a neighbor of CA, and vice versa, if BA and CB intersect or 
% BC and CA intersect. There is, however, exception to this rule: if the
% gap between the sets is too large, then they are not neighbors.
%
% The point cloud is first filtered, if RFIL1 > 0, to remove isolated 
% points and point groups. This is done as follows: First RFIL1-ball is
% generated for each point and then if a ball contains less than NFIL1 
% points, the point is removed. Next the remaining point cloud is covered 
% with RFIL2-balls, and the connected components are determined. If 
% a component has less than NFIL2 cover sets, then the points in the 
% component is removed.
%
% Inputs: 
% P         (Un)filtered point cloud, (m_points x 3)-matrix, the rows
%               give the coordinates of the points.
%               The order of the points is not meaningful
% dmin      Minimum distance between centers of cover sets; i.e. the
%               minimum diameter of a cover set
% rcov      Radius of the balls used to generate the cover sets, these 
%               balls are also used to determine the neighbors and the 
%               cover set characteristics              
% nmin      Minimum number of points in a rcov-ball
% string    Name string for saving output files
% rfil1     Radius of cover sets used in the first filtering process
% nfil1     Minimum number of points in the cover sets passing the first filtering
% rfil2     Radius of cover sets used in the second filtering process
% nfil2     Minimum number of cover sets in components passing the second filtering
% 
% Cylinder model outputs:
% Sta       Starting points of the cylinders, matrix
% Axe       Axes of the cylinders, matrix
% Rad       Radii of the cylinders, vector
% Len       Lengths of the cylinders, vector
% CPar      Parent cylinder of each cylinder, vector
% CExt      Extension cylinder of each cylinder, vector
% BoC       Branch of the cylinder, vector
% BOrd      Branch order, vector
% BPar      Parent branch, vector
% BVol      Volumes of the branches, vector
% BLen      Lengths of the branches, vector
% BAng      Branching angles of the branches, vector
% FCB       First cylinders of the branches, vector
%
% Additional outputs:
% TreeData  Vector containing basic tree attributes from the model
% BSeg      Segment of the branch, vector  (not every segment forms a branch)
% BChi      Child branches, cell-array
% CiB       Cylinders in the branches, cell-array
% CChi      Children cylinders of each cylinder, cell array
% CiS       Cylinders forming each segment, cell array
% Added     Logical vector indicating cylinders that are added to fill gaps
% P         Filtered point cloud, matrix
% Bal       Cover sets, cell array
% Cen       Center points of the cover sets, vector
% Bale      Extended cover sets, cell array
% Nei       Neighboring cover sets, cell array
% Segs      Tree segments, cell array
% SPar      Parent segment of each segment, vector
% SChi      Child segments of each segment, cell array
% SDir      Direction lines of the segment bases, matrix
% SoC       Segments the cylinders belong, vector
% Dim       Dimensionality values of the cover sets, matrix
% Dir       Direction lines of the conver sets, matrix
% Prin      Principal components of the cover sets, matrix
% Base      Base of the tree
% Forb      Cover sets not part of the tree, vector

% Names of the steps to display
name = ['Filtering  '; 
        'Cover sets ';
        'Tree sets  '; 
        'Segments   ';
        'Cylinders  '; 
        'Gap filling'];    
    
disp('Progress:')
tot = 0;

% only 3-dimensional data
if size(P,2) > 3
    P = P(:,1:3);
end


%% Filtering
if nargin > 5
    tic
    disp('FILTERING...')
    P = filtering(P,rfil1,nfil1,rfil2,nfil2);
    t = toc;
    tot = tot+t;
    str = [name(1,:),' ',num2str(t),' sec.',', total: ',num2str(tot),' sec.'];
    disp(str)
end


%% Cover sets
tic
disp('GENERATING COVER...')
[Bal,Cen,Nei,Bale,Dim,Dir,Prin] = cover_sets(P,dmin,rcov,nmin);
t = toc;
tot = tot+t;
str = [name(2,:),' ',num2str(t),' sec.',', total: ',num2str(tot),' sec.'];
disp(str)


%% Tree components and their bases
tic
disp('DETERMINING TREE SETS...')
[Base,Forb,Bale,Nei] = tree_sets(P,Cen,Bale,Nei,Dir,Dim,Prin,dmin,Bal);
t = toc;
tot = tot+t;
str = [name(3,:),' ',num2str(t),' sec.',', total: ',num2str(tot),' sec.'];
disp(str)


%% Segmenting
tic
disp('SEGMENTING...')
[Segs,SPar,SChi,SDir] = segments(P,Bal,Nei,Cen,Dir,Prin,Base,Forb,dmin);
[Segs,SPar,SChi,SDir] = correct_segments(P,Cen,Segs,SPar,SChi,SDir);
t = toc;
tot = tot+t;
str = [name(4,:),' ',num2str(t),' sec.',', total: ',num2str(tot),' sec.'];
disp(str)


%% Cylinders
tic
disp('CONSTRUCTING THE CYLINDER MODEL...')
[Rad,Len,Axe,Sta,CPar,CExt,CChi,SoC,CiS,Segs,SPar,SChi] = ...
    cylinders(P,Bal,Cen,Nei,Segs,SPar,SDir,SChi,dmin);
t = toc;
tot = tot+t;
str = [name(5,:),' ',num2str(t),' sec.',', total: ',num2str(tot),' sec.'];
disp(str)

%% Gap filling
tic
%Added = false(length(Rad),1);
Again = false;
disp('COMPLETING THE CYLINDER MODEL: FILLING GAPS...')
[Segs,SPar,SChi,CiS,Sta,Axe,Rad,Len,CPar,CExt,CChi,SoC,SDir,Added,Again]...
    = fill_gaps(Segs,SPar,SChi,CiS,Sta,Axe,Rad,Len,CPar,CExt,CChi,SoC,SDir,dmin,Again);
t = toc;
tot = tot+t;
str = [name(6,:),' ',num2str(t),' sec.',', total: ',num2str(tot),' sec.'];
disp(str)


if Again
    %% Fit cylinders again
    tic
    disp('CONSTRUCTING THE CYLINDER MODEL...')
    [Rad,Len,Axe,Sta,CPar,CExt,CChi,SoC,CiS,Segs,SPar,SChi] = ...
        cylinders(P,Bal,Cen,Nei,Segs,SPar,SDir,SChi,dmin);
    t = toc;
    tot = tot+t;
    str = [name(5,:),' ',num2str(t),' sec.',', total: ',num2str(tot),' sec.'];
    disp(str)
    
    
    %% Fill gaps again
    tic
    %Added = false(length(Rad),1);
    disp('COMPLETING THE CYLINDER MODEL: FILLING GAPS...')
    [Segs,SPar,SChi,CiS,Sta,Axe,Rad,Len,CPar,CExt,CChi,SoC,SDir,Added,Again]...
        = fill_gaps(Segs,SPar,SChi,CiS,Sta,Axe,Rad,Len,CPar,CExt,CChi,SoC,SDir,dmin,Again);
    t = toc;
    tot = tot+t;
    str = [name(6,:),' ',num2str(t),' sec.',', total: ',num2str(tot),' sec.'];
    disp(str)
    
end


%% Determine the branches
[BoC,BOrd,BPar,BChi,CiB,BVol,BLen,BAng,FCB,BSeg] = branches(Segs,SChi,CiS,CChi,CPar,Rad,Len,Axe,Added);


%% Plot models and compute and display model attributes
% plot pointcloud and cylinder model ADDED BY KIM
scatter3(P(:,1),P(:,2),P(:,3),0.01,'green')
hold on
plot_cylinder_model(Rad,Len,Axe,Sta,1,20,0.3)
hold on
hold off











%% Save the cylinder, branch, and tree data in text-files
CylData = [Rad Len Sta Axe CPar CExt BoC Added];
BranchData = [BOrd BPar BVol BLen BAng];

str = ['cyl_data_',string,'.txt'];
fid = fopen(str, 'wt');
fprintf(fid, [repmat('%g\t', 1, size(CylData,2)-1) '%g\n'], CylData.');
fclose(fid);

str = ['branch_data_',string,'.txt'];
fid = fopen(str, 'wt');
fprintf(fid, [repmat('%g\t', 1, size(BranchData,2)-1) '%g\n'], BranchData.');
fclose(fid);

str = ['tree_data_',string,'.txt'];
fid = fopen(str, 'wt');
fprintf(fid, [repmat('%g\t', 1, size(TreeData,2)-1) '%g\n'], TreeData.');
fclose(fid);

%save(string,'Sta','Axe','Rad','Len','CPar','CExt','BoC','BOrd','BPar','BVol',...
%    'BLen','BAng','BSeg','FCB','TreeData','BChi','CiB','CChi','CiS','Added','P',...
%    'Bal','Cen','Bale','Nei','Segs','SPar','SChi','SDir','SoC','Dim','Dir','Prin','Base','Forb');