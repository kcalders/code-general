function [TreeData, Sta, Axe, Rad, Len, CPar, CExt, BoC, BOrd, BPar, ...
    BVol, BLen, BAng, BSeg, FCB, BChi, CiB, CChi, CiS, Added, P, Bal, ...
    Cen, Bale, Nei, Segs, SPar, SChi, SDir, SoC, Base, Forb]...
    = cylinder_tree_model(P, dmin, rcov, nmin, lcyl, NoGround, string,...
        rfil1, nfil1, rfil2, nfil2)

% ---------------------------------------------------------------------
% CYLINDER_TREE_MODEL.M     Creates cylinder tree models from point clouds 
%                           scanned from trees.
%
% Version 1.2
% Author        Pasi Raumonen
% Created       2 Dec 2013
%
% This software cannot be distributed without the permission of P. Raumonen
% and it is only for non-commercial use. These restrictions holds also for 
% the modules or subprograms of the software:
% FILTERING.M, COVER_SETS.M, TREE_SETS.M, SEGMENTS.M, CORRECT_SEGMENTS.M,
% CYLINDERS.M, FILL_GAPS.M, BRANCHES.M, TREE_DATA.M
% ---------------------------------------------------------------------

% Kim & Andy, new features and changes:
% 1) COVER_SETS don't calculate dimensionality, direction, principal
%       components of the sets anymore, therefore
% 2) TREE_SETS finds the trunk differently, and
% 3) if the ground is not part of the point cloud, the base of the trunk
%       can be determined faster and more reliably just by taking the bottom 
%       part of the point cloud. This means that the input parameter NOGROUND
%       is "true", if it is "false", then tries to separate the ground.
% 4) SEGMENTS do not threat the trunk separately, instead
% 5) CORRECT_SEGMENTS defines the trunk from top to bottom, and in this
%       process changes the segments accordingly.
% 6) CORRECT_SEGMENTS removes lot of sets from the parent segment near the
%       base of the child segment. This way particularly the forks, but also
%       other branches do not generate too big cylinders in the parent
%       segment where the child segments originate.
% 7) Another new input parameter, LCYL, defines the relative length of the
%       cylinders, LCYL = (cylinder lenght)/(cylinder radius). If this
%       parameter is under 6, then at least for trunk the cylinders are
%       fitted again in the subcode "refinement_fitting" so that the
%       relative length of the cylinders is about LCYL. Suitable value
%       usually is about 3-5.
% 8) BRANCHES will modify each branch so that the radius of the cylinders
%       along the branch is always decreasing.
% 9) TREE_DATA caculates for the trunk also the volume and DBH from
%       triangulation. The triangulation is only for the part with radius
%       larger than 33% of the first cylinder and small angles between
%       cylinders. The volume, DBH, and length of the triangulated part are
%       in the TREEDATA vector, as is also for comparison the corresponding
%       cylinder model values.
% 10) "CylData" text-output file contains now the following info:
%       Radius  Length  Starting_point(x,y,z)  Axis(x,y,z)  Parent  Extension  Branch  Branch_order  Index_inside_branch  Added        


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
% lcyl      Cylinder length/radius ratio
% NoGround  Logical value, true if no ground in the point cloud, in which
%               case defines the base of the trunk as the lowest part the
%               cloud
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
% Base      Base of the tree
% Forb      Cover sets not part of the tree, vector

% Names of the steps to display
name = ['Filtering  ';
        'Cover sets ';
        'Tree sets  ';
        'Segments   ';
        'Cylinders  '];
    
disp('---------------')
disp(string)
str = ['dmin = ',num2str(dmin),', rcov = ',num2str(rcov),', nmin = ',num2str(nmin),', lcyl = ',num2str(lcyl)];
disp(str)
disp('Progress:')
tot = 0;

% only 3-dimensional data
if size(P,2) > 3
    P = P(:,1:3);
end


%% Filtering
if nargin > 7
    tic
    disp('FILTERING...')
    P = filtering(P,rfil1,nfil1,rfil2,nfil2);
    t = toc;
    tot = tot+t;
    [tmin,tsec] = sec2min(t);
    [Tmin,Tsec] = sec2min(tot);
    str = [name(1,:),' ',num2str(tmin),' min ',num2str(tsec),' sec',', total: ',num2str(Tmin),' min ',num2str(Tsec),' sec'];
    disp(str)
end


%% Cover sets
tic
disp('GENERATING COVER...')
[Bal,Cen,Nei,Bale] = cover_sets(P,dmin,rcov,nmin);
t = toc;
tot = tot+t;
[tmin,tsec] = sec2min(t);
[Tmin,Tsec] = sec2min(tot);
str = [name(2,:),' ',num2str(tmin),' min ',num2str(tsec),' sec',', total: ',num2str(Tmin),' min ',num2str(Tsec),' sec'];
disp(str)


%% Tree sets
tic
disp('DETERMINING TREE SETS...')
[Base,Forb,Bale,Nei] = tree_sets(P,Cen,Bal,Bale,Nei,dmin,NoGround);
t = toc;
tot = tot+t;
[tmin,tsec] = sec2min(t);
[Tmin,Tsec] = sec2min(tot);
str = [name(3,:),' ',num2str(tmin),' min ',num2str(tsec),' sec',', total: ',num2str(Tmin),' min ',num2str(Tsec),' sec'];
disp(str)


%% Segmenting
tic
disp('SEGMENTING...')
[Segs,SPar,SChi,SDir] = segments(P,Bal,Nei,Cen,Base,Forb,dmin);

[Segs,SPar,SChi,SDir] = correct_segments(P,Cen,Segs,SPar,SChi,SDir,Bal,dmin);
t = toc;
tot = tot+t;
[tmin,tsec] = sec2min(t);
[Tmin,Tsec] = sec2min(tot);
str = [name(4,:),' ',num2str(tmin),' min ',num2str(tsec),' sec',', total: ',num2str(Tmin),' min ',num2str(Tsec),' sec'];
disp(str)


%% Cylinders
tic
disp('CONSTRUCTING THE CYLINDER MODEL...')
[Rad,Len,Axe,Sta,CPar,CExt,CChi,SoC,CiS,Segs,SPar,SChi,Again1] = ...
    cylinders(P,Bal,Cen,Segs,SPar,SDir,SChi,dmin,lcyl);

[Segs,SPar,SChi,CiS,Sta,Axe,Rad,Len,CPar,CExt,CChi,SoC,SDir,Added,Again2]...
    = fill_gaps(Segs,SPar,SChi,CiS,Sta,Axe,Rad,Len,CPar,CExt,CChi,SoC,SDir,dmin);
t = toc;
tot = tot+t;
[tmin,tsec] = sec2min(t);
[Tmin,Tsec] = sec2min(tot);
str = [name(5,:),' ',num2str(tmin),' min ',num2str(tsec),' sec',', total: ',num2str(Tmin),' min ',num2str(Tsec),' sec'];
disp(str)

Again = Again1 || Again2;
if Again
    disp('CHANGES IN SEGMENTATION, FIT CYLINDERS AGAIN')
    %% Fit cylinders again
    tic
    disp('CONSTRUCTING THE CYLINDER MODEL...')
    [Rad,Len,Axe,Sta,CPar,CExt,CChi,SoC,CiS,Segs,SPar,SChi] = ...
    cylinders(P,Bal,Cen,Segs,SPar,SDir,SChi,dmin,lcyl);

    [Segs,SPar,SChi,CiS,Sta,Axe,Rad,Len,CPar,CExt,CChi,SoC,SDir,Added]...
        = fill_gaps(Segs,SPar,SChi,CiS,Sta,Axe,Rad,Len,CPar,CExt,CChi,SoC,SDir,dmin);
    t = toc;
    tot = tot+t;
    [tmin,tsec] = sec2min(t);
    [Tmin,Tsec] = sec2min(tot);
    str = [name(5,:),' ',num2str(tmin),' min ',num2str(tsec),' sec',', total: ',num2str(Tmin),' min ',num2str(Tsec),' sec'];
    disp(str)
    
end


%% Determine the branches
[Rad,BoC,BOrd,BPar,BChi,CiB,BVol,BLen,BAng,FCB,BSeg] = ...
    branches(Segs,SChi,CiS,CChi,CPar,Rad,Len,Axe,Added);


%% Compute and display model attributes
T = Segs{1};
T = vertcat(T{:});
T = vertcat(Bal{T});
Trunk = P(T,:);
[TreeData,Vert,Facets,fvd,Trunk] = tree_data(Rad,Len,Sta,Axe,BOrd,CiB,Trunk);


%% Plot models
plot_branch_structure(P,Bal,Segs,BSeg,BChi,1,3,0,1);
plot_tree_structure(P,Bal,Segs,SChi,4,4,0,1)
plot_cylinder_model(Rad,Len,Axe,Sta,2,20,0.3)

% plot the segments and cylinders in the same figure
plot_branch_structure(P,Bal,Segs,BSeg,BChi,3,3,0,1);
hold on
plot_cylinder_model(Rad,Len,Axe,Sta,3,20,0.3)
hold off

% plot part of the trunk segment and its triangulation 
if TreeData(19) > 0
    point_cloud_plotting(Trunk,3,5)
    patch('Vertices',Vert,'Faces',Facets,'FaceVertexCData',fvd,'FaceColor','flat')
    axis equal
    alpha(0.4)
end



%% Save the cylinder, branch, and tree data in text-files
rad = round(10000*Rad)/10000;
len = round(10000*Len)/10000;
sta = round(10000*Sta)/10000;
axe = round(10000*Axe)/10000;
bvol = round(1000*BVol)/1000;
blen = round(1000*BLen)/1000;
bang = round(1000*BAng)/1000;
CylData = [rad len sta axe CPar CExt BoC Added];
BranchData = [BOrd BPar bvol blen bang];

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

% save(string,'TreeData','Sta','Axe','Rad','Len','CPar','CExt','BoC','BOrd','BPar', ...
%     'BVol','BLen','BAng','BSeg','FCB','BChi','CiB','CChi','CiS','Added','P','Bal', ...
%     'Cen','Bale','Nei','Segs','SPar','SChi','SDir','SoC','Base','Forb');