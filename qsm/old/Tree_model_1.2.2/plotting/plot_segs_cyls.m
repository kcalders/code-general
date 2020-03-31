function plot_segs_cyls(P,Bal,Segs,SChi,Rad,Len,Axe,Sta,fig,ms,gen,segnum,ns,alp)

% Plots a segmented point cloud and the corresponding cylinder model into same plot

% Inputs:
% P 		Point cloud
% Bal 		Cover sets
% Segs 		Segments
% SChi		Child segments
% Rad 		Radii of the cylinders
% Len		Lengths of the cylinders
% Axe 		Axes of the cylinders
% Sta		Starting points of the cylinders
% fig 		Figure number
% ms 		Marker size for point cloud points
% gen 		Number of branch orders used for plotting: 0 == all, 1 = trunk, 
%				2 = trunk & 1st-order branches, 3 = trunk & 1st- & 2nd-order branches, etc.
% segnum 	Segment number from which the segmented point cloud starts, 1 = trunk
% ns 		Number of "sides" used for plotting the largest cylinder (e.g. 30), decreases
% 				linearly for radii but minimum of 6 sides for cylinder plotting
% alp 		"Alpha" or how opaque the cylinders are (between 0 and 1 = cannot see through)

% Plots all the cylinders regardles of the plotted segmens


plot_tree_structure(P,Bal,Segs,SChi,fig,ms,gen,segnum)
hold on
plot_cylinder_model(Rad,Len,Axe,Sta,fig,ns,alp)
hold off