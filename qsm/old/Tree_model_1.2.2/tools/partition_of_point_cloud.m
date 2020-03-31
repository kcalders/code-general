function [Partition,CubeCoord,Info] = partition_of_point_cloud(P,EdgeLength)

% Partitions the point cloud into cubes.
% 
% Inputs:
% P             Point cloud, (n_points x 3)-matrix
% EdgeLength    Length of the cube edges
%
% Outputs:              
% Partition     Point cloud partitioned into cubical cells,
%                   (nx x ny x nz)-cell, where nx,ny,nz are the number
%                   of cubes in x,y,z-directions, respectively
% CubeCoord     (n_points x 3)-matrix whose rows are the cube coordinates 
%                   of each point: x,y,z-coordinates
% Info          The minimum coordinate values and number of cubes in each
%                   coordinate direction


% The vertices of the big cube containing P
x_min = min(P(:,1));
x_max = max(P(:,1));
y_min = min(P(:,2));
y_max = max(P(:,2));
z_min = min(P(:,3));
z_max = max(P(:,3));

np = length(P(:,1));    % number of points in the point cloud

% Number of cubes with edge length "EdgeLength" in the sides 
% of the big cube
nx = ceil((x_max-x_min)/EdgeLength)+3;
ny = ceil((y_max-y_min)/EdgeLength)+3;
nz = ceil((z_max-z_min)/EdgeLength)+3;

while 8*nx*ny*nz > 4e9
    EdgeLength = 1.1*EdgeLength;
    nx = ceil((x_max-x_min)/EdgeLength)+3;
    ny = ceil((y_max-y_min)/EdgeLength)+3;
    nz = ceil((z_max-z_min)/EdgeLength)+3;
end

Info = [x_min y_min z_min nx ny nz EdgeLength];

% Calculates the cube-coordinates of the points
px = (P(:,1)-x_min)/EdgeLength;
py = (P(:,2)-y_min)/EdgeLength;
pz = (P(:,3)-z_min)/EdgeLength;

px = floor(px)+2;
py = floor(py)+2;
pz = floor(pz)+2;

% Sorts the points according a lexicographical order
index = [px py-1 pz-1]*[1 nx nx*ny]';
[S,I] = sort(index);

% Define the second output
CubeCoord = [px py pz];

% Initialization of "Partition"
Partition = cell(nx,ny,nz);

p = 1; % The index of the point under comparison
while p <= np
    t = 1;
    while (p+t <= np) && (S(p) == S(p+t))
        t = t+1;
    end
    q = I(p);
    Partition{px(q),py(q),pz(q)} = I(p:p+t-1);
    p = p+t;
end
