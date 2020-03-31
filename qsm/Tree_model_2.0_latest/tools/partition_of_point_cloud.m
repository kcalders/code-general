function [Partition,CC,Info] = partition_of_point_cloud(P,EL)

% Partitions the point cloud into cubes.
%
% Version 1.3
% Author        Pasi Raumonen
% Created       21 February 2014
% 
% Inputs:
% P             Point cloud, (n_points x 3)-matrix
% EL            Length of the cube edges
%
% Outputs:              
% Partition     Point cloud partitioned into cubical cells,
%                   (nx x ny x nz)-cell, where nx,ny,nz are the number
%                   of cubes in x,y,z-directions, respectively
% CC            (n_points x 3)-matrix whose rows are the cube coordinates 
%                   of each point: x,y,z-coordinates
% Info          The minimum coordinate values and number of cubes in each
%                   coordinate direction


% The vertices of the big cube containing P
Min = double(min(P));
Max = double(max(P));

% Number of cubes with edge length "EdgeLength" in the sides 
% of the big cube
N = double(ceil((Max-Min)/EL)+7);
while 8*N(1)*N(2)*N(3) > 4e9
    EL = 1.1*EL;
    N = double(ceil((Max-Min)/EL)+7);
end
Info = [Min N EL];

% Calculates the cube-coordinates of the points
CC = floor([P(:,1)-Min(1) P(:,2)-Min(2) P(:,3)-Min(3)]/EL)+4;

% Sorts the points according a lexicographical order
S = [CC(:,1) CC(:,2)-1 CC(:,3)-1]*[1 N(1) N(1)*N(2)]';
CC = uint16(CC);
[S,I] = sort(S);
I = int32(I);
S = int32(S);

% Define "Partition"
Partition = cell(N(1),N(2),N(3));
np = size(P,1);     % number of points
p = 1;              % The index of the point under comparison
while p <= np
    t = 1;
    while (p+t <= np) && (S(p) == S(p+t))
        t = t+1;
    end
    q = I(p);
    Partition{CC(q,1),CC(q,2),CC(q,3)} = I(p:p+t-1);
    p = p+t;
end
