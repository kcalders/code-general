function [Partition,CubeCoord,CubeNei,Info] = ...
    partition_of_point_cloud(P,EdgeLength,Fast)

% Partitions the point cloud into cubes.
% 
% Inputs:
% P             Point cloud, (n_points x 3)-matrix
% EdgeLength    Length of the cube edges
% Fast          If true, uses always normal approach where the EdgeLength
%                   is increased so long as the total number of cubes is 
%                   small enough. If not true, uses always the given 
%                   EdgeLength and may use "sparse" approach, which is 
%                   much slower 
%
% Outputs:              
% Partition     Point cloud partitioned into cubical cells,
%                   (nx x ny x nz)-cell, where nx,ny,nz are the number
%                   of cubes in x,y,z-directions, respectively
% CubeCoord     (n_points x 4)-matrix whose rows are the cube coordinates 
%                   of each point: the three first numbers are the 
%                   x,y,z-coordinates and the fourth number is the ordinal 
%                   of the point in that cube
% CubeNei       Neighboring nonempty cubes for each cube in the "sparse"
%                   approach
% Info          The minimum coordinate values and number of cubes in each
%                   coordinate direction


if nargin == 2
    Fast = true;
end

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

while Fast && 8*nx*ny*nz > 3e9
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

if (nx*ny*nz < 4e8) || Fast
    % Small number of cubes in the whole cube space
    
    % Define the second output
    CubeCoord = [px py pz];
    CubeNei = zeros(0);

    % Initialization of "Partition"
    Partition = cell(nx,ny,nz);
    
    p = 1; % The index of the point under comparison
    while p <= np
        t = 1;
        while (p+t <= np) && (S(p) == S(p+t))
            t = t+1;
        end
        Partition{px(I(p)),py(I(p)),pz(I(p))} = I(p:p+t-1);
        p = p+t;
    end
    
else
    % Large number of cubes in the whole cube space, use "sparse" approach
    disp('N O T I C E !!!!!!!!')
    disp('Uses "sparse" approach because memory requirements')
    disp(['would be otherwise at least ' num2str(ceil(8*nx*ny*nz/1e6)) ' Megabytes'])
    
    % Define "Partition" and "CubeCoord" and other thing needed for
    % "CubeNei"
    Partition = cell(np,1);
    CubeCoord = zeros(np,1);
    Sheets = cell(nz,1); % xy-sheets of the cubes 
    Cubes = zeros(np,2); % info about nonempty cubes
    p = 1; % The index of the point under comparison
    q = 0;
    nxy = 20*ceil(np/nz);
    Sh = zeros(nxy,2);
    sn = pz(I(1));
    se = 0;
    while p <= np
        t = 1;
        while (p+t <= np) && (S(p) == S(p+t))
            t = t+1;
        end
        q = q+1;
        if sn == pz(I(p))
            se = se+1;
            Sh(se,:) = [S(p) q];
        else
            Sheets{sn} = Sh(1:se,:);
            Sh = zeros(nxy,2);
            Sh(1,:) = [S(p) q];
            se = 1;
            sn = pz(I(p));
        end
        Partition{q} = I(p:p+t-1);
        Cubes(q,:) = [pz(I(p)) S(p)];
        CubeCoord(I(p:p+t-1),1) = q;
        p = p+t;
    end
    Partition = Partition(1:q);
    Cubes = Cubes(1:q,:);
    
    CubeNei = cell(q,1);
    for i = 1:q
        C = Cubes(i,:);
        N = neighbors(C(2),nx,ny);
        M = zeros(27,1);
        Sh = Sheets{C(1)-1};
        if ~isempty(Sh)
            k = 1;
            n = length(Sh(:,1));
            for j = 1:9
                while (k <= n) && (Sh(k,1) ~= N(j))
                    k = k+1;
                end
                if k <= n
                   M(j) = Sh(k,2);
                end
            end
        end
        Sh = Sheets{C(1)};
        if ~isempty(Sh)
            k = 1;
            n = length(Sh(:,1));
            for j = 1:9
                while (k <= n) && (Sh(k,1) ~= N(j))
                    k = k+1;
                end
                if k <= n
                   M(j) = Sh(k,2);
                end
            end
        end
        Sh = Sheets{C(1)+1};
        if ~isempty(Sh)
            k = 1;
            n = length(Sh(:,1));
            for j = 1:9
                while (k <= n) && (Sh(k,1) ~= N(j))
                    k = k+1;
                end
                if k <= n
                   M(j) = Sh(k,2);
                end
            end
        end
        J = M > 0;
        CubeNei{i} = M(J);
    end
    
end
end % End of main function


function N = neighbors(I,nx,ny)

a = nx*ny;
N = [I-1-nx-a;
    I-nx-a;
    I+1-nx-a;
    I-1-a;
    I-a;
    I+1-a;
    I-1+nx-a;
    I+nx-a;
    I+1+nx-a;
    I-1-nx;
    I-nx;
    I+1-nx;
    I-1;
    I;
    I+1;
    I-1+nx;
    I+nx;
    I+1+nx;
    I-1-nx+a;
    I-nx+a;
    I+1-nx+a;
    I-1+a;
    I+a;
    I+1+a;
    I-1+nx+a;
    I+nx+a;
    I+1+nx+a];

end