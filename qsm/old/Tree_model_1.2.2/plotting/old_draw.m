function draw(radius,length,axis,start_point,varargin)
% Draw a branch (cylinder). The additional arguments in
% varargin are passed to the cylinder-command which computes
% points for a unit cylinder. The first argument is assumed
% to be the number of points to be used in the circles of
% the cylinder. The second is a generator curve which is
% used for the profile of the cylinder. For more details
% about these arguments SEE CYLINDER.


% Pass arguments to the cylinder-command.
if nargin == 6
    [x y z] = cylinder(1,varargin{2});
elseif nargin == 7
    [x y z] = cylinder(varargin{3},varargin{2});
else
    [x y z] = cylinder();
end

% Scale the radius and height of the unit cylinder with the
% corresponding parameters of the branch.
x = radius*x;
y = radius*y;
z = length*z;

% Axis of the unit cylinder.
u = [0 0 1];

% Compute rotation axis and normalize it.
raxis = cross(u,axis);
raxis = raxis/norm(raxis);
if isnan(raxis)
    raxis = [1 0 0];
end

% Compute rotation angle.
angle = acos(dot(u,axis)/norm(axis));

% Store size of the cylinder points for later reshaping.
s = size(x);

% Vercat point matrices for faster computations. Multiply
% with the rotation matrix.
X = [x(:) y(:) z(:)]*rotation_matrix(raxis,angle)';

% Reshape to size which is required by the surf-command.
x = reshape(X(:,1),s);
y = reshape(X(:,2),s);
z = reshape(X(:,3),s);

% Transfer computed cylinder to the starting point of the
% branch.
x = x + start_point(1);
y = y + start_point(2);
z = z + start_point(3);

% Draw the surface of the branch.
surf(x,y,z)
