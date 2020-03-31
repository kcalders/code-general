function [D,V] = dimensionality(points,varargin)

% Calculates the dimensionality estimates of "points" using principal 
% components. Returns also the corresponding principal direction vectors.
% The "points" used in the calculations must be a finite subset of R^3 and 
% this can be given directly as the input or as the indexes of points in 
% a point cloud or as the indexes of the cover sets in a point cloud.
%
% Inputs:
% points    Subset of R^3, or the indexes of point given in the optional
%           second input, or the indexes of the cover sets in the optional
%           third input
% P         Point cloud, optional second input
% Bal       Cover sets of P, optional third input
%
% Outputs;
% D         Dimensionality estimates of the "points", (1x3)-vector
% V         Principal component directions of the "points", (1x9)-vector


% Define the subset of R^3 used for calculations, if necessary 
if nargin == 2
    % second input is the point cloud "P" and "points" are the indexes 
    % of the points in "P"
    P = varargin{1};    
    points = P(points,:);
elseif nargin == 3
    % second and third inputs are the point cloud "P" and cover sets "Bal".
    % "points" are the indexes of the cover sets. 
    P = varargin{1};
    Bal = varargin{2};
    I = unique(vertcat(Bal{points}));
    points = P(I,:);
end


% Calculate the principal components using the eigenvectors and eigenvalues 
% of the covariance matrix
X = cov(points);
[U,S,~] = svd(X);

% Dimensionality estimates from the eigenvalues
D = [(S(1,1)-S(2,2))/S(1,1) (S(2,2)-S(3,3))/S(1,1) S(3,3)/S(1,1)];

% Principal component direction vectors, the eigenvectors
if nargout == 2
    V = [U(:,1)' U(:,2)' U(:,3)'];
end