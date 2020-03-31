function v = optimal_orthogonal_vector(V)

% For a given set of unit vectors (the rows of the matrix "V"),
% returns a unit vector ("v") that is the most orthogonal to them all
% in the sense that the sum of squared dot products of v with the
% vectors of V is minimized.
 
A = V'*V;
[U,~,~] = svd(A);
v = U(:,3)';