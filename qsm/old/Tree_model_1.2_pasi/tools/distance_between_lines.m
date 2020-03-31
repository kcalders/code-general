function [d,C,E] = distance_between_lines(P,D,Q,V)

% Calculates the distance from line (Q,V) (Q is point in the line and V 
% is a direction vector V) to the lines (P,D).
% P and D may be vectors whose rows define the lines.


C = cross_prod(D,V); % cross products between the direction vectors
l = sqrt(sum(C.*C,2)); % lengths of the cross products
C = [C(:,1)./l C(:,2)./l C(:,3)./l];  % make unit vectors

A = mat_vec_subtraction(P,Q); % vectors between line points
B = sum(A.*C,2); % dot product between corresponding rows of A and C

d = abs(B); % the distances between lines

n = size(P,1);
E = zeros(n,3);
for i = 1:n
    A = [V' -D(i,:)' C(i,:)'];
    b = A\(P(i,:)-Q)';
    E(i,:) = b';
end

disp([d E])