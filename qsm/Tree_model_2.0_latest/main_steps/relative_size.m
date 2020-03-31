function [P,RS] = relative_size(P,Cen,Bal,Segs,SChi)

% ---------------------------------------------------------------------
% RELATIVE_SIZE.M   Determines relative cover set size for points in new covers
%
% Version 1.0
% Author        Pasi Raumonen
% Created       15 March 2014
%
% This software cannot be distributed without the permission of P. Raumonen
% and it is only for non-commercial use.
% ---------------------------------------------------------------------
% 
% Uses existing segmentation and its branching structure to determine
% relative size of the cover sets distributed over new covers. The idea is 
% to decrease the relative size as the branch size decreases. This is 
% realised so that the relative size decreases when the branch order 
% increses, the height of the branch increases, and when we approach the
% tip of the branch. Maximum relative size is 256 at the bottom of the
% stem and the minimum is 1 at the tip of every branch.
%
% Output:
% P     Point cloud containing only the segmented points
% RS    Relative size (1-256), uint8-vector, (n_points x 1)


np = size(P,1);     % number of points
ns = size(Segs,1);  % number of segments

% Determine the branch orders of the segments
Ord = zeros(ns,1);
C = SChi{1};
i = 0;
while ~isempty(C)
    i = i+1;
    Ord(C) = i;
    C = vertcat(SChi{C});
end
maxO = i+1; % maximum branching order (plus one)

% Determine tree height
Top = max(P(Cen,3));
Bot = min(P(Cen,3));
H = Top-Bot;

TS = 1;
RS = zeros(np,1,'uint8');
for i = 1:ns
    O = Ord(i);
    S = Segs{i};
    Q = S{1};
    bot = min(P(Cen(Q),3));
    h = bot-Bot;
    BS = ceil(256*(maxO-O)/maxO*(H-h)/H);
    s = size(S,1);
    for j = 1:s
        Q = S{j};
        RS(vertcat(Bal{Q})) = BS-(BS-TS)*(j-1)/s;
    end
end

% Use only the segmented points
I = RS > 0;
P = P(I,:);
RS = RS(I,:);