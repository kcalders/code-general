function [Components,CompSize] = connected_components(Nei,Sub,MinSize,Fal)

% Determines connected components of the subset of cover sets defined
% by "Sub" such that each component has at least "MinSize"
% number of cover sets.
%
% Inputs:
% Nei       Neighboring cover sets of each cover set, (n_sets x 1)-cell
% Sub       Subset whose components are determined,
%               length(Sub) < 2 means no subset and thus the whole point cloud
%               "Sub" may be also a vector of cover set indexes in the subset
%               or a logical (n_sets)-vector, where n_sets is the number of
%               all cover sets
% MinSize   Minimum number of cover sets in an acceptable component
%
% Outputs:
% Components    Connected components, (n_comp x 1)-cell
% CompSize      Number of sets in the components, (n_comp x 1)-vector

if isempty(Sub)
    Components = cell(0,1);
    CompSize = zeros(0,1);
elseif length(Sub) <= 3 && ~islogical(Sub) && Sub(1) > 0
    % Very small subset, i.e. at most 3 cover sets
    n = length(Sub);
    if n == 1
        Components = cell(1,1);
        Components{1} = Sub;
        CompSize = 1;
    elseif n == 2
        I = Nei{Sub(1)} == Sub(2);
        if any(I)
            Components = cell(1,1);
            Components{1} = Sub;
            CompSize = 1;
        else
            Components = cell(2,1);
            Components{1} = Sub(1);
            Components{2} = Sub(2);
            CompSize = [1 1];
        end
    elseif n == 3
        I = Nei{Sub(1)} == Sub(2);
        J = Nei{Sub(1)} == Sub(3);
        K = Nei{Sub(2)} == Sub(3);
        if any(I)+any(J)+any(K) >= 2
            Components = cell(1,1);
            Components{1} = Sub;
            CompSize = 1;
        elseif any(I)
            Components = cell(2,1);
            Components{1} = Sub(1:2);
            Components{2} = Sub(3);
            CompSize = [2 1];
        elseif any(J)
            Components = cell(2,1);
            Components{1} = Sub([1 3]);
            Components{2} = Sub(2);
            CompSize = [2 1];
        elseif any(K)
            Components = cell(2,1);
            Components{1} = Sub(2:3);
            Components{2} = Sub(1);
            CompSize = [2 1];
        else
            Components = cell(3,1);
            Components{1} = Sub(1);
            Components{2} = Sub(2);
            Components{3} = Sub(3);
            CompSize = [1 1 1];
        end
    end
    
else
    nb = size(Nei,1);
    if length(Sub) == 1 && Sub == 0
        % All the cover sets
        ns = nb;
        if nargin == 3
            Sub = true(nb,1);
        else
            Sub = ~Fal;
        end
    elseif ~islogical(Sub)
        % Subset of cover sets
        ns = length(Sub);
        if nargin == 3
            sub = false(nb,1);
        else
            sub = Fal;
        end
        sub(Sub) = true;
        Sub = sub;
    else
        % Subset of cover sets
        ns = nnz(Sub);
    end
    
    Components = cell(ns,1);
    CompSize = zeros(ns,1);
    nc = 0;      % number of components found
    m = 1;
    while ~Sub(m)
        m = m+1;
    end
    i = 0;
    while i < ns
        Added = Nei{m};
        I = Sub(Added);
        Added = Added(I);
        a = length(Added);
        Comp = zeros(ns,1);
        Comp(1) = m;
        Sub(m) = false;
        t = 1;
        while a > 0
            Comp(t+1:t+a) = Added;
            Sub(Added) = false;
            t = t+a;
            Ext = unique(vertcat(Nei{Added}));
            I = Sub(Ext);
            Added = Ext(I);
            a = length(Added);
        end
        i = i+t;
        if t >= MinSize
            nc = nc+1;
            Components{nc} = Comp(1:t);
            CompSize(nc) = t;
        end
        if i < ns
            while m <= nb && Sub(m) == false
                m = m+1;
            end
        end
    end
    Components = Components(1:nc);
    CompSize = CompSize(1:nc);
end