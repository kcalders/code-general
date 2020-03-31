function plot_segs(P,comps,fig,Bal)

% Plots the segments given in the cell array "comps".
% If 3 input arguments, cells contain the point indexes, but if 4 input
% arguments, cells contain the indexes of the cover sets.
% "fig" is the figure number.

col = [
	0.00  0.00  1.00
	0.00  0.50  0.00
	1.00  0.00  0.00
	0.00  0.75  0.75
	0.75  0.00  0.75
	0.75  0.75  0.00
	0.25  0.25  0.25
	0.75  0.25  0.25
	0.95  0.95  0.00
	0.25  0.25  0.75
	0.75  0.75  0.75
	0.00  1.00  0.00
	0.76  0.57  0.17
	0.54  0.63  0.22
	0.34  0.57  0.92
	1.00  0.10  0.60
	0.88  0.75  0.73
	0.10  0.49  0.47
	0.66  0.34  0.65
	0.99  0.41  0.23];

% col = [0 0 1;
%        0 1 0;
%        1 0 0;
%        0 1 1;
%        1 1 0;
%        1 1 1];
 
col = rand(10000,3);

S = comps{1};
if iscell(S)
    n = size(comps,1);
    for i = 1:n
        S = comps{i};
        if ~isempty(S)
            S = vertcat(S{:});
            comps{i} = S;
        else
            comps{i} = zeros(0,1);
        end
    end
end

ms = 8; % Markersize

if nargin == 3
    % Define the colors for plotting segments
    n = max(size(comps));
    col(n,:) = zeros(1,3);
    if n > 20
        col = rand(n,3);
    end
    
    % Plot the segments
    figure(fig)
    C = comps{1};
    plot3(P(C,1),P(C,2),P(C,3),'.','Color',col(1,:),'Markersize',ms)
    hold on
    for i = 2:n
        C = comps{i};
        plot3(P(C,1),P(C,2),P(C,3),'.','Color',col(i,:),'Markersize',ms)
    end
    axis equal
    hold off
    pause(0.1)
    
else
    
    n = max(size(comps));
    col(n,:) = zeros(1,3);
    if n > 20
        col = rand(n,3);
    end
    
    D = zeros(0,1);    
    C = unique(vertcat(Bal{comps{1}}));
    figure(fig)
    plot3(P(C,1),P(C,2),P(C,3),'.','Color',col(1,:),'Markersize',ms)
    hold on
    for i = 2:n
        C = unique(vertcat(Bal{comps{i}}));
        C = setdiff(C,D);
        D = union(D,C);
        plot3(P(C,1),P(C,2),P(C,3),'.','Color',col(i,:),'Markersize',ms)
    end
    hold off
    axis equal
    pause(0.1)
end