function point_cloud_plotting(P,fig,ms,Bal,Sub)

% Plots the given point cloud "P". With additional inputs we can plot only
% those points that are included in the cover sets "Bal" or in the
% subcollection of cover sets "Sub"
% "fig" and "ms" are the figure number and marker size.

if nargin == 2
    ms = 3;
elseif ms == 0
    ms = 3;
end

if nargin < 4
    figure(fig)
    plot3(P(:,1),P(:,2),P(:,3),'.b','Markersize',ms)
    axis equal

elseif nargin == 4
    I = vertcat(Bal{:});
    figure(fig)
    plot3(P(I,1),P(I,2),P(I,3),'.b','Markersize',ms)
    axis equal

else
    if iscell(Sub)
        S = vertcat(Sub{:});
        Sub = vertcat(S{:});
        I = vertcat(Bal{Sub});
        figure(fig)
        plot3(P(I,1),P(I,2),P(I,3),'.b','Markersize',ms)
        axis equal
    else
        I = vertcat(Bal{Sub});
        figure(fig)
        plot3(P(I,1),P(I,2),P(I,3),'.b','Markersize',ms)
        axis equal
    end
end