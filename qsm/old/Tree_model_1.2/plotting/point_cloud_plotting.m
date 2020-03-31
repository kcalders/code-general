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

if size(P,2) == 7
    I = P(:,5) > P(:,6);
    J = P(:,5) > P(:,7);
    I1 = I&J;
    I = P(:,6) > P(:,5);
    J = P(:,6) > P(:,7);
    I2 = I&J;
    I = P(:,7) > P(:,5);
    J = P(:,7) > P(:,6);
    I3 = I&J;
    disp([nnz(I1) nnz(I2) nnz(I3)])
    
    figure(fig)
    plot3(P(I1,1),P(I1,2),P(I1,3),'.r','Markersize',ms)
    hold on
    plot3(P(I2,1),P(I2,2),P(I2,3),'.g','Markersize',ms)
    plot3(P(I3,1),P(I3,2),P(I3,3),'.b','Markersize',ms)
    axis equal
    hold off
    
    I = P(:,5) > P(:,6);
    J = P(:,7) > P(:,6);
    I = I|J;
    disp([nnz(I) nnz(~I)])
    
    figure(fig+1)
    plot3(P(I,1),P(I,2),P(I,3),'.r','Markersize',ms)
    hold on
    plot3(P(~I,1),P(~I,2),P(~I,3),'.g','Markersize',ms)
    axis equal
    hold off
else
    
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
end