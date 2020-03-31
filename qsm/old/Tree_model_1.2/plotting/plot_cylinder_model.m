function plot_cylinder_model(Rad,Len,Axe,Sta,fig,ne,alp,P)

% Plots the cylinders defined by Rad, Len, Axe, Sta.
% "fig" is figure number, "ne" is the number of facets in the cyliders, 
% "alp" is the alpha value, and "P" is the point cloud.

if nargin == 5
    ne = 10;
    alp = 1;
elseif nargin == 6
    alp = 1;
end

m = max(Rad);
figure(fig)
draw(Rad(1),Len(1),Axe(1,:),Sta(1,:),1,ne)
hold on
for j = 2:length(Rad)
    n = ceil(Rad(j)/m*ne);
    n = max(n,6);
    draw(Rad(j),Len(j),Axe(j,:),Sta(j,:),1,n)
end
axis equal
hold off
alpha(alp)

if nargin == 8
    hold on
    plot3(P(:,1),P(:,2),P(:,3),'.r','Markersize',8)
    axis equal
    hold off
end