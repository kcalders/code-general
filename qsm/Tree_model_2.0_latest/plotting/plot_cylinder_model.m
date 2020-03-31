function plot_cylinder_model(Rad,Len,Axe,Sta,fig,ne,alp,C)

% Plots the cylinders defined by Rad, Len, Axe, Sta.
% "fig" is figure number, "ne" is the number of facets in the cyliders, 
% "alp" is the alpha value, and "C" is the color.
grid off
if nargin == 5
    ne = 10;
    alp = 1;
elseif nargin == 6
    alp = 1;
end

if nargin <= 7
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
else
    m = max(Rad);
    figure(fig)
    draw(Rad(1),Len(1),Axe(1,:),Sta(1,:),1,ne,C)
    hold on
    for j = 2:length(Rad)
        n = ceil(Rad(j)/m*ne);
        n = max(n,6);
        draw(Rad(j),Len(j),Axe(j,:),Sta(j,:),1,n,C)
    end
    axis equal
    hold off
    alpha(alp)
end
