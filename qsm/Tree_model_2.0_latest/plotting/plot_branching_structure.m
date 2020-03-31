function plot_branching_structure(BoC,Rad,Len,Axe,Sta,fig,ne,alp)
grid off
BO = max(BoC(:,2));
for i = 0:BO
    I = BoC(:,2) == i;
    if i >= 1
        hold on
    end
    plot_cylinder_model(Rad(I),Len(I),Axe(I,:),Sta(I,:),fig,ne,alp,i)
end
hold off
axis equal