function h = plot_dist(D,fig,strtit,strx,stry)

if size(D,2) == 1
    figure(fig)
    h = bar(D);
    t = title(strtit);
    x = xlabel(strx);
    y = ylabel(stry);
    set(gca,'fontsize',12)
    set(gca,'FontWeight','bold')
    set(t,'fontsize',12)
    set(t,'FontWeight','bold')
    set(x,'fontsize',12)
    set(x,'FontWeight','bold')
    set(y,'fontsize',12)
    set(y,'FontWeight','bold')
else
    figure(fig)
    h = bar(D(:,1),D(:,2));
    t = title(strtit);
    x = xlabel(strx);
    y = ylabel(stry);
    set(gca,'fontsize',12)
    set(gca,'FontWeight','bold')
    set(t,'fontsize',12)
    set(t,'FontWeight','bold')
    set(x,'fontsize',12)
    set(x,'FontWeight','bold')
    set(y,'fontsize',12)
    set(y,'FontWeight','bold')
end
grid on