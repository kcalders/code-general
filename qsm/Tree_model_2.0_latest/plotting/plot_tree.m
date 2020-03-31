function plot_tree(P,Bal,Segs,SChi,Rad,Sta,Len,Axe,CiS,fig,gen,segnum,ms,alp)

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
col = repmat(col,[100,1]);

ne = 30;
R = max(Rad);

if iscell(Segs{1})
    n = max(size(Segs));
    segs = cell(n,1);
    for i = 1:n
        m = size(Segs{i},1);
        S = zeros(0);
        for j = 1:m
            s = Segs{i}(j);
            s = s{:};
            S = [S; s];
        end
        segs{i} = S;
    end
else
    segs = Segs;
end

if gen > 0
    % plot the first segment
    S = unique(vertcat(Bal{segs{segnum}}));
    figure(fig)
    plot3(P(S,1),P(S,2),P(S,3),'.','Color',col(1,:),'Markersize',ms)
    axis equal
    hold on
    forb = S;
    
    % plot the cylinders of the first segment
    C = CiS{segnum};
    if ~isempty(C)
        figure(fig)
        for j = 1:length(C)
            a = max(6,ceil(ne*Rad(C(j))/R));
            draw(Rad(C(j)),Len(C(j)),Axe(C(j),:),Sta(C(j),:),1,a)
        end
        axis equal
        alpha(alp)
    end
    
    if gen > 1
        pause
        c = SChi{segnum};
        i = 2;
        while (i <= gen) && (~isempty(c))
            C = unique(vertcat(Bal{unique(vertcat(segs{c}))}));
            C = setdiff(C,forb);
            figure(fig)
            plot3(P(C,1),P(C,2),P(C,3),'.','Color',col(i,:),'Markersize',ms)
            axis equal
            forb = union(forb,C);
            
            C = vertcat(CiS{c});
            if ~isempty(C)
                figure(fig)
                a = max(6,ceil(ne*Rad(C(1))/R));
                draw(Rad(C(1)),Len(C(1)),Axe(C(1),:),Sta(C(1),:),1,a)
                for j = 2:length(C)
                    a = max(6,ceil(ne*Rad(C(j))/R));
                    draw(Rad(C(j)),Len(C(j)),Axe(C(j),:),Sta(C(j),:),1,a)
                end
                axis equal
                alpha(alp)
            end
            
            c = unique(vertcat(SChi{c}));
            i = i+1;
            pause
        end
        hold off
    end
end