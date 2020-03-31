function plot_branch_structure(P,Bal,Segs,BSeg,BChi,fig,ms,gen,branum)

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
col = repmat(col,[1000,1]);

if iscell(Segs{1})
    n = size(Segs,1);
    segs = cell(n,1);
    for i = 1:n
        m = size(Segs{i},1);
        B = zeros(0);
        for j = 1:m
            s = Segs{i}(j);
            s = s{:};
            B = [B; s];
        end
        segs{i} = B;
    end
else
    segs = Segs;
end

m = size(P,1); 
forb = false(m,1);
if gen > 0
    B = unique(vertcat(Bal{segs{BSeg(branum)}}));
    figure(fig)
    plot3(P(B,1),P(B,2),P(B,3),'.','Color',col(1,:),'Markersize',ms)
    axis equal
    forb(B) = true;
    if gen > 1
        pause
        hold on
        c = BChi{branum};
        i = 2;
        while (i <= gen) && (~isempty(c))
            C = unique(vertcat(Bal{unique(vertcat(segs{BSeg(c)}))}));
            I = forb(C);
            C = C(~I);
            figure(fig)
            plot3(P(C,1),P(C,2),P(C,3),'.','Color',col(i,:),'Markersize',ms)
            axis equal
            c = unique(vertcat(BChi{c}));
            i = i+1;
            forb(C) = true;
            if i <= gen
                pause
            end
        end
        hold off
    end
else
    B = unique(vertcat(Bal{segs{BSeg(branum)}}));
    figure(fig)
    plot3(P(B,1),P(B,2),P(B,3),'.','Color',col(1,:),'Markersize',ms)
    forb(B) = true;
    hold on
    c = BChi{branum};
    i = 2;
    while ~isempty(c)
        C = unique(vertcat(Bal{unique(vertcat(segs{BSeg(c)}))}));
        I = forb(C);
        C = C(~I);
        figure(fig)
        plot3(P(C,1),P(C,2),P(C,3),'.','Color',col(i,:),'Markersize',ms)
        c = unique(vertcat(BChi{c}));
        i = i+1;
        forb(C) = true;
    end
    hold off
    axis equal
end