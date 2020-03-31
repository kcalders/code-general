function plot_tree_structure(P,Bal,Seg,SChi,fig,ms,gen,segnum)

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

if iscell(Seg{1})
    n = max(size(Seg));
    Segs = cell(n,1);
    for i = 1:n
        m = size(Seg{i},1);
        S = zeros(0);
        for j = 1:m
            s = Seg{i}(j);
            s = s{:};
            S = [S; s];
            %pause
        end
        Segs{i} = S;
    end
else
    Segs = Seg;
end

%ms = 1;
if gen > 0
    S = unique(vertcat(Bal{Segs{segnum}}));
    figure(fig)
    plot3(P(S,1),P(S,2),P(S,3),'.','Color',col(1,:),'Markersize',ms)
    axis equal
    forb = S;
    if gen > 1
        pause
        hold on
        c = SChi{segnum};
        i = 2;
        while (i <= gen) && (~isempty(c))
            C = unique(vertcat(Bal{unique(vertcat(Segs{c}))}));
            C = setdiff(C,forb);
            figure(fig)
            plot3(P(C,1),P(C,2),P(C,3),'.','Color',col(i,:),'Markersize',ms)
            axis equal
            c = unique(vertcat(SChi{c}));
            i = i+1;
            forb = union(forb,C);
            if i <= gen
                pause
            end
        end
        hold off
    end
    i-1
else
    S = unique(vertcat(Bal{Segs{1}}));
    figure(fig)
    plot3(P(S,1),P(S,2),P(S,3),'.','Color',col(1,:),'Markersize',ms)
    hold on
    n = max(size(Segs));
    segs = zeros(n,1);
    segs(1) = 1;
    segnum = 1;
    i = 2;
    while nnz(segs) < n
        c = SChi{segnum};
        while ~isempty(c)
            C = unique(vertcat(Bal{unique(vertcat(Segs{c}))}));
            segs(c) = 1;
            plot3(P(C,1),P(C,2),P(C,3),'.','Color',col(i,:),'Markersize',ms)
            c = unique(vertcat(SChi{c}));
            I = segs(c);
            if sum(I) == length(c)
                c = zeros(0);
            end
            i = i+1;
        end
        if nnz(segs) < n
            segnum = 2;
            while segs(segnum) == 1
                segnum = segnum+1;
            end
            segs(segnum) = 1;
            S = unique(vertcat(Bal{Segs{segnum}}));
            plot3(P(S,1),P(S,2),P(S,3),'.','Color',col(1,:),'Markersize',ms)
            i = 2;
        end
    end
    axis equal
    hold off
end