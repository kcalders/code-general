function [FeatureValues,FeatureNames] = compute_features(QSMs)

Nmodels = max(size(QSMs));
F = zeros(20000,Nmodels);
FN = cell(20000,1);

% TotalVolume
% TrunkVolume
% BranchVolume
% TreeHeight
% TrunkLength
% BranchLength
% TotalLength
% NumberBranches    Total number of branches
% MaxBranchOrder 
% TrunkArea 
% BranchArea 
% TotalArea 
% DBHqsm        From the cylinder of the QSM at the right heigth
% DBHcyl        From the cylinder fitted to the section 1.1-1.5m
% CrownDiamAve
% CrownDiamMax
% CrownAreaConv
% CrownAreaAlpha
% CrownBaseHeight
% CrownLength
% CrownRatio
% CrownVolumeConv
% CrownVolumeAlpha
% location      (x,y,z)-coordinates of the base of the tree
% StemTaper     Stem taper function/curve from the QSM
% VerticalProfile
% spreads
% VolCylDia     Distribution of the total volume in diameter classes
% AreCylDia     Distribution of the total area in diameter classes
% LenCylDia     Distribution of the total length in diameter classes
% VolCylHei     Distribution of the total volume in height classes
% AreCylHei     Distribution of the total area in height classes
% LenCylHei     Distribution of the total length in height classes
% VolCylAzi     Distribution of the total volume in azimuth angle classes
% AreCylAzi     Distribution of the total area in azimuth angle classes
% LenCylAzi     Distribution of the total length in azimuth angle classes
% VolCylZen     Distribution of the total volume in zenith angle classes
% AreCylZen     Distribution of the total area in zenith angle classes
% LenCylZen     Distribution of the total length in zenith angle classes
% VolBranchOrd     Branch volume per branching order
% AreBranchOrd       Branch area per branching order
% LenBranchOrd     Branch length per branching order
% NumBranchOrd     Number of branches per branching order
% VolBranchDia
% VolBranch1Dia
% AreBranchDia
% AreBranch1Dia
% LenBranchDia
% LenBranch1Dia
% NumBranchDia
% NumBranch1Dia
% VolBranchAng
% VolBranch1Ang
% AreBranchAng
% AreBranch1Ang
% LenBranchAng
% LenBranch1Ang
% NumBranchAng
% NumBranch1Ang
% VolBranchAzi 
% VolBranch1Azi
% AreBranchAzi
% AreBranch1Azi
% LenBranchAzi
% LenBranch1Azi
% NumBranchAzi
% NumBranch1Azi
% VolBranchHei
% VolBranch1Hei
% AreBranchHei
% AreBranch1Hei
% LenBranchHei
% LenBranch1Hei
% NumBranchHei
% NumBranch1Hei
% VolBranchZen 
% VolBranch1Zen 
% AreBranchZen 
% AreBranch1Zen 
% LenBranchZen 
% LenBranch1Zen 
% NumBranchZen 
% NumBranch1Zen 

NameB = fieldnames(QSMs(1).branch);
mb = size(NameB,1);

NameT = fieldnames(QSMs(1).treedata);
mt = size(NameT,1);
for i = 1:mt
    if strcmp(NameT{i},'CrownVolumeAlpha')
        nf = i;
    end
    if strcmp(NameT{i},'VolCylDia')
        ncs = i;
    end
    if strcmp(NameT{i},'LenCylZen')
        nce = i;
    end
    if strcmp(NameT{i},'VolBranchOrd')
        nbs = i;
    end
    if strcmp(NameT{i},'NumBranch1Zen')
        nbe = i;
    end
end

% Shedding ratios
SR = [0 1 2 3 4 0 0 0 0 1 1 1 2 2 3]/5;
ER = [1 2 3 4 5 2 3 4 5 3 4 5 4 5 5]/5;

% relative height/diameter/zenith of bottom 5, 10, 15, etc.
% volume/area/length
RHN = (5:5:95); % for names
RH = RHN/100; % for computations

% Cylinder height/diameter classes (Sc (m/cm) - Ec (m/cm)) and branch
% orders (Sb - Eb)
Sc = [1 1 1 1 1  3 3 3 3  5 5 5  7 7  9]/10;
Ec = [2 4 6 8 10 4 6 8 10 6 8 10 8 10 10]/10;
Sb = [1 2 3 4 5 1 1 1 1 2 2 2 3 3 4];
Eb = [1 2 3 4 5 2 3 4 5 3 4 5 4 5 5];

% Uniform and step distributions
Su = [0 0 0 0 0 1 2 3 4 1 2 3 1 1 2]; % relative start (0/5, 1/5, 2/5, etc)
Eu = [1 2 3 4 5 5 5 5 5 4 4 4 3 2 3]; % relative ends (1/5, 2/5, 5/5, etc)

% normal distributions (mean,std)
Mn = [1 1 2 2 3 3]; % means (relative, quadrants, i.e. 1/4 2/4, etc)
Sn = [1 2 1 2 1 2]; % standard deviations (relative, 1/5 and 2/5)
        
% Triangle distributions (start, mid/high point, end points):
% start points:
St = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 2 2 2 2 2 3];
% mid/high points
Mt = [0 0 0 0 1 1 1 1 2 2 2 3 3 4 1 1 1 2 2 2 3 3 4 2 2 3 3 4 4];
% end points
Et = [1 2 3 4 1 2 3 4 2 3 4 3 4 4 2 3 4 2 3 4 3 4 4 3 4 3 4 4 4];

%% Treedata features
for i = 1:Nmodels
    f = 0; % feature index
    t = QSMs(i).treedata; % treedata of the model
    if ~isempty(t)
    
    %% Markku's features from the 2017 RSE-paper
    % Stem branch angle (median branching angle of 1st-ord. branches)
    b = QSMs(i).branch;
    Ord1 = b.order == 1;
    nb = length(b.order);
    indb = (1:1:nb)';
    Ord1 = indb(Ord1);
    m1 = length(Ord1);
    f = f+1;
    F(f,i) = median(b.angle(Ord1));
    if i == 1
    FN{f} = 'StemBranchAngle';
    end
    % Stem branch cluster size (Mean number of 1branches inside 40cm height 
    % interval. Each branch can only belong to one interval).
    if m1 > 0
        h1 = b.height(Ord1);
        Int = zeros(m1,1);
        k = 10;
        j = 1;
        while j <= m1
            h0 = h1(j);
            jk = j;
            while jk <= m1 && h1(jk) < h0+0.4
                jk = jk+1;
            end
            jk = jk-1;
            k = k+1;
            Int(k) = jk-j+1;
            j = jk+1;
        end
        CS = mean(Int(1:k));
    else
        CS = 0;
    end
    f = f+1;
    F(f,i) = CS;
    if i == 1
    FN{f} = 'StemBranchClusterSize';
    end
    % Stem branch radius (Mean ratio between the 10 largest 1branches and 
    % stem radius at respective height).
    r1 = b.diameter(Ord1)/2;
    [r,I] = sort(r1,'descend');
    if m1 >= 10
        r = r(I(1:10));
        Ind = Ord1(I(1:10));
    else
        Ind = Ord1;
    end
    m = length(r);
    c = QSMs(i).cylinder;
    R = r;
    for j = 1:m
        a = c.branch == Ind(j) & c.PositionInBranch == 1;
        p = c.parent(a);
        R(j) = c.radius(p);
    end
    f = f+1;
    F(f,i) = mean(r./R);
    if i == 1
    FN{f} = 'StemBranchRadius';
    end
    % Stem branch length (Mean length of 1branches normalized by DBH)
    f = f+1;
    F(f,i) = mean(b.length(Ord1))/t.DBHcyl;
    if i == 1
    FN{f} = 'StemBranchLength';
    end  
    % Stem branch distance (Mean distance between 1branches from moving 
    % average with 1m window. If empty window, set value to half the
    % window width. Normilize with dbh).
    if m1 > 0
        D = zeros(m1^2,1);
        k = 1;
        for j = 1:m1
            d = abs(h1-h1(j));
            I = d <= 0.5;
            I(j) = false;
            d = d(I);
            m = length(d);
            if m > 0
                D(k:k+m-1) = d;
            else
                D(k) = 0.5;
                m = 1;
            end
            k = k+m;
        end
        k = k-m;
        D = mean(D(1:k));
    else
        D = 0;
    end
    f = f+1;
    F(f,i) = D/t.DBHcyl;
    if i == 1
    FN{f} = 'StemBranchDistance';
    end 
    % Crown start height (Height of first stem branch in crown relative to 
    % tree height)
    f = f+1;
    F(f,i) = t.CrownBaseHeight/t.TreeHeight;
    if i == 1
    FN{f} = 'CrownStartHeight';
    end 
    % Crown height (Crown length divided by the tree height)
    f = f+1;
    F(f,i) = t.CrownRatio;
    if i == 1
    FN{f} = 'CrownHeight';
    end 
    % Crown evenness (Crown cylinders divided into 8 angular bins. Ratio 
    % between extreme minimum heights in bins)
    nc = length(c.radius);
    indc = (1:1:nc)';
    bot = min(c.start(:,3));
    hcs = c.start(:,3)-bot;
    E = c.start+[c.length.*c.axis(:,1) c.length.*c.axis(:,2) c.length.*c.axis(:,3)];
    hce = E(:,3)-bot;
    I = (hcs >= t.CrownBaseHeight | hce >= t.CrownBaseHeight) & c.BranchOrder > 0;
    crown = indc(I);
    if ~isempty(crown)
        hs = hcs(crown);
        he = hce(crown);
        Cen = mean(E(crown,:),1);
        V = mat_vec_subtraction(E(crown,:),Cen);
        ang = atan2(V(:,2),V(:,1))+pi;
        H = zeros(8,1);
        for j = 1:8
            I = ang >= (j-1)*pi/4 & ang < j*pi/4;
            if any(I)
                H(j) = min(min(hs(I)),min(he(I)));
            end
        end
    end
    f = f+1;
    F(f,i) = min(H(H > 0))/max(H);
    if i == 1
    FN{f} = 'CrownEvenness';
    end 
    % Crown diameter/height (Ratio between crown diameter and height)
    f = f+1;
    F(f,i) = t.CrownDiamAve/t.CrownLength;
    if i == 1
    FN{f} = 'CrownDiameter/Height';
    end
    % DBH/height ratio (Ratio between DBH and tree height)
    f = f+1;
    F(f,i) = t.DBHcyl/t.TreeHeight;
    if i == 1
    FN{f} = 'DBH/TreeHeight';
    end
    % DBH/tree volume ratio (Ratio between DBH and total tree volume)
    f = f+1;
    F(f,i) = t.DBHcyl/t.TotalVolume;
    if i == 1
    FN{f} = 'DBH/TotalVolume';
    end
    % DBH/minimum tree radius (Ratio between DBH and the minimum of the 
    % vertical bins radius estimate)
    height= max(hce);
    CylVol = 1000*pi*c.radius.^2.*c.length;
    R = zeros(3,1);
    for j = 1:3
        I = hce >= (j-1)*height/3 & hce < j*height/3;
        L = indc(I);
        V = CylVol(L);
        stem = c.branch(L) == 1;
        stem = L(stem);
        if ~isempty(stem)
            Cen = mean(c.start(stem,:),1);
        end
        d = distances_to_line(E(L,:),[0 0 1],Cen);
        [d,I] = sort(d);
        v = cumsum(V(I))/sum(V);
        k = 1;
        while v(k) < 0.9
            k = k+1;
        end
        R(j) = 2*d(k);
    end
    f = f+1;
    F(f,i) = t.DBHcyl/min(R);
    if i == 1
    FN{f} = 'DBH/MinimumTreeRadius';
    end
    % Volume below 55% of height (Relative volume below 55%)
    I = hce/height < 0.55;
    f = f+1;
    F(f,i) = sum(CylVol(I))/t.TotalVolume;
    if i == 1
    FN{f} = 'VolumeBelow55%';
    end
    % Cylinder length/tree volume (Ratio between total length and total volume)
    f = f+1;
    F(f,i) = t.TotalLength/t.TotalVolume;
    if i == 1
    FN{f} = 'TotalLength/TotalVolume';
    end
    % Shedding ratio (Number of branches without children divided by the number
    % of branches in the bottom third)
    B = indb(b.height < t.TreeHeight/3 & b.order == 1);
    m = length(B);
    k = 0;
    for j = 1:m
        if ~any(b.parent == B(j))
            k = k+1;
        end
    end
    f = f+1;
    if m > 0
    F(f,i) = k/m;
    end
    if i == 1
    FN{f} = 'SheddingRatio';
    end
    
    
    % Shedding ratios (Number of branches without children divided by the number
    % of branches in different height layers)
    h = b.height;
    o = b.order;
    th = t.TreeHeight;
    for l = 1:15
        B = indb(h < ER(l)*th & h >= SR(l)*th & o == 1);
        m = length(B);
        k = 0;
        for j = 1:m
            if ~any(b.parent == B(j))
                k = k+1;
            end
        end
        f = f+1;
        if m > 0
        F(f,i) = k/m;
        end
        if i == 1
        FN{f} = ['SheddingRatio_',num2str(SR(l)),'_',num2str(ER(l))];
        end
    end
        
        
    %% Treedata and treedata divided by other treedata
    for j = 1:nf
        f = f+1;
        F(f,i) = t.(NameT{j});
        if i == 1
        FN{f} = NameT{j};
        end
        for k = 1:nf
            f = f+1;
            if t.(NameT{k}) ~= 0
                F(f,i) = t.(NameT{j})/t.(NameT{k});
            end
            if i == 1
            FN{f} = [NameT{j},'/',NameT{k}];
            end
        end
    end
    
    
    %% Tree attributes divided by treedata
    % Cylinder volumes, areas and lengths between certain diameter and 
    % height classes and branch orders divided by treedata
    for j = 1:nf
        a = t.(NameT{j});
        for l = ncs:ncs+5
            D = t.(NameT{l});
            if length(D) < 10
                D(10) = 0; 
            end
            b = length(D);
            for k = 1:length(Sc)
                f = f+1;
                F(f,i) = sum(D(ceil(Sc(k)*b):floor(Ec(k)*b)))/a; % Vol S-E Diam / a
                if i == 1
                FN{f} = [NameT{l},'_',num2str(Sc(k)-0.1),'_',num2str(Ec(k)),'/',NameT{j}];
                end
            end
        end
        for l = nbs:nbs+3
            D = t.(NameT{l});
            if length(D) < 5
                D(5) = 0;  
            end
            for k = 1:length(Sb)
                f = f+1;
                F(f,i) = sum(D(Sb(k):Eb(k)))/a; % #branch_k / a
                if i == 1
                FN{f} = [NameT{l},'_',num2str(Sb(k)),'_',num2str(Eb(k)),'/',NameT{j}];
                end
            end
        end
    end

    
    %% Tree segments (cylinder) distributions
    % Volume, area, length of cylinder as functions of dia, hei, azi, zen
    for j = ncs:nce
        d = t.(NameT{j}); % distribution
        if strcmp(NameT{j}(1:3),'Vol')
            a = t.TotalVolume;
        elseif strcmp(NameT{j}(1:3),'Are')
            a = t.TotalArea;
        elseif strcmp(NameT{j}(1:3),'Len')
            a = t.TotalLength;
        end
        dr = d/a; % relative distribution
        dc = cumsum(dr); % cumulative relative distribution
        
        % relative height/diameter/zenith of 5, 10, 15 etc % of bottom 
        % volume/area/length
        N = relative_height(dc,RH);
        for k = 1:length(RH)
            f = f+1;
            F(f,i) = N(k);
            if i == 1
            FN{f} = ['Rel_cyl_',NameT{j}(end-2:end),'_bottom_',...
                num2str(RHN(k)),'%_',NameT{j}(1:3)];
            end
        end
        
        % distribution comparisons:
        m = length(dr);
        % differences to the triangle distributions (start,mid,end)
        % start points:
        Sti = ceil(St/4*m);
        Mti = ceil(Mt/4*m);
        Eti = ceil(Et/4*m);
        for k = 1:length(St)
            di = abs(dr-triad(Sti(k),Mti(k),Eti(k),m));
            % mean and max difference to the comparison distribution:
            f = f+1;
            if ~isempty(di)
            F(f,i) = mean(di);
            end
            if i == 1
            FN{f} = [NameT{j},'_triad_',num2str(St(k)),'_',...
                num2str(Mt(k)),'_',num2str(Et(k)),'_mean'];
            end
            f = f+1;
            if ~isempty(di)
            F(f,i) = max(di);
            end
            if i == 1
            FN{f} = [NameT{j},'_triad_',num2str(St(k)),'_',...
                num2str(Mt(k)),'_',num2str(Et(k)),'_max'];
            end
            
%             if strcmp(NameT{j}(1:3),'Len')
%             g = figure(1);
%             bar((1:1:m),[dr; triad(Sti(k),Mti(k),Eti(k),m); di]',0.7)
%             grid on
%             title([NameT{j},' vs. triangle distribution. Mean difference: ',...
%                 num2str(F(f-1,i))])
%             legend('empirical','theoretical','difference')
%             saveas(g,FN{f-1},'png')
%             pause
%             end
        end
        
        % differences to the normal distributions (mean,std)
        Mni = ceil(Mn/4*m);
        Sni = ceil(Sn/5*m);
        for k = 1:length(Mn)
            di = abs(dr-normd(Mni(k),Sni(k),m));
            % mean and max difference to the comparison distribution:
            f = f+1;
            if ~isempty(di)
            F(f,i) = mean(di);
            end
            if i == 1
            FN{f} = [NameT{j},'_normd_',num2str(Mn(k)),'_',num2str(Sn(k)),'_mean'];
            end
            f = f+1;
            if ~isempty(di)
            F(f,i) = max(di);
            end
            if i == 1
            FN{f} = [NameT{j},'_normd_',num2str(Mn(k)),'_',num2str(Sn(k)),'_max'];
            end
        end
        
        % difference to the uniform and step distributions
        Sui = ceil(Su/5*m);
        Eui = ceil(Eu/5*m);
        for k = 1:length(Su)
            di = abs(dr-unid(Sui(k),Eui(k),m));
            % mean and max difference to the comparison distribution:
            f = f+1;
            if ~isempty(di)
            F(f,i) = mean(di);
            end
            if i == 1
            FN{f} = [NameT{j},'_unid_',num2str(Su(k)),'_',num2str(Eu(k)),'_mean'];
            end
            f = f+1;
            if ~isempty(di)
            F(f,i) = max(di);
            end
            if i == 1
            FN{f} = [NameT{j},'_unid_',num2str(Su(k)),'_',num2str(Eu(k)),'_max'];
            end
        end
    end
    
    
    %% Branch distributions
    % Volume, area, length, number (of 1st-order) of branches as functions 
    % of height, diameter, angle
    for j = nbs:nbe
        d = t.(NameT{j}); % distribution
        if strcmp(NameT{j}(1:3),'Vol')
            a = t.BranchVolume;
            if strcmp(NameT{j}(end-3),'1')
                a = t.VolBranchOrd(1);
            end
        elseif strcmp(NameT{j}(1:3),'Are')
            a = t.BranchArea;
            if strcmp(NameT{j}(end-3),'1')
                a = t.AreBranchOrd(1);
            end
        elseif strcmp(NameT{j}(1:3),'Len')
            a = t.BranchLength;
            if strcmp(NameT{j}(end-3),'1')
                a = t.LenBranchOrd(1);
            end
        elseif strcmp(NameT{j}(1:3),'Num')
            a = t.NumberBranches;
            if strcmp(NameT{j}(end-3),'1')
                a = t.NumBranchOrd(1);
            end
        end
        dr = d/a; % relative distribution
        dc = cumsum(dr); % cumulative relative distribution
        
        % relative height/diameter/zenith of 5, 10, 15 etc % of bottom 
        % volume/area/length
        N = relative_height(dc,RH);
        for k = 1:19
            f = f+1;
            F(f,i) = N(k);
            if i == 1
                if strcmp(NameT{j}(end-3),'h')
                    FN{f} = ['Rel_branch_',NameT{j}(end-2:end),...
                        '_bottom_',num2str(RHN(k)),'%_',NameT{j}(1:3)];
                else
                    FN{f} = ['Rel_branch1_',NameT{j}(end-2:end),...
                        '_bottom_',num2str(RHN(k)),'%_',NameT{j}(1:3)];
                end
            end
%             if strcmp('Rel_branch_Hei_bottom_20%_Are',FN{f})
%                 disp(F(f,i))
%             end
%             if strcmp('Rel_branch1_Hei_bottom_20%_Are',FN{f})
%                 disp(F(f,i))
%             end
        end
        
        % distribution comparisons:
        m = length(dr);
        % differences to the triangle distributions (start,mid,end)
        Sti = ceil(St/4*m);
        Mti = ceil(Mt/4*m);
        Eti = ceil(Et/4*m);
        for k = 1:length(St)
            di = abs(dr-triad(Sti(k),Mti(k),Eti(k),m));
            % mean and max difference to the comparison distribution:
            f = f+1;
            if ~isempty(di)
            F(f,i) = mean(di);
            end
            if i == 1
            FN{f} = [NameT{j},'_triad_',num2str(St(k)),...
                '_',num2str(Mt(k)),'_',num2str(Et(k)),'_mean'];
            end
            f = f+1;
            if ~isempty(di)
            F(f,i) = max(di);
            end
            if i == 1
            FN{f} = [NameT{j},'_triad_',num2str(St(k)),'_',...
                num2str(Mt(k)),'_',num2str(Et(k)),'_max'];
            end
        end
        
        % differences to the normal distributions (mean,std)
        Mni = ceil(Mn/4*m);
        Sni = ceil(Sn/5*m);
        for k = 1:length(Mn)
            di = abs(dr-normd(Mni(k),Sni(k),m));
            % mean and max difference to the comparison distribution:
            f = f+1;
            if ~isempty(di)
            F(f,i) = mean(di);
            end
            if i == 1
            FN{f} = [NameT{j},'_normd_',num2str(Mn(k)),'_',num2str(Sn(k)),'_mean'];
            end
            f = f+1;
            if ~isempty(di)
            F(f,i) = max(di);
            end
            if i == 1
            FN{f} = [NameT{j},'_normd_',num2str(Mn(k)),'_',num2str(Sn(k)),'_max'];
            end
        end
        
        % difference to the uniform and step distributions
        Sui = ceil(Su/5*m);
        Eui = ceil(Eu/5*m);
        for k = 1:length(Su)
            di = abs(dr-unid(Sui(k),Eui(k),m));
            % mean and max difference to the comparison distribution:
            f = f+1;
            if ~isempty(di)
            F(f,i) = mean(di);
            end
            if i == 1
            FN{f} = [NameT{j},'_unid_',num2str(Su(k)),'_',num2str(Eu(k)),'_mean'];
            end
            f = f+1;
            if ~isempty(di)
            F(f,i) = max(di);
            end
            if i == 1
            FN{f} = [NameT{j},'_unid_',num2str(Su(k)),'_',num2str(Eu(k)),'_max'];
            end
        end
    end
    
    
    %% Branch-data
    % Medians, averages, minimums, maximums and their ratios of branch-data
    b = QSMs(i).branch;
    I1 = b.order == 1;
    I2 = b.order == 2;
    I3 = b.order == 3;
    for j = [1 3:mb]
        d = double(b.(NameB{j}));
        for o = 0:3
            if o == 0
                a = d; % all branches
            elseif o == 1
                a = d(I1); % 1st-order branches
            elseif o == 2
                a = d(I2); % 2nd-order branches
            elseif o == 3
                a = d(I3); % 3rd-order branches
            end
            f = f+1;
            if ~isempty(a)
            F(f:f+6,i) = [median(a); mean(a); min(a); max(a); ...
                mean(a)/max(a); min(a)/max(a); min(a)/mean(a)];
            end
            if i == 1
                FN{f} = ['branch_or',num2str(o),'_',NameB{j},'_median'];
                FN{f+1} = ['branch_or',num2str(o),'_',NameB{j},'_mean'];
                FN{f+2} = ['branch_or',num2str(o),'_',NameB{j},'_min'];
                FN{f+3} = ['branch_or',num2str(o),'_',NameB{j},'_max'];
                FN{f+4} = ['branch_or',num2str(o),'_',NameB{j},'_mean/max'];
                FN{f+5} = ['branch_or',num2str(o),'_',NameB{j},'_min/max'];
                FN{f+6} = ['branch_or',num2str(o),'_',NameB{j},'_min/mean'];
            end
            f = f+6;
            for k = [1 3:mb]
                if k ~= j
                    c = double(b.(NameB{k}));
                    for or = 0:3
                        if or == 0
                            a = c;
                        elseif or == 1
                            a = c(I1);
                        elseif or == 2
                            a = c(I2);
                        elseif or == 3
                            a = c(I3);
                        end
                        f = f+1;
                        if ~isempty(a) && ~isempty(c)
                        F(f,i) = mean(c)/mean(a);
                        end
                        if i == 1
                        FN{f} = ['branch_or',num2str(or),'_',NameB{k},...
                            '_mean/','or',num2str(o),'_',NameB{j},'_mean'];
                        end
                    end
                end
            end
        end
    end

    
    %% Branch azimuth
    % distribution comparisons
    m = t.NumberBranches;
    a = t.NumBranchAzi/m; % relative branch-angle distribution
    di = abs(a-unid(0,36,36)); % difference to the 0-360 uniform distribution
    f = f+1;
    if ~isempty(di)
    F(f,i) = mean(di);   % mean difference to the uniform angle distribution
    end
    if i == 1
    FN{f} = 'NumBranchAzi_unid_0_360_mean';
    end
    f = f+1;
    if ~isempty(di)
    F(f,i) = max(di);   % maximum difference to the uniform angle distribution
    end
    if i == 1
    FN{f} = 'NumBranchAzi_unid_0_360_max';
    end
    
    m1 = t.NumBranchOrd(1);
    a = t.NumBranch1Azi/m1; % relative branch-angle distribution (1st-branches)
    di = abs(a-unid(0,36,36));   % difference to the 0-360 uniform distribution
    f = f+1;
    if ~isempty(di)
    F(f,i) = mean(di); % mean difference to the uniform angle distribution
    end
    if i == 1
    FN{f} = 'NumBranch1Azi_unid_0_360_mean';
    end
    f = f+1;
    if ~isempty(di)
    F(f,i) = max(di); % maximum difference to the uniform angle distribution
    end
    if i == 1
    FN{f} = 'NumBranch1Azi_unid_0_360_max';
    end
    
    end
end
f
FeatureValues = F(1:f,:);
FeatureNames = FN(1:f);

end % End of main function



function d = normd(m,v,n)
x = (1/n:1/n:1);
G = normpdf(x,m/n,v/n);
d = G/sum(G);
end

function d = triad(a,c,b,n)
x = (0.5:1:n-0.5);
T = zeros(1,n);
T(a+1:c) = 2*(x(a+1:c)-a)/(c-a)/(b-a);
T(c+1:b) = 2*(b-x(c+1:b))/(b-a)/(b-c);
d = T/sum(T);
end

function d = unid(a,b,n)
c = b-a;
d = [zeros(1,a) 1/c*ones(1,c) zeros(1,n-a-c)];
end

function N = relative_height(d,N)
h = length(d);
k = 1;
n = length(N);
for j = 1:h
    for i = k:n
        if k <= n && d(j) > N(k)
            N(k) = j;
            k = k+1;
        end
    end
end
N = N/h;
end
