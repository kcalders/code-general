
load('results/data3_Tree4_M.mat')
strR = 'M';

n = size(Sum,2);
np = max(size(dmin1)); % number of parameters
ns = round(n/np/N/5); % number of scans
scans = [3; 4; 6; 9];

rad = 41;
% Define diameter-volume distribution
vol = zeros(rad,ns*np*N);
t = 0;
for i = 1:ns % number of scans
    for j = 1:np % number of parameters
        for k = 1:N % number of models
            string = [strT,num2str(scans(i)),'_dmin_',num2str(1000*dmin1(j)),'_',num2str(k)];
            str = ['results/cyl_data_',string,'.mat'];
            load(str)
            R = CylData(:,1);
            L = CylData(:,2);
            t = t+1;
            for a = 1:rad
                if a <= rad-1
                    I = R <= a/200;
                    J = R > (a-1)/200;
                    I = I&J;
                else
                    I = R > 0.2;
                end
                vol(a,t) = 1000*pi*sum(R(I).^2.*L(I));
            end
        end
    end
end

% Determine means and std:s
m = zeros(33,ns*np);
s = zeros(33,ns*np);
v = zeros(rad,ns*np);
n = 5*np;
for i = 1:ns
    for j = 1:np
        m(:,(i-1)*np+j) = Sum(:,(i-1)*n+(j-1)*5+1);
        s(:,(i-1)*np+j) = Sum(:,(i-1)*n+(j-1)*5+4);
        v(:,(i-1)*np+j) = mean(vol(:,(i-1)*n+(j-1)*5+1:(i-1)*n+j*5),2);
    end
end

leg = num2str(dmin1);
TotV = zeros(ns,np);
TruV = zeros(ns,np);
BraV = zeros(ns,np);
TruL = zeros(ns,np);
BraL = zeros(ns,np);
TotVE = zeros(ns,np);
TruVE = zeros(ns,np);
BraVE = zeros(ns,np);
TruLE = zeros(ns,np);
BraLE = zeros(ns,np);
for i = 1:ns
    TotV(i,:) = m(1,(i-1)*np+1:i*np);
    TruV(i,:) = m(2,(i-1)*np+1:i*np);
    BraV(i,:) = m(3,(i-1)*np+1:i*np);
    TruL(i,:) = m(5,(i-1)*np+1:i*np);
    BraL(i,:) = m(6,(i-1)*np+1:i*np);
    TotVE(i,:) = s(1,(i-1)*np+1:i*np);
    TruVE(i,:) = s(2,(i-1)*np+1:i*np);
    BraVE(i,:) = s(3,(i-1)*np+1:i*np);
    TruLE(i,:) = s(5,(i-1)*np+1:i*np);
    BraLE(i,:) = s(6,(i-1)*np+1:i*np);
end

X = repmat(1000*dmin1',[4 1]);
%string = ['_tree2_',strR];

str = ['Total volume, ',strR,'-res'];
h = plot2d(X,TotV,1,str,'dmin (mm)','volume (L)',['3'; '4'; '6'; '9'],TotVE);
% str = ['total_volume_scans',string];
% saveas(h,str,'fig')
% saveas(h,str,'epsc')
% saveas(h,str,'png')

str = ['Trunk volume, ',strR,'-res'];
h = plot2d(X,TruV,2,str,'dmin (mm)','volume (L)',['3'; '4'; '6'; '9'],TruVE);
% str = ['trunk_volume_scans',string];
% saveas(h,str,'fig')
% saveas(h,str,'epsc')
% saveas(h,str,'png')

str = ['Branch volume, ',strR,'-res'];
h = plot2d(X,BraV,3,str,'dmin (mm)','volume (L)',['3'; '4'; '6'; '9'],BraVE);
% str = ['branch_volume_scans',string];
% saveas(h,str,'fig')
% saveas(h,str,'epsc')
% saveas(h,str,'png')

str = ['Trunk length, ',strR,'-res'];
h = plot2d(X,TruL,4,str,'dmin (mm)','length (m)',['3'; '4'; '6'; '9'],TruLE);
% str = ['trunk_length_scans',string];
% saveas(h,str,'fig')
% saveas(h,str,'epsc')
% saveas(h,str,'png')

str = ['Branch length, ',strR,'-res'];
h = plot2d(X,BraL,5,str,'dmin (mm)','length (m)',['3'; '4'; '6'; '9'],BraLE);
% str = ['branch_length_scans',string];
% saveas(h,str,'fig')
% saveas(h,str,'epsc')
% saveas(h,str,'png')


V = v;
P = zeros(ns,rad);
Q = P;
for j = 1:rad
    for i = 1:ns
        P(i,j) = sum(V(j:rad,(i-1)*np+1));
        Q(i,j) = sum(V(j:rad,i*np));
    end
end
R = zeros(np,rad);
S = R;
for j = 1:rad
    for i = 1:np
        R(i,j) = sum(V(j:rad,i));
        S(i,j) = sum(V(j:rad,(ns-1)*np+i));
    end
end

V = V';
C = cumsum(V,2);

str = ['Volume distribution, 3',strR];
h = plot2d([],V(1:np,:),6,str,'diameter (cm)','volume (L)',leg);
% str = ['volume_distri_dmin_3',string];
% saveas(h,str,'fig')
% saveas(h,str,'epsc')
% saveas(h,str,'png')

str = ['Cumulative volume, 3',strR];
h = plot2d([],C(1:np,:),7,str,'diameter (cm)','volume (L)',leg);
% str = ['cumulative_volume_dmin_3',string];
% saveas(h,str,'fig')
% saveas(h,str,'epsc')
% saveas(h,str,'png')

str = ['Difference in cumulative volume, 3',strR];
D = mat_vec_subtraction(C,C(1,:));
h = plot2d([],D(1:np,:),8,str,'diameter (cm)','volume (L)',leg);
% str = ['diff_cum_vol_dmin_3',string];
% saveas(h,str,'fig')
% saveas(h,str,'epsc')
% saveas(h,str,'png')

str = ['Volume distribution, 9',strR];
h = plot2d([],V((ns-1)*np+1:ns*np,:),9,str,'diameter (cm)','volume (L)',leg);
% str = ['volume_distri_dmin_9',string];
% saveas(h,str,'fig')
% saveas(h,str,'epsc')
% saveas(h,str,'png')

str = ['Cumulative volume, 9',strR];
h = plot2d([],C((ns-1)*np+1:ns*np,:),10,str,'diameter (cm)','volume (L)',leg);
% str = ['cumulative_volume_dmin_9',string];
% saveas(h,str,'fig')
% saveas(h,str,'epsc')
% saveas(h,str,'png')

str = ['Difference in cumulative volume, 9',strR];
D = mat_vec_subtraction(C,C((ns-1)*np+1,:));
h = plot2d([],D((ns-1)*np+1:ns*np,:),11,str,'diameter (cm)','volume (L)',leg);
% str = ['diff_cum_vol_dmin_9',string];
% saveas(h,str,'fig')
% saveas(h,str,'epsc')
% saveas(h,str,'png')

str = ['Over x-diameter part volumes, dmin = smallest, ',strR];
h = plot2d(repmat(1:rad,[4 1]),P,12,str,'diameter (cm)','volume (L)',['3'; '4';'6'; '9']);
% str = ['over_x_diam_dmin_small_scans',string];
% saveas(h,str,'fig')
% saveas(h,str,'epsc')
% saveas(h,str,'png')

str = ['Over x-diameter part volumes, dmin = largest, ',strR];
h = plot2d(repmat(1:rad,[4 1]),Q,13,str,'diameter (cm)','volume (L)',['3'; '4';'6'; '9']);
% str = ['over_x_diam_dmin_large_scans',string];
% saveas(h,str,'fig')
% saveas(h,str,'epsc')
% saveas(h,str,'png')

str = ['Over x-diameter part volumes, #scans = 3, ',strR];
h = plot2d(repmat(1:rad,[np 1]),R,14,str,'diameter (cm)','volume (L)',leg);
% str = ['over_x_diam_scan_3_dmin',string];
% saveas(h,str,'fig')
% saveas(h,str,'epsc')
% saveas(h,str,'png')

str = ['Over x-diameter part volumes, #scans = 9, ',strR];
h = plot2d(repmat(1:rad,[np 1]),S,15,str,'diameter (cm)','volume (L)',leg);
% str = ['over_x_diam_scan_9_dmin',string];
% saveas(h,str,'fig')
% saveas(h,str,'epsc')
% saveas(h,str,'png')

str = ['Differences over x-diam parts, dmin = smallest, ',strR];
D = mat_vec_subtraction(P,P(1,:));
h = plot2d(repmat(1:rad,[4 1]),D,16,str,'diameter (cm)','volume (L)',['3'; '4';'6'; '9']);
% str = ['diff_over_x_diam_dmin_small_scans',string];
% saveas(h,str,'fig')
% saveas(h,str,'epsc')
% saveas(h,str,'png')

str = ['Differences over x-diam parts, dmin = largest, ',strR];
D = mat_vec_subtraction(Q,Q(1,:));
h = plot2d(repmat(1:rad,[4 1]),D,17,str,'diameter (cm)','volume (L)',['3'; '4';'6'; '9']);
% str = ['diff_over_x_diam_dmin_large_scans',string];
% saveas(h,str,'fig')
% saveas(h,str,'epsc')
% saveas(h,str,'png')

str = ['Differences over x-diam parts, #scans = 3, ',strR];
D = mat_vec_subtraction(R,R(1,:));
h = plot2d(repmat(1:rad,[np 1]),D,18,str,'diameter (cm)','volume (L)',leg);
% str = ['diff_over_x_diam_scan_3_dmin',string];
% saveas(h,str,'fig')
% saveas(h,str,'epsc')
% saveas(h,str,'png')

str = ['Differences over x-diam parts, #scans = 9, ',strR];
D = mat_vec_subtraction(S,S(1,:));
h = plot2d(repmat(1:rad,[np 1]),D,19,str,'diameter (cm)','volume (L)',leg);
% str = ['diff_over_x_diam_scan_9_dmin',string];
% saveas(h,str,'fig')
% saveas(h,str,'epsc')
% saveas(h,str,'png')
