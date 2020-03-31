
load('results/data.mat')

h = plot_dist(Taper(:,1:2),1,'Trunk taper','Length (m)','Diameter (cm)');
%str = ['results/taper_',string];
%saveas(h,str,'epsc')
pause(0.01)
h = plot_dist([Taper(:,1) Taper(:,3)],2,'Cumulative trunk volume','Length (m)','Volume (L)');
%str = ['results/cum_vol_trunk_abs_',string];
%saveas(h,str,'epsc')
pause(0.01)
h = plot_dist([Taper(:,1) Taper(:,4)],3,'Cumulative trunk volume','Length (m)','Percent of total');
%str = ['results/cum_vol_trunk_rel_',string];
%saveas(h,str,'epsc')
pause(0.01)