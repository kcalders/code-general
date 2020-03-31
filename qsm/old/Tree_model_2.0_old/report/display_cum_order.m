
load('results/data.mat')

h = plot_dist([CODist(:,1) CODist(:,2)],1,'Cumulative volume by order','Branch order','Volume');
%str = ['results/cum_vol_ord_abs_',string];
%saveas(h,str,'epsc')
pause(0.01)
h = plot_dist([CODist(:,1) CODist(:,4)],2,'Cumulative volume by order','Branch order','Percent of total');
%str = ['results/cum_vol_ord_abs_',string];
%saveas(h,str,'epsc')
pause(0.01)
h = plot_dist([CODist(:,1) CODist(:,3)],3,'Cumulative length by order','Branch order','Volume');
%str = ['results/cum_vol_ord_abs_',string];
%saveas(h,str,'epsc')
pause(0.01)
h = plot_dist([CODist(:,1) CODist(:,5)],4,'Cumulative length by order','Branch order','Percent of total');
%str = ['results/cum_vol_ord_abs_',string];
%saveas(h,str,'epsc')
pause(0.01)