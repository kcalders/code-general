
load('results/data.mat')

h = plot_dist(BODist(:,3),1,'Number of branches by order','Branch order','Number');
%str = ['results/fre_ord_abs_',string];
%saveas(h,str,'epsc')
pause(0.01)
h = plot_dist(BODist(:,6),2,'Number of branches by order','Branch order','Percent of total');
%str = ['results/fre_ord_rel_',string];
%saveas(h,str,'epsc')
pause(0.01)