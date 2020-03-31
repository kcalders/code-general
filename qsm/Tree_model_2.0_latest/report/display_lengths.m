
load('results/data.mat')

h = plot_dist(BSDist(:,2),1,'Branch length by size',strX,'Length (m)');
%str = ['results/len_size_abs_',string];
%saveas(h,str,'epsc')
pause(0.01)
h = plot_dist(BSDist(:,4),2,'Branch length by size',strX,'Percent of total');
%str = ['results/len_size_rel_',string];
%saveas(h,str,'epsc')
pause(0.01)
h = plot_dist(BODist(:,2),3,'Branch length by order','Branch order','Length (m)');
%str = ['results/len_ord_abs_',string];
%saveas(h,str,'epsc')
pause(0.01)
h = plot_dist(BODist(:,5),4,'Branch length by order','Branch order','Percent of total');
%str = ['results/len_ord_rel_',string];
%saveas(h,str,'epsc')
pause(0.01)