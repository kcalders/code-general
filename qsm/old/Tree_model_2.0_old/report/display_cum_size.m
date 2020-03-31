
load('results/data.mat')

h = plot_dist(CSDist(:,1),1,'Cumulative volume by size',strS,'Volume');
%str = ['results/cum_vol_size_abs_',string];
%saveas(h,str,'epsc')
pause(0.01)
h = plot_dist(CSDist(:,3),2,'Cumulative volume by size',strS,'Percent of total');
%str = ['results/cum_vol_size_abs_',string];
%saveas(h,str,'epsc')
pause(0.01)
h = plot_dist(CSDist(:,2),3,'Cumulative length by size',strS,'Volume');
%str = ['results/cum_len_size_abs_',string];
%saveas(h,str,'epsc')
pause(0.01)
h = plot_dist(CSDist(:,4),4,'Cumulative length by size',strS,'Percent of total');
%str = ['results/cum_len_size_abs_',string];
%saveas(h,str,'epsc')
pause(0.01)