
h = plot_dist(BSDist(:,1),1,'Branch volume by size',strX,'Volume (L)');
str = ['vol_size_abs_',string];
saveas(h,str,'epsc')
pause(0.1)
h = plot_dist(BSDist(:,3),2,'Branch volume by size',strX,'Percent of total');
str = ['vol_size_rel_',string];
saveas(h,str,'epsc')
pause(0.1)
h = plot_dist(BODist(:,1),3,'Branch volume by order','Branch order','Volume (L)');
str = ['vol_ord_abs_',string];
saveas(h,str,'epsc')
pause(0.1)
h = plot_dist(BODist(:,4),4,'Branch volume by order','Branch order','Percent of total');
str = ['vol_ord_rel_',string];
saveas(h,str,'epsc')
pause(0.1)