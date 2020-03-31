
load('results/data2.mat')

disp(['Tree: ' str])
disp(datestr(clock))
disp(['dmin0 = ',num2str(dmin0),', rcov0 = ',num2str(rcov0),', nmin0 = ',num2str(nmin0)])
disp(['dmin = ',num2str(dmin),', rcov = ',num2str(rcov),...
    ', nmin = ',num2str(nmin),', lcyl = ',num2str(lcyl)])
disp(['Number of models: ',num2str(N)])
pause(0.01)
