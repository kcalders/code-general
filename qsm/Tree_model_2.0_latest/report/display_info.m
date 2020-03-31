
load('results/data.mat')

disp(['Tree: ' string])
disp(datestr(clock))
disp(['Modeling time: ' num2str(Tmin),' min ',num2str(Tsec),' sec'])
disp(['dmin0 = ',num2str(dmin0),', rcov0 = ',num2str(rcov0),', nmin0 = ',num2str(nmin0)]);
disp(['dmin = ',num2str(dmin),', rcov = ',num2str(rcov),...
    ', nmin = ',num2str(nmin),', lcyl = ',num2str(lcyl)]);
pause(0.01)
