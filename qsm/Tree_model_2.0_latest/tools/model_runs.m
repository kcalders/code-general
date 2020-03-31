function Summary = model_runs(P,dmin0,rcov0,nmin0,dmin,rcov,nmin,lcyl,NoGround,string,N)

% N         Number of modelling runs
% Summary   Mean, minimum, maximum and standard deviations of tree data

Data = zeros(33,N);
for i = 1:N
    str = [string,'_',num2str(i)];
    Data(:,i) = qsm_tree(P,dmin0,rcov0,nmin0,dmin,rcov,nmin,lcyl,NoGround,str);
    pause(2)
end

% Define and round off Summary
Summary = [mean(Data,2) min(Data,[],2) max(Data,[],2) std(Data,[],2)...
    std(Data,[],2)./mean(Data,2)*100];
for i = 1:33
    for j = 1:5
        D = Summary(i,j);
        if D > 100
            D = round(D);
        elseif D > 10
            D = round(10*D)/10;
        elseif D > 1
            D = round(100*D)/100;
        else
            D = round(1000*D)/1000;
        end
        Summary(i,j) = D;
    end
end

str = ['results/Summary_',string,'.mat'];
save(str,'Summary','Data')