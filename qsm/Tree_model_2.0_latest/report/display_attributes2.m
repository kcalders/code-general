
load('results/data2.mat')


disp('Tree attributes:  Mean,  Minimum,  Maximum')
disp(['Total volume = ' num2str(Summary(1,1:3)) ' L'])
disp(['Trunk volume = ' num2str(Summary(2,1:3)) ' L'])
disp(['Branch volume = ' num2str(Summary(3,1:3)) ' L'])
disp(['Total height = ' num2str(Summary(4,1:3)) ' m'])
disp(['Trunk length = ' num2str(Summary(5,1:3)) ' m'])
disp(['Branch length = ' num2str(Summary(6,1:3)) ' m'])
disp(['Number of branches = ' num2str(Summary(7,1:3))])
disp(['Maximum branch order = ' num2str(Summary(8,1:3))])
disp(['Total cylinder area = ' num2str(Summary(9,1:3)) ' m^2'])
disp(['Dbh (cylinder) = ' num2str(Summary(10,1:3)) ' cm'])
disp(['Dbh (triangulation) = ' num2str(Summary(11,1:3)) ' cm'])
disp(['Trunk volume (cylinders) = ' num2str(Summary(12,1:3)) ' L'])
disp(['Trunk volume (triangulation) = ' num2str(Summary(13,1:3)) ' L'])
disp(['Trunk length (cylinders) = ' num2str(Summary(14,1:3)) ' m'])
disp(['Trunk length (triangulation) = ' num2str(Summary(15,1:3)) ' m'])
disp('-----')
disp('Branch order data:')
S = ['1st'; '2nd'; '3rd'; '4th'; '5th'; '6th'];
for i = 1:min(6,BO)
    str = ['Number of ',S(i),'-order branches = ' num2str(Summary(15+i,1:3))];
    disp(str)
end
for i = 1:min(6,BO)
    str = ['Volume of ',S(i),'-order branches = ' num2str(Summary(21+i,1:3)),' L'];
    disp(str)
end
for i = 1:min(6,BO)
    str = ['Length of ',S(i),'-order branches = ' num2str(Summary(27+i,1:3)),' m'];
    disp(str)
end

disp('')
disp('------------')
disp('')

disp('Tree attributes:  STD,  STD/Mean (%)')
disp(['Total volume = ' num2str(Summary(1,4:5))])
disp(['Trunk volume = ' num2str(Summary(2,4:5))])
disp(['Branch volume = ' num2str(Summary(3,4:5))])
disp(['Total height = ' num2str(Summary(4,4:5))])
disp(['Trunk length = ' num2str(Summary(5,4:5))])
disp(['Branch length = ' num2str(Summary(6,4:5))])
disp(['Number of branches = ' num2str(Summary(7,4:5))])
disp(['Maximum branch order = ' num2str(Summary(8,4:5))])
disp(['Total cylinder area = ' num2str(Summary(9,4:5))])
disp(['Dbh (cylinder) = ' num2str(Summary(10,4:5))])
disp(['Dbh (triangulation) = ' num2str(Summary(11,4:5))])
disp(['Trunk volume (cylinders) = ' num2str(Summary(12,4:5))])
disp(['Trunk volume (triangulation) = ' num2str(Summary(13,4:5))])
disp(['Trunk length (cylinders) = ' num2str(Summary(14,4:5))])
disp(['Trunk length (triangulation) = ' num2str(Summary(15,4:5))])
disp('-----')
disp('Branch order data:')
S = ['1st'; '2nd'; '3rd'; '4th'; '5th'; '6th'];
for i = 1:min(6,BO)
    str = ['Number of ',S(i),'-order branches = ' num2str(Summary(15+i,4:5))];
    disp(str)
end
for i = 1:min(6,BO)
    str = ['Volume of ',S(i),'-order branches = ' num2str(Summary(21+i,4:5))];
    disp(str)
end
for i = 1:min(6,BO)
    str = ['Length of ',S(i),'-order branches = ' num2str(Summary(27+i,4:5))];
    disp(str)
end
pause(0.01)