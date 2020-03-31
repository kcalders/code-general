
disp('Tree attributes:')
disp(['Total volume = ' num2str(TreeData(1)) ' liters'])
disp(['Trunk volume = ' num2str(TreeData(2)) ' liters'])
disp(['Branch volume = ' num2str(TreeData(3)) ' liters'])
disp(['Total height = ' num2str(TreeData(4)) ' meters'])
disp(['Trunk length = ' num2str(TreeData(5)) ' meters'])
disp(['Branch length = ' num2str(TreeData(6)) ' meters'])
disp(['Number of branches = ' num2str(TreeData(7))])
disp(['Maximum branch order = ' num2str(TreeData(8))])
disp(['Total cylinder area = ' num2str(TreeData(9)) ' square meters'])
disp(['Dbh (cylinder) = ' num2str(TreeData(10)) ' centimeters'])
disp(['Dbh (triangulation) = ' num2str(TreeData(11)) ' centimeters'])
disp(['Trunk volume (cylinders) = ' num2str(TreeData(12)) ' liters'])
disp(['Trunk volume (triangulation) = ' num2str(TreeData(13)) ' liters'])
disp(['Trunk length (cylinders) = ' num2str(TreeData(14)) ' meters'])
disp(['Trunk length (triangulation) = ' num2str(TreeData(15)) ' meters'])
disp('-----')
disp('Branch order data:')
S = ['1st'; '2nd'; '3rd'; '4th'; '5th'; '6th'];
for i = 1:min(6,BO)
    str = ['Number of ',S(i),'-order branches = ' num2str(TreeData(15+i))];
    disp(str)
end
for i = 1:min(6,BO)
    str = ['Volume of ',S(i),'-order branches = ' num2str(TreeData(21+i)),' liters'];
    disp(str)
end
for i = 1:min(6,BO)
    str = ['Length of ',S(i),'-order branches = ' num2str(TreeData(27+i)),' meters'];
    disp(str)
end