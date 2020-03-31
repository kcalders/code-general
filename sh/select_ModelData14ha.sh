mkdir selection

cat trees14ha.txt | while read line
do
	location=$(find . -name $line -print -quit)
	mv $location selection/
done

