find * -name "*" -and -not -name "*txt" > in.txt

cat in.txt | while read line
do
	mv $line $line.mat
done


rm in.txt