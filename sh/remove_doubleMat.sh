find * -name "*mat.mat" > in.txt

cat in.txt | while read line
do
	newline=$(echo $line | rev| cut -c 5- | rev)
	mv $line $newline
done


rm in.txt