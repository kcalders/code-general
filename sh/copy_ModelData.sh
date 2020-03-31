mkdir notfound
mkdir trackers

cat trees.txt | while read line
do
	echo "Analysing " $line
	if [ ! -f ../1.4ha_model/ModelData/$line ] && [ ! -f trackers/$line ]
	then
		touch trackers/$line #tmp tracker
		echo "Finding file " $line
		location=$(find . -name $line -not -path "./trackers/*" -print -quit)
		if [[ $location ]]  # var is set and it is not empty
		then
			echo 'copying file' $location
			cp $location ../1.4ha_model/ModelData/
		else
			echo 'File does not exist'
			touch notfound/$line
		fi
#		rm trackers/$line
	else
		echo "File already copied"
		
	fi
done

