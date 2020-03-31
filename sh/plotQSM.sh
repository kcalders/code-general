#!/bin/bash

# in.txt is the file with all the paths of the files that have to processed
mkdir trackers
mkdir images

ls *txt > infiles.dat

cat infiles.dat | while read xyz
do
	echo processing $xyz
	name=$(echo $xyz | cut -d .  -f1)
	if [ ! -f  images/$name\_clouds_trim.png ] && [ ! -f trackers/$name ]
	then
		touch trackers/$name #temp tracker
		plotTrees.py -c $xyz -out images/$name -az 30 -el 30
		convert -trim images/$name\_clouds.png images/$name\_clouds_trim.png
		rm images/$name\_clouds.png
		rm trackers/$name
	fi

	
	ls opt_1mod/${name}*mat > tmp  #list all the generated models
	cat tmp | while read cyl
	do
		echo plotting model $cyl
		namemod=$(echo $cyl |  cut -c 10- | rev | cut -c 5- | rev)
		echo $namemod
		if [ ! -f  images/$namemod\_models_trim.png ] && [ ! -f trackers/$namemod ]
		then
			touch trackers/$namemod	#temp tracker
			plotTrees.py -m $cyl -out images/$namemod -az 30 -el 30
			convert -trim images/$namemod\_models.png images/$namemod\_models_trim.png
			rm images/$namemod\_models.png
			plotTrees.py -m $cyl -c $xyz -ol -out images/$namemod -az 30 -el 30
			convert -trim images/$namemod\_cloudsmodels.png images/$namemod\_cloudsmodels_trim.png
			rm images/$namemod\_cloudsmodels.png
			rm trackers/$namemod
		fi
	done
done


