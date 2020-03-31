#!/usr/bin/env python

import os
import glob
import shutil

# fnames = glob.glob("stem_*.pcd")
fnames = glob.glob("stem*.pcd")
for i in xrange(len(fnames)):
	tmp1 = fnames[i].split("_")
	tmp2 = tmp1[1].split(".")
	tmp3 = tmp2[0]
 	cyl_name = "cylinder_" + tmp3 + ".pcd"
 	clus_name = "cluster_" + tmp3 + ".pcd"
 	run_name = "/Users/kcalders/code/cloudcompare/CloudCompare/CloudCompare/CloudCompare.app/Contents/MacOS/CloudCompare "+fnames[i] + " " + cyl_name +  " " + clus_name
	print run_name
	os.system(run_name)
#	os.remove(clus_name)
# 	os.remove(cyl_name)
#	shutil.move(fnames[i], "../../cleaned_summer_stems/"+fnames[i])

