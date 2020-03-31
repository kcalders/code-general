#!/usr/bin/env python

import os
import glob
import shutil

# fnames = glob.glob("stem_*.pcd")
fnames = glob.glob("wytham*.pcd")
for i in xrange(len(fnames)):
	tmp1 = fnames[i].split("_")
	tmp2 = tmp1[2].split(".")
	tmp3 = tmp2[0]
# 	cyl_name = "cylinder_" + tmp3 + ".pcd"
# 	clus_name = "cluster_" + tmp3 + ".pcd"
	slice_name = "slices_" + tmp3 + ".pcd"
# 	run_name = "/Users/kimcalders/code/cloudcompare/CloudCompare/CloudCompare/CloudCompare.app/Contents/MacOS/CloudCompare "+fnames[i] + " " + cyl_name +  " " + clus_name
	run_name = "/Users/kcalders/code/cloudcompare/CloudCompare/CloudCompare/CloudCompare.app/Contents/MacOS/CloudCompare "+fnames[i] + " " +  " " + slice_name
	print run_name
	os.system(run_name)
# 	os.remove(clus_name)
# 	os.remove(cyl_name)
 	shutil.move(fnames[i], "../../processed/"+fnames[i])
 	shutil.move(slice_name, "../../processed/"+slice_name)
#	shutil.move(fnames[i], "../cleaned/cluster_"+fnames[i])

