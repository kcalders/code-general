#!/usr/bin/env python

import sys
import open3d as o3d
import os.path

for txt in sys.argv[1:]:
    cloud=o3d.io.read_point_cloud(txt)
    outname=txt.split('.')[0]+'.ply'
    if not os.path.isfile(outname): #check if file is already converted
    	o3d.io.write_point_cloud(outname,cloud)
