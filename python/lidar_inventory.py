#!/usr/bin/env python
"""
Author: Kim Calders, kim.calders@wur.nl
Feb 2014

read in textfile with listed individual trees (x,y,z in col 1,2,3)
calcualte height and DBH + xy location for each tree

output: TreeID, x, y, DBH, H

"""

import sys
import optparse
import os
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import optimize
from matplotlib.patches import Circle

def PointsInCircum(xc,yc,r,n):
    return [(xc+math.cos(2*math.pi/n*x)*r,yc+math.sin(2*math.pi/n*x)*r) for x in xrange(0,n+1)]

def load_tree(f): 
  fin = open(f,'r')
  num_lines = sum(1 for line in fin)
  fin.close()
  points = np.zeros((num_lines,3),dtype=float)
  # now read
  fin = open(f,'r')
  i = 0
  for line in fin.xreadlines():
    rec = line.split(' ')
    x = float(rec[0])
    y = float(rec[1])
    z = float(rec[2])
    points[i] = (x,y,z)
    i += 1
  fin.close()
  return points

def DBH_plot(points,center_coordinates,radius,tree):
  fig = plt.figure()
  ax = plt.subplot(111,aspect='equal') #the aspect keeps the ratio 1:1 between the axis
  circle = Circle((center_coordinates[0],center_coordinates[1]),radius,color='red',fc='w')
  circle.set_linewidth(5)
  ax.add_patch(circle)
  plt.title('Tree ' + str(tree))
  plt.plot(points[0],points[1],'o')
  plt.plot(center_coordinates[0],center_coordinates[1],'o')
  plt.xlim(center_coordinates[0]-0.5,center_coordinates[0]+0.5)
  plt.ylim(center_coordinates[1]-0.5,center_coordinates[1]+0.5)
  plt.savefig('plot/DBH_circle_fitting_tree'+str(tree)+'.png')
  plt.close()
  # plt.show()



# main program
def main(cmdargs):
	# Check input file exists
	if not os.path.exists(cmdargs.inFile):
		print "%s does not exist locally" % cmdargs.inFile
   		sys.exit()
	else:
	        print "List of trees in file %s" % cmdargs.inFile

	f_out = open('lidar_inventory.csv','w')
	f_out.write('File_Name,Tree_ID,Center_x_lidar_m,Center_y_lidar_m,DBH_lidar_cm,hgt_lidar_m\n')
	

	trees = open(cmdargs.inFile, "r")
	line = trees.readline().strip()
	while line: # process individual trees here
		print "Processing tree %s" % line
		pts=load_tree(line)
		#tree height		
		hgt=max(pts[:,2])-min(pts[:,2])
		if hgt > 1.3:
			#dbh From http://wiki.scipy.org/Cookbook/Least_Squares_Circle, algebraic approach
			dbh_pts=pts[np.where( (pts[:,2] > min(pts[:,2])+1.27) & (pts[:,2] < min(pts[:,2])+1.33) )] #limits according to Tansey et al 2009
			fdbh="dbh/dbh_"+line
			np.savetxt(fdbh, dbh_pts, delimiter=',') 
	#		plt.scatter(dbh_pts[:,0],dbh_pts[:,1])
	#		plt.show()
			x=dbh_pts[:,0]
			y=dbh_pts[:,1]
			z=dbh_pts[:,2]
			x_m = x.mean()
			y_m = y.mean()
			# calculation of the reduced coordinates
			u = x - x_m
			v = y - y_m
			# linear system defining the center (uc, vc) in reduced coordinates:
			#    Suu * uc +  Suv * vc = (Suuu + Suvv)/2
			#    Suv * uc +  Svv * vc = (Suuv + Svvv)/2
			Suv  = sum(u*v)
			Suu  = sum(u**2)
			Svv  = sum(v**2)
			Suuv = sum(u**2 * v)
			Suvv = sum(u * v**2)
			Suuu = sum(u**3)
			Svvv = sum(v**3)
			# Solving the linear system
			A = np.array([ [ Suu, Suv ], [Suv, Svv]])
			B = np.array([ Suuu + Suvv, Svvv + Suuv ])/2.0
			uc, vc = np.linalg.solve(A, B)
			xc_1 = x_m + uc
			yc_1 = y_m + vc
			# Calcul des distances au centre (xc_1, yc_1)
			Ri_1     = np.sqrt((x-xc_1)**2 + (y-yc_1)**2)
			r_final      = np.mean(Ri_1) # this is the final radius!
			residu_1 = sum((Ri_1-r_final)**2)

			#plot results
			DBH_plot((x,y),(xc_1,yc_1),r_final,line.split('.')[0])
			#write results
			f_out.write(str(line.split('.')[0]) + ',' + str(line.split('.')[0].split('_')[2]) + ',' + str(np.round(xc_1,2)) + ',' + str(np.round(yc_1,2)) + ',' +\
			str(np.round(r_final*2*100,2)) + ',' + str(np.round(hgt,2)) + '\n')
			#export circlepoints to plot circle
			circlepoints=PointsInCircum(xc_1,yc_1,r_final,200)
			fcircle="circle/circle_"+line
			np.savetxt(fcircle, circlepoints, delimiter=',') 

		else:  #this happened when using stems instead of whole trees
			#dbh From http://wiki.scipy.org/Cookbook/Least_Squares_Circle, algebraic approach
			lower_lim=hgt-0.08 #we still want a 6cm slice, starting from 2cm below the height of stem until 8cm below: 8-2=6
			upper_lim=hgt-0.02
			dbh_pts=pts[np.where( (pts[:,2] > min(pts[:,2])+lower_lim) & (pts[:,2] < min(pts[:,2])+upper_lim) )]
			fdbh="dbh/dbh_"+line
			np.savetxt(fdbh, dbh_pts, delimiter=',') 
	#		plt.scatter(dbh_pts[:,0],dbh_pts[:,1])
	#		plt.show()
			x=dbh_pts[:,0]
			y=dbh_pts[:,1]
			z=dbh_pts[:,2]
			x_m = x.mean()
			y_m = y.mean()
			# calculation of the reduced coordinates
			u = x - x_m
			v = y - y_m
			# linear system defining the center (uc, vc) in reduced coordinates:
			#    Suu * uc +  Suv * vc = (Suuu + Suvv)/2
			#    Suv * uc +  Svv * vc = (Suuv + Svvv)/2
			Suv  = sum(u*v)
			Suu  = sum(u**2)
			Svv  = sum(v**2)
			Suuv = sum(u**2 * v)
			Suvv = sum(u * v**2)
			Suuu = sum(u**3)
			Svvv = sum(v**3)
			# Solving the linear system
			A = np.array([ [ Suu, Suv ], [Suv, Svv]])
			B = np.array([ Suuu + Suvv, Svvv + Suuv ])/2.0
			uc, vc = np.linalg.solve(A, B)
			xc_1 = x_m + uc
			yc_1 = y_m + vc
			# Calcul des distances au centre (xc_1, yc_1)
			Ri_1     = np.sqrt((x-xc_1)**2 + (y-yc_1)**2)
			r_final      = np.mean(Ri_1) # this is the final radius!
			residu_1 = sum((Ri_1-r_final)**2)

			#plot results
			DBH_plot((x,y),(xc_1,yc_1),r_final,line.split('.')[0])
			#write results
			f_out.write(str(line.split('.')[0].split('_')[2]) + ',' + str(np.round(xc_1,2)) + ',' + str(np.round(yc_1,2)) + ',' +\
			str(np.round(r_final*2*100,2)) + ',' + str(np.round(hgt,2)) + '\n')
			#export circlepoints to plot circle
			circlepoints=PointsInCircum(xc_1,yc_1,r_final,200)
			fcircle="circle/circle_"+line
			np.savetxt(fcircle, circlepoints, delimiter=',') 

		line = trees.readline().strip()
	f_out.close()

# Command arguments
class CmdArgs:
  def __init__(self):
    p = optparse.OptionParser()
    p.add_option("-i","--inFile", dest="inFile", default=None, help="Input text file with names of individual trees textfiles")    
    (options, args) = p.parse_args()
    self.__dict__.update(options.__dict__)
    
    if self.inFile is None:
        p.print_help()
        print "Input text filename must be set."
        sys.exit()

# Run the script
if __name__ == "__main__":
    cmdargs = CmdArgs()
    main(cmdargs)



