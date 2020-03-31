#!/usr/bin/env python

import optparse
import scipy.io
import sys
import os
import numpy as np
import glob


def averageQSMs(tree,f_out):
	print "Analysing %s" %tree
	file_list = glob.glob("%s*" %(tree))
  	height = ([])
  	blength = ([])
  	vol = ([])
  	dbh = ([])
#  	dbh_tr = ([])
  	x = ([])
  	y = ([])
  	for file in file_list:
#  		print "opening %s" %file
  		data = scipy.io.loadmat(file)
  		i=0
		while sum(data['Len'][0:i+1]) < 1.3:
			i = i+1
		x.append(data['Sta'][i][0])
		y.append(data['Sta'][i][1])
  		height.append(np.round(data['TreeData'][3][0],2))
  		blength.append(np.sum(data['BLen']))
		vol.append(data['TreeData'][0][0]/1000)
		dbh.append(np.round(data['TreeData'][9][0],2))
#		dbh_tr.append(np.round(data['TreeData'][10][0],2))
	qsm_iterations=len(file_list)
	height = np.real(height)
  	blength = np.real(blength)
  	vol = np.real(vol)
  	dbh = np.real(dbh)
	x_avg = np.mean(x,axis=0)
	y_avg = np.mean(y,axis=0)
	vol_avg = np.mean(vol,axis=0)
	vol_sd = np.std(vol, ddof=1,axis=0) #the same as excell STDEV
	blength_avg = np.mean(blength,axis=0)
	blength_sd = np.std(blength, ddof=1,axis=0) #the same as excell STDEV
	height_avg = np.mean(height,axis=0)
	height_sd = np.std(height, ddof=1,axis=0) #the same as excell STDEV
	dbh_avg = np.mean(dbh,axis=0)
	dbh_sd = np.std(dbh, ddof=1,axis=0) #the same as excell STDEV
#	dbh_tr_avg = np.mean(dbh_tr,axis=0)
#	dbh_tr_sd = np.std(dbh_tr, ddof=1,axis=0) #the same as excell STDEV
	
	params=str(file.split('-',1)[1].rsplit('-',1)[0])
	
#	f_out.write(tree + ',' + str(qsm_iterations) + ',' + str(x_avg) + ','+ str(y_avg) + ',' + str(vol_avg)+ ',' + str(vol_sd)+ ',' + str(np.round(blength_avg,4))  +\
#	',' + str(np.round(blength_sd,4)) +	',' + str(height_avg) +	',' + str(height_sd)+ ',' + str(dbh_avg)+ ',' + str(dbh_sd)+ ',' + str(dbh_tr_avg)+ ',' + str(dbh_tr_sd)+ ',' + params +'\n')
	
#	return tree,x_avg, x_avg, vol_avg,vol_sd,blength_avg,blength_sd,height_avg,height_sd,dbh_avg,dbh_sd,dbh_tr_avg,dbh_tr_sd

	f_out.write(tree + ',' + str(qsm_iterations) + ',' + str(x_avg) + ','+ str(y_avg) + ',' + str(vol_avg)+ ',' + str(vol_sd)+ ',' + str(np.round(blength_avg,4))  +\
	',' + str(np.round(blength_sd,4)) +	',' + str(height_avg) +	',' + str(height_sd)+ ',' + str(dbh_avg)+ ',' + str(dbh_sd)+ ',' + params +'\n')
	
	return tree,x_avg, x_avg, vol_avg,vol_sd,blength_avg,blength_sd,height_avg,height_sd,dbh_avg,dbh_sd


 	


# main program
def main(cmdargs):
	tree_list = [x.split("-")[0] for x in glob.glob("*mat")]	
	myset = set(tree_list)
#	print type(list(myset)[0])
	tree_list=sorted(list(myset))
	f_out = open('QSMs_summary.csv','w')
#	f_out.write('Tree_ID,#iterations,x,y,Volume_mean_m3,Volume_sd_m3,Branch_Length_mean_[m],Branch_Length_sd_[m],Height_mean_m,Height_sd_m,DBH_mean_cm,DBH_sd_cm,DBH_triangulation_mean_cm,DBH_triangulation_sd_cm,qsm_parameters\n')
	f_out.write('Tree_ID,#iterations,x,y,Volume_mean_m3,Volume_sd_m3,Branch_Length_mean_[m],Branch_Length_sd_[m],Height_mean_m,Height_sd_m,DBH_mean_cm,DBH_sd_cm,qsm_parameters\n')
	for tree in tree_list:
		averageQSMs(tree,f_out)
	f_out.close()




# Command arguments
class CmdArgs:
  def __init__(self):
    p = optparse.OptionParser()
#    p.add_option("-i","--inFile", dest="inFile", default=None, help="Input text file with names of individual QSM *mat models")    
    (options, args) = p.parse_args()
    self.__dict__.update(options.__dict__)
    
#    if self.inFile is None:
#        p.print_help()
#        print "Input text filename must be set."
#        sys.exit()

# Run the script
if __name__ == "__main__":
    cmdargs = CmdArgs()
    main(cmdargs)
