#!/usr/bin/env python

import optparse
import scipy.io
import sys
import os
import numpy as np

# main program
def main(cmdargs):
	# Check input file exists
	if not os.path.exists(cmdargs.inFile):
		print "%s does not exist locally" % cmdargs.inFile
   		sys.exit()
	else:
	        print "List of QSMs in file %s" % cmdargs.inFile

	f_out = open('QSMs_summary.csv','w')
	f_out.write('Tree_ID,Volume_m3,Branch_Length_[m],Height_m,DBH_cm,DBH_triangulation_cm\n')
	

	trees = open(cmdargs.inFile, "r")
	line = trees.readline().strip()
	while line: # process individual trees here
		print "Processing QSM %s" % line
		data = scipy.io.loadmat(line)
		blength= np.sum(data['BLen'])
		vol=data['TreeData'][0][0]/1000
		height=np.round(data['TreeData'][3][0],2)
		dbh=np.round(data['TreeData'][9][0],2)
		dbh_tr=np.round(data['TreeData'][10][0],2)
		#write results
		f_out.write(str(line.split('-')[0].split('_')[2]) + ',' + str(vol)+ ',' + str(np.round(blength,2))  +\
		',' + str(height) + ',' + str(dbh)+ ',' + str(dbh_tr)+'\n')
		line = trees.readline().strip()
	f_out.close()



# Command arguments
class CmdArgs:
  def __init__(self):
    p = optparse.OptionParser()
    p.add_option("-i","--inFile", dest="inFile", default=None, help="Input text file with names of individual QSM *mat models")    
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
