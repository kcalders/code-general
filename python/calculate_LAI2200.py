#!/usr/bin/env python
"""
Author: Kim Calders, kim.calders@npl.co.uk
August 2015

run as ./calculate_LAI2200.py

"""

import sys
import optparse
import csv
import bisect as b
import numpy as np
import time, datetime
import math
import glob

	
def timestamp2second(t):
	"""
	timestamp is hh:mm:ss in string format
	"""
	x=time.strptime(t,'%H:%M:%S')
	sec=datetime.timedelta(hours=x.tm_hour,minutes=x.tm_min,seconds=x.tm_sec).total_seconds()
	return sec
	
def nearest(s,ts,date):
    # Given a presorted list of timestamps:  s = sorted(index) [list], ts = timestamp to compare to (string)
    # from http://stackoverflow.com/questions/8162379/python-locating-the-closest-timestamp
    i = b.bisect_left(s,ts)
    if abs(timestamp2second(ts)-timestamp2second(s[i])) > abs(timestamp2second(ts)-timestamp2second(s[i-1])):
    	i=i-1
#     print "matched Below time %s with Above time %s on %s" % (ts,s[i],date)
    return i
    	
def calcualte_Pgap(belowfile):
	"""
	Calculate Pgap value for a single measurement. Read in Below file and time-match with Above reading
	Pgap defined here a B/A
	"""
	below=csv.reader(open(belowfile),delimiter='\t')
	Breadings=([])
	tmp=([])
	for row in below:
		tmp.append(row)
 	val_index=tmp.index(['### Observations'])+1
	for i in range(val_index,len(tmp)-1): #readings start at index 31
		if tmp[i][0]=='B':
			Breadings.append(tmp[i])
	Breadings=np.asarray(Breadings)

	date=Breadings[:][0,2].split(" ")[0]
	Btimestamps=([])
	for i in range(0,len(Breadings)):
		Btimestamps.append(Breadings[:][i,2].split(' ')[1])
	
	abovefile=date[2:8]+"_lai2200.txt"
	above=csv.reader(open(abovefile), delimiter='\t')
	Areadings=([])
	tmp=([])
	for row in above:
		tmp.append(row)
	for i in range(36,len(tmp)): #readings start at index 36
		if tmp[i][0]=='A':
			Areadings.append(tmp[i])
	Areadings=np.asarray(Areadings)
	Atimestamps=([])
	for i in range(0,len(Areadings)):
		Atimestamps.append(Areadings[:][i,2].split(' ')[1])
	
	idxA=([])
	for i in range(0,len(Btimestamps)):
		idxA.append(nearest(Atimestamps,Btimestamps[i],date))
	# calculate Pgap
	Pgap=np.zeros(shape=(len(idxA),5),dtype=float) #5 rings, most of the time 4 B reading, sometimes 5
	for i in range(0,len(Breadings)):
		for j in range(0,5):
			Pgap[i,j]=float(Breadings[i,j+4])/float(Areadings[idxA[i]][j+4])
	#CALCULATE MEAN B/A FOR EACH RING. WE DID THOSE N MEASUREMENTS TO APPROX THE LIGHT AT THAT SPECIFIC TIME, AVERAGE BEFORE CALCULATIONS
 	Pgap_psm_avg = np.mean(Pgap,axis=0)
 	Pgap_psm_sd = np.std(Pgap,axis=0)
 	Pgap_psm_sd_pct=Pgap_psm_sd/Pgap_psm_avg*100
	np.savetxt("output/Pgap_%s_all.csv" %(belowfile.split('.')[0]), Pgap, header='R1,R2,R3,R4,R5', delimiter=",")
	np.savetxt("output/Pgap_%s_avg.csv" %(belowfile.split('.')[0]), [Pgap_psm_avg], header='R1,R2,R3,R4,R5', delimiter=",")
	np.savetxt("output/Pgap_%s_sd.csv" %(belowfile.split('.')[0]), [Pgap_psm_sd], header='R1,R2,R3,R4,R5', delimiter=",")
  	np.savetxt("output/Pgap_%s_sd_pct.csv" %(belowfile.split('.')[0]), [Pgap_psm_sd_pct], header='R1,R2,R3,R4,R5', delimiter=",")

	
# 	return(Pgap)

	
def calculatePAIe(p,s):
	"""
	Calculates the PAIe for each VALERI plot
	PAIe = clumping index X PAI 
	See Ryu et. al (AFM, 2010) on the correct estimation of effictive leaf area index.
	
	Take ln of average Pgap (average over 12 VALERI samples)
		
	dist1=1.008;dist2=1.087;dist3=1.270;dist4=1.662;dist5=2.670 #values from LICOR2200 manual
	w1=0.041;w2=0.131;w3=0.201;w4=0.290;w5=0.337 #values from LICOR2200 manual

	In this script p=plot, s=subplot (i.e. one VALERI plot)
	"""
	print "Processing PAIe VALERI plot P%iS%i" % (p,s)
	
#	weight=[0.041,0.131,0.201,0.290,0.337]
# 	distance=[1.008,1.087,1.270,1.662,2.670]
	dist1=1.008;dist2=1.087;dist3=1.270;dist4=1.662;dist5=2.670 #values from LICOR2200 manual
# 	w1=0.041;w2=0.131;w3=0.201;w4=0.290;w5=0.337 #values from LICOR2200 manual
	w1=0.0461;w2=0.1431;w3=0.2085;w4=0.2655;w5=0.3368 #correct weights

	#AVERAGE PGAPS
	file_list = glob.glob("output/Pgap_P%iS%iM*_avg.csv" %(p,s))
#  	print file_list
  	Pgap_ps = ([])
  	for file_path in file_list:
  		Pgap_ps.append(np.genfromtxt(file_path, delimiter=',', skip_header=1))
	Pgap_ps = np.asarray(Pgap_ps)
 	Pgap_ps_avg = np.mean(Pgap_ps,axis=0)
  	Pgap_ps_sd = np.std(Pgap_ps,axis=0)	
 	np.savetxt("output/Pgap_P%iS%i.csv" %(p,s), Pgap_ps, header='R1,R2,R3,R4,R5', delimiter=",", fmt="%s")
 	np.savetxt("output/Pgap_P%iS%i_avg.csv" %(p,s), Pgap_ps_avg, header='R1,R2,R3,R4,R5', delimiter=",", fmt="%s")
 	np.savetxt("output/Pgap_P%iS%i_sd.csv" %(p,s), Pgap_ps_sd, header='R1,R2,R3,R4,R5', delimiter=",", fmt="%s")

 	
 	#CALCULATE PAIe
  	PAIe=2*(-math.log(Pgap_ps_avg[0])/dist1*w1-math.log(Pgap_ps_avg[1])/dist2*w2-math.log(Pgap_ps_avg[2])/dist3*w3-math.log(Pgap_ps_avg[3])/dist4*w4-math.log(Pgap_ps_avg[4])/dist5*w5)
	return PAIe

 	
def calculatePAI(p,s):
	"""
	Calculates the PAIe for each VALERI plot
	PAIe = clumping index X PAI 
	See Ryu et. al (AFM, 2010) on the correct estimation of effictive leaf area index.
	
	Take average of individual ln Pgap (average over 12 VALERI samples)
		
	dist1=1.008;dist2=1.087;dist3=1.270;dist4=1.662;dist5=2.670 #values from LICOR2200 manual
	w1=0.041;w2=0.131;w3=0.201;w4=0.290;w5=0.337 #values from LICOR2200 manual

	In this script p=plot, s=subplot (i.e. one VALERI plot)
	"""
	print "Processing PAI VALERI plot P%iS%i" % (p,s)
	
#	weight=[0.041,0.131,0.201,0.290,0.337]
# 	distance=[1.008,1.087,1.270,1.662,2.670]
	dist1=1.008;dist2=1.087;dist3=1.270;dist4=1.662;dist5=2.670 #values from LICOR2200 manual
# 	w1=0.041;w2=0.131;w3=0.201;w4=0.290;w5=0.337 #values from LICOR2200 manual
	w1=0.0461;w2=0.1431;w3=0.2085;w4=0.2655;w5=0.3368 #correct weights


	file_list = glob.glob("output/Pgap_P%iS%iM*_avg.csv" %(p,s))
#  	print file_list
  	Pgap_ps = ([])
  	for file_path in file_list:
  		Pgap_ps.append(np.genfromtxt(file_path, delimiter=',', skip_header=1))
#	Pgap_ps = np.asarray(Pgap_ps)
# 	Pgap_ps_avg = np.mean(Pgap_ps,axis=0)
#  	Pgap_ps_sd = np.std(Pgap_ps,axis=0)	
# 	np.savetxt("output/Pgap_P%iS%i.csv" %(p,s), Pgap_ps, header='R1,R2,R3,R4,R5', delimiter=",", fmt="%s")
# 	np.savetxt("output/Pgap_P%iS%i_avg.csv" %(p,s), Pgap_ps_avg, header='R1,R2,R3,R4,R5', delimiter=",", fmt="%s")
# 	np.savetxt("output/Pgap_P%iS%i_sd.csv" %(p,s), Pgap_ps_sd, header='R1,R2,R3,R4,R5', delimiter=",", fmt="%s")
	
	PAI_m = ([])
	for i in range(0,12):	
 		#CALCULATE PAI per measurement
  		PAI_m.append(2*(-math.log(Pgap_ps[i][0])/dist1*w1-math.log(Pgap_ps[i][1])/dist2*w2-math.log(Pgap_ps[i][2])/dist3*w3-math.log(Pgap_ps[i][3])/dist4*w4-math.log(Pgap_ps[i][4])/dist5*w5))
	PAI=np.mean(PAI_m)
	return PAI


def main(cmdargs):	
 	#CALCULATE ALL Pgap_psm values
 	for p in (1,2,5,6,7,8):
		for s in range(1,26):
			for m in range(1,13):
				belowfile="P%iS%iM%i.TXT" %(p,s,m)
				print "Processing below reading file %s" % belowfile
				Pgap=calcualte_Pgap(belowfile)
 	
 	#CALCULATE ALL PAI values per VALERI plot; this is effective PAI
 	VALERI=([])
 	PAIe_vals=([])
 	PAI_vals=([])
  	for p in (1,2,5,6,7,8):
 		for s in range(1,26):	
			PAIe_vals.append([calculatePAIe(p,s)])
			PAI_vals.append([calculatePAI(p,s)])
			VALERI.append(["P%iS%i" %(p,s)])
 	PAI_VALERI=np.hstack([VALERI,PAIe_vals,PAI_vals])
 	np.savetxt("SUMMARY_PAI_VALERI.csv", PAI_VALERI, header='VALERI_PLOT,PAIe,PAI', fmt="%s", delimiter=",")


# Command arguments
class CmdArgs:
  def __init__(self):
    p = optparse.OptionParser()
    (options, args) = p.parse_args()
    self.__dict__.update(options.__dict__)
    
#    if self.inFile is None:
#        p.print_help()
#        print "Input filename must be set."
#        sys.exit()
#	if self.outFile is None:
#        p.print_help()
#        print "Output filename must be set."
#        sys.exit()

# Run the script
if __name__ == "__main__":
    cmdargs = CmdArgs()
    main(cmdargs)
