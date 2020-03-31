#!/usr/bin/env python
"""
This will take an *obj output from meshlab (during export, deselect all normals and colors etc.) and turn it into a librat *obj

The librat *.obj is also optimised with bounding boxes

In Meshlab: use filters > remeshing, simplification and reconstruction > quadric edge collapse decimation. Use percentage reduction 0.5 (i.e. 50%)

Kim Calders - kim.calders@npl.co.uk
Andrew Burt - a.burt.12@ucl.ac.uk
June 2016

"""

import sys
import argparse
import numpy as np

class mesh2librat(object):
	
	def __init__(self,file_in,file_out,file_out_bound,fpb,name):
		self.f1 = file_in  #meshlab obj input
		self.f2 = file_out #obj output librat
		self.f3 = file_out_bound #obj output librat with bounding boxes
		self.fpb = int(fpb)
		self.name = name
		self.H = 1
		self.CreateObj()
    
	def CreateObj(self):
		"""
		converts meshlab obj > librat obj
		"""	
		f1 = open(self.f1,"r")
		tmp_lines=f1.readlines()
		self.vertices=[x for x in tmp_lines if x.startswith('v')]
		self.f=[x for x in tmp_lines if x.startswith('f')]
		self.facets=[]
		self.v = np.zeros((len(self.f)*3,3))
		#create facets and extract coordinates for BB calculations
		i=0    
		for line in self.f:
			id1=int(line.split(" ")[1])-1 #indices in meshlab obj file are 1 based, not 0 based
			id2=int(line.split(" ")[2])-1
			id3=int(line.split(" ")[3])-1
			#create facet
			line="usemtl %s\n" % (self.name)
			self.facets.append(line)
			sign_digits=args.rounding
			str1='v %s %s %s\n' % (round(float(self.vertices[id1].split(" ")[1]),sign_digits),round(float(self.vertices[id1].split(" ")[2]),sign_digits),round(float(self.vertices[id1].split(" ")[3]),sign_digits))
			str2='v %s %s %s\n' % (round(float(self.vertices[id2].split(" ")[1]),sign_digits),round(float(self.vertices[id2].split(" ")[2]),sign_digits),round(float(self.vertices[id2].split(" ")[3]),sign_digits))
			str3='v %s %s %s\n' % (round(float(self.vertices[id3].split(" ")[1]),sign_digits),round(float(self.vertices[id3].split(" ")[2]),sign_digits),round(float(self.vertices[id3].split(" ")[3]),sign_digits))
			self.facets.append(str1)
			self.facets.append(str2)
			self.facets.append(str3)
			self.facets.append("f -3 -2 -1\n")
			#extract coordinates
			self.v[i][0] = round(float(self.vertices[id1].split(" ")[1]),sign_digits)
			self.v[i][1] = round(float(self.vertices[id1].split(" ")[2]),sign_digits)
			self.v[i][2] = round(float(self.vertices[id1].split(" ")[3]),sign_digits)
			self.v[i+1][0] = round(float(self.vertices[id2].split(" ")[1]),sign_digits)
			self.v[i+1][1] = round(float(self.vertices[id2].split(" ")[2]),sign_digits)
			self.v[i+1][2] = round(float(self.vertices[id2].split(" ")[3]),sign_digits)
			self.v[i+2][0] = round(float(self.vertices[id3].split(" ")[1]),sign_digits)
			self.v[i+2][1] = round(float(self.vertices[id3].split(" ")[2]),sign_digits)
			self.v[i+2][2] = round(float(self.vertices[id3].split(" ")[3]),sign_digits)
			i += 3
	
	def WriteStdObj(self):
		f2 = open(self.f2,'w')
		f2.write('!{\n')
		f2.write('g '+ self.name+' 0\n')
		f2.write('#define\n')
		i = 0
		while i < len(self.facets):
			f2.write(self.facets[i])
			i += 1
		f2.write('!}\n')
		f2.close()

		

	def DefineBBox(self):
		Max = self.v.max(axis=0)
		Min = self.v.min(axis=0)	
		xMax = Max[0] * 1.01
		xMin = Min[0] * 1.01
		yMax = Max[1] * 1.01
		yMin = Min[1] * 1.01
		zMax = Max[2] * 1.01
		zMin = Min[2] * 1.01
		i = 0	
		while i < 1:
			i = 1
			pwr = 2**self.H # number of bounding boxes in one dimension, so pwr BB in x dir, pwr BB in y dir and pwr BB in z dir
			self.Boxes = np.zeros((pwr+1,pwr+1,pwr+1,self.fpb)) # self.facets stores facets in quintets - vertices are on L1,2,3
			self.nBoxes = np.zeros((pwr+1,pwr+1,pwr+1))  # this counts the number of elements in a bounding box, if exceeds fpb it will throw an IndexError and the while loop will start again
			a = 0
			while a < len(self.facets):	
				xBox = max(0,min(pwr,int((float(self.facets[a+1].rstrip().split(" ")[1])-xMin)/((1.0/pwr)*(xMax-xMin))),int((float(self.facets[a+2].rstrip().split(" ")[1])-xMin)/((1.0/pwr)*(xMax-xMin))),int((float(self.facets[a+3].rstrip().split(" ")[1])-xMin)/((1.0/pwr)*(xMax-xMin)))))
				yBox = max(0,min(pwr,int((float(self.facets[a+1].rstrip().split(" ")[2])-yMin)/((1.0/pwr)*(yMax-yMin))),int((float(self.facets[a+2].rstrip().split(" ")[2])-yMin)/((1.0/pwr)*(yMax-yMin))),int((float(self.facets[a+3].rstrip().split(" ")[2])-yMin)/((1.0/pwr)*(yMax-yMin)))))
				zBox = max(0,min(pwr,int((float(self.facets[a+1].rstrip().split(" ")[3])-zMin)/((1.0/pwr)*(zMax-zMin))),int((float(self.facets[a+2].rstrip().split(" ")[3])-zMin)/((1.0/pwr)*(zMax-zMin))),int((float(self.facets[a+3].rstrip().split(" ")[3])-zMin)/((1.0/pwr)*(zMax-zMin)))))			
				try:
					self.Boxes[xBox][yBox][zBox][int(self.nBoxes[xBox][yBox][zBox])] = a  # put a (i.e. index of self.facets in self.Boxes, appending in its 4th dimension using the count of self.nboxes
					self.nBoxes[xBox][yBox][zBox] += 1
					a += 5
				except IndexError: # too many facets in one single BB
					if (self.H < 12):
						self.H += 1
						i -= 1
						break  # start while loop again with increased H
					else:
						sys.exit("H>7, increase -f")
		print 'H = '+str(self.H)
# 		print self.Boxes
# 		print np.sum(self.nBoxes)
		self.Final = []
		# header of BB Obj file
		self.Final.append('!{\n')
		self.Final.append('g '+ self.name+' 0\n')
		self.Final.append('#define\n')
		self.Level = np.zeros((self.H,3))

				
	def CreateBBox(self,level):
		IJK = np.zeros((3))
		for i in xrange(0,2):
			for j in xrange(0,2):
				for k in xrange(0,2):
					self.Final.append('!{')
					self.Level[level][0]=i;
					self.Level[level][1]=j;
					self.Level[level][2]=k;
					if (level==self.H-1):
						for l in xrange(0,3):
							IJK[l]=0    # reset the IJK matrix to 0
						for l in xrange(0,level+1):
							base=2**(level-l)
							for m in xrange(0,3):
								IJK[m]+=self.Level[l][m]*base
						self.Final.append('g box 1'+str(int(IJK[0]))+str(int(IJK[1]))+str(int(IJK[2])))
						for l in xrange(0,int(self.nBoxes[IJK[0]][IJK[1]][IJK[2]])):
							n = self.Boxes[IJK[0]][IJK[1]][IJK[2]][l]
							self.Final.append(self.facets[int(n)])
							self.Final.append(self.facets[int(n+1)])
							self.Final.append(self.facets[int(n+2)])
							self.Final.append(self.facets[int(n+3)])
							self.Final.append(self.facets[int(n+4)])
					else:
						self.CreateBBox(level+1)
					self.Final.append('!}')


	def WriteBBoxLibrat(self):

		self.Final.append('!}') #this is the final closure of the obj file
# 		print self.Final
		r__ = [0] * len(self.Final)
		p1__ = [0] * len(self.Final)
		r = []
		p1 = []
		__NR = -1
		for i in xrange(len(self.Final)):  # removing empties
			__NR += 1
			tmp = self.Final[i]
			tmp_s = tmp.split()
			p1__[__NR] = tmp_s[0]
			r__[__NR] = self.Final[i].strip()
			if(p1__[__NR] == '!}' and p1__[__NR-1][0] == 'g' and p1__[__NR-2] == '!{'):
				__NR -= 3  # remove 3 lines, empty BB
			elif(p1__[__NR] == '!}' and p1__[__NR-1] == '!{'):
				__NR -= 2  # remove 2 lines, empty BB
		for i in xrange(0,__NR):
			if (r__[i] != 0):
				r.append(r__[i])
				p1.append(p1__[i])
		__NR = len(r)
		NR=__NR
		nChanges = 1
		pass1 = 0
		_p1 = [0] * __NR
		_r = [0] * __NR
		while (nChanges > 0):
			pass1 += 1
			nChanges = 0
			_NR = -1
			i = 0
			while (i < NR):
				if (p1[i] == '!{' and p1[i+1] == '!}'):
					nChanges += 1
					i += 1
				elif (p1[i] =='!{' and p1[i+1] =='B' and p1[i+2] == '!}'):
					_NR += 1
					_p1[_NR] =p1[i+1]
					_r[_NR]=r[i+1]
					i+=2
					nChanges+=1
				elif (p1[i]=="!{" and p1[i+1]=="!}"):
       					nChanges+=1
					i+=1
				elif(p1[i]=="!{" and p1[i+1]=="g" and p1[i+2]=="!}"):
					nChanges+=1
					i+=2
				elif(p1[i]=="!{" and p1[i+1]=="g" and p1[i+2]!="!{" ):
					_NR += 1
					_r[_NR] = str(r[i])+"\n"+str(r[i+1])
# 					_r[_NR] = str(r[i])+str(r[i+1])
					_p1[_NR]="B"
					for j in xrange(i+2,NR+1):
						_r[_NR] = str(_r[_NR])+"\n"+str(r[j])
# 						_r[_NR] = str(_r[_NR])+str(r[j])						
						if(p1[j]=="!}"):
							i=j
							nChanges+=1
							break
				else:
					_NR += 1
					_p1[_NR] = p1[i]
					_r[_NR] = r[i]
				i += 1
			for i in xrange(0,_NR+1):
				p1[i] = _p1[i]
				r[i] = _r[i]
			NR = _NR + 1
		f3 = open(self.f3,'w')
		for i in xrange(0,NR):
#			print str(r[i])
#			print i
			f3.write(str(r[i]))
		f3.write('!}\n')
		f3.close()

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--infile',default=False,help='meshlab OBJ file: .obj format. Deselect all optional exports in meshlab export. Only v and f required.')
	parser.add_argument('-o', '--outfile',default=False,help='librat outfile in .obj format')
	parser.add_argument('-b', '--outfile_bound',default=False,help='bound librat outfile in .obj.bbox format')
	parser.add_argument('-f', '--fpb',default=100,help='Maximum number of facets per bounding box (default 100)')
	parser.add_argument('-n', '--name',default='g plant 0',help='Object name ie g plant 0')
	parser.add_argument('-r', '--rounding',default=3,help='Significant digits. Default: rounds the vertices using 3 significant digits (the nearest mm if object is in meter units)')
	args = parser.parse_args()
	if(args.infile == False):
		sys.exit('No input file specified. -h for usage')
	elif (args.outfile != False and args.outfile_bound != False):
		example = mesh2librat(args.infile,args.outfile,args.outfile_bound,args.fpb,args.name)
		example.WriteStdObj()
		example.DefineBBox()
		example.CreateBBox(level=0)
		example.WriteBBoxLibrat()
	elif (args.outfile != False and args.outfile_bound == False):
		example = mesh2librat(args.infile,args.outfile,args.outfile_bound,args.fpb,args.name)
		example.WriteStdObj()
	elif (args.outfile == False and args.outfile_bound != False):
		example = mesh2librat(args.infile,args.outfile,args.outfile_bound,args.fpb,args.name)
		example.DefineBBox()
		example.CreateBBox(level=0)
		example.WriteBBoxLibrat()
