#!/usr/bin/env python

#Andrew Burt - a.burt.12@ucl.ac.uk

import argparse
import numpy
import math
import scipy.io
import scipy.linalg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import mpl_toolkits.mplot3d
import sys

def cylinder(p0,p1,R):
	v = p1 - p0
	mag = scipy.linalg.norm(v)
	v = v / mag
	not_v = numpy.array([1, 0, 0])
	if (v == not_v).all():
	    not_v = numpy.array([0, 1, 0])
	n1 = numpy.cross(v, not_v)
	n1 /= scipy.linalg.norm(n1)
	n2 = numpy.cross(v, n1)
	t = numpy.linspace(0, mag, 100)
	theta = numpy.linspace(0, 2 * numpy.pi, 100)
	t,theta = numpy.meshgrid(t, theta)
	X,Y,Z = [p0[i] + v[i] * t + R * numpy.sin(theta) * n1[i] + R * numpy.cos(theta) * n2[i] for i in [0, 1, 2]]
	return X,Y,Z

def roundd(num,divisor,up=True):
	if(up == True):
		return num + (num%divisor)
	else:
		return num - (num%divisor)

def plotClouds(clouds,azimuth,elevation,pdf,out):

# 	fig_width_pt = 418.82394
# 	inches_per_pt = 1.0/72.27
# 	golden_mean = (numpy.sqrt(5.0)-1.0)/2.0
# 	fig_width = fig_width_pt*inches_per_pt*1
# 	fig_height = fig_width*golden_mean
# 	fig_size = [fig_width,fig_height]
# 	params = {
# 			"backend": "pdflatex",
# 			"text.usetex": True,
# 			"font.family": "serif",
# 			"font.serif": [],
# 			"font.sans-serif": [],
# 			"font.monospace": [],
# 			"axes.labelsize": 10,
# 			"text.fontsize": 10,
# 			"legend.fontsize": 8,
# 			"xtick.labelsize": 8,
# 			"ytick.labelsize": 8,
# 			"figure.figsize": fig_size
# 		}
# 	matplotlib.rcParams.update(params)
	fig = matplotlib.pyplot.figure(frameon=False)
	ax = fig.add_subplot(111, projection='3d')
	ax.view_init(azim=azimuth,elev=elevation)
	matplotlib.pyplot.axis('off')
	x_min = float("inf")
	x_max = 0
	y_min = float("inf")
	y_max = 0
	z_min = float("inf")
	z_max = 0
	for i in xrange(len(clouds)):
		data = numpy.loadtxt(clouds[i])
		if(len(clouds) > 1):
			c = numpy.random.rand(3,1)
		else:
			c = 'k'
		ax.scatter(data[:,0],data[:,1],data[:,2],alpha=1,s=0.01,lw=0,color=c)
		for j in xrange(len(data)):
			if(data[j][0] < x_min):
				x_min = data[j][0]
			if(data[j][0] > x_max):
				x_max = data[j][0]
			if(data[j][1] < y_min):
				y_min = data[j][1]
			if(data[j][1] > y_max):
				y_max = data[j][1]
			if(data[j][2] < z_min):
				z_min = data[j][2]
			if(data[j][2] > z_max):
				z_max = data[j][2]
	max_range = numpy.array([x_max-x_min,y_max-y_min,z_max-z_min]).max()
	mean_x = (x_max+x_min)/2
	mean_y = (y_max+y_min)/2
	mean_z = (z_max+z_min)/2
	ax.set_xlim(roundd(mean_x-max_range/2,5,False),roundd(mean_x+max_range/2,5,True))
	ax.set_ylim(roundd(mean_y-max_range/2,5,False),roundd(mean_y+max_range/2,5,True))
	ax.set_zlim(roundd(mean_z-max_range/2,5,False),roundd(mean_z+max_range/2,5,True))
	if(pdf == False):
		matplotlib.pyplot.savefig(out+'_clouds.png',dpi=2400,bbox_inches='tight')
	else:
		matplotlib.pyplot.savefig(out+'_clouds.pdf')
		
		
def plot2Clouds(clouds,wood,leaf,azimuth,elevation,pdf,out):

	fig = matplotlib.pyplot.figure(frameon=False)
	ax = fig.add_subplot(111, projection='3d')
	ax.view_init(azim=azimuth,elev=elevation)
	matplotlib.pyplot.axis('off')
	x_min = float("inf")
	x_max = 0
	y_min = float("inf")
	y_max = 0
	z_min = float("inf")
	z_max = 0
	for i in xrange(len(clouds)):
		data = numpy.loadtxt(clouds[i])
		dwood = numpy.loadtxt(wood[i])
		dleaf = numpy.loadtxt(leaf[i])
		c = 'k'
		ax.scatter(dwood[:,0],dwood[:,1],dwood[:,2],alpha=0.6,s=0.01,lw=0,color=c)
		c = 'g'
		ax.scatter(dleaf[:,0],dleaf[:,1],dleaf[:,2],alpha=0.6,s=0.01,lw=0,color=c)
		for j in xrange(len(data)):
			if(data[j][0] < x_min):
				x_min = data[j][0]
			if(data[j][0] > x_max):
				x_max = data[j][0]
			if(data[j][1] < y_min):
				y_min = data[j][1]
			if(data[j][1] > y_max):
				y_max = data[j][1]
			if(data[j][2] < z_min):
				z_min = data[j][2]
			if(data[j][2] > z_max):
				z_max = data[j][2]
	max_range = numpy.array([x_max-x_min,y_max-y_min,z_max-z_min]).max()
	mean_x = (x_max+x_min)/2
	mean_y = (y_max+y_min)/2
	mean_z = (z_max+z_min)/2
	ax.set_xlim(roundd(mean_x-max_range/2,5,False),roundd(mean_x+max_range/2,5,True))
	ax.set_ylim(roundd(mean_y-max_range/2,5,False),roundd(mean_y+max_range/2,5,True))
	ax.set_zlim(roundd(mean_z-max_range/2,5,False),roundd(mean_z+max_range/2,5,True))
	if(pdf == False):
		matplotlib.pyplot.savefig(out+'_separated_clouds.png',dpi=2400,bbox_inches='tight')
	else:
		matplotlib.pyplot.savefig(out+'_separated_clouds.pdf')
		fig = matplotlib.pyplot.figure(frameon=False)
		
	ax = fig.add_subplot(111, projection='3d')
	ax.view_init(azim=azimuth,elev=elevation)
	matplotlib.pyplot.axis('off')
	x_min = float("inf")
	x_max = 0
	y_min = float("inf")
	y_max = 0
	z_min = float("inf")
	z_max = 0
	for i in xrange(len(clouds)):
		data = numpy.loadtxt(clouds[i])
		dwood = numpy.loadtxt(wood[i])
		dleaf = numpy.loadtxt(leaf[i])
		c = 'k'
		ax.scatter(dwood[:,0],dwood[:,1],dwood[:,2],alpha=0.6,s=0.01,lw=0,color=c)
		for j in xrange(len(data)):
			if(data[j][0] < x_min):
				x_min = data[j][0]
			if(data[j][0] > x_max):
				x_max = data[j][0]
			if(data[j][1] < y_min):
				y_min = data[j][1]
			if(data[j][1] > y_max):
				y_max = data[j][1]
			if(data[j][2] < z_min):
				z_min = data[j][2]
			if(data[j][2] > z_max):
				z_max = data[j][2]
	max_range = numpy.array([x_max-x_min,y_max-y_min,z_max-z_min]).max()
	mean_x = (x_max+x_min)/2
	mean_y = (y_max+y_min)/2
	mean_z = (z_max+z_min)/2
	ax.set_xlim(roundd(mean_x-max_range/2,5,False),roundd(mean_x+max_range/2,5,True))
	ax.set_ylim(roundd(mean_y-max_range/2,5,False),roundd(mean_y+max_range/2,5,True))
	ax.set_zlim(roundd(mean_z-max_range/2,5,False),roundd(mean_z+max_range/2,5,True))
	if(pdf == False):
		matplotlib.pyplot.savefig(out+'_wood_clouds.png',dpi=2400,bbox_inches='tight')
	else:
		matplotlib.pyplot.savefig(out+'_wood_clouds.pdf')

		


def plotModels(models,azimuth,elevation,pdf,out):

# 	fig_width_pt = 418.82394
# 	inches_per_pt = 1.0/72.27
# 	golden_mean = (numpy.sqrt(5.0)-1.0)/2.0
# 	fig_width = fig_width_pt*inches_per_pt*1
# 	fig_height = fig_width*golden_mean
# 	fig_size = [fig_width,fig_height]
# 	params = {
# 			"backend": "pdflatex",
# 			"text.usetex": True,
# 			"font.family": "serif",
# 			"font.serif": [],
# 			"font.sans-serif": [],
# 			"font.monospace": [],
# 			"axes.labelsize": 10,
# 			"text.fontsize": 10,
# 			"legend.fontsize": 8,
# 			"xtick.labelsize": 8,
# 			"ytick.labelsize": 8,
# 			"figure.figsize": fig_size
# 		}
# 	matplotlib.rcParams.update(params)
	colors = "bgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykw"
	fig = matplotlib.pyplot.figure(frameon=False)
	ax = fig.add_subplot(111, projection='3d')
	ax.view_init(azim=azimuth,elev=elevation)
	matplotlib.pyplot.axis('off')
	x_min = float("inf")
	x_max = 0
	y_min = float("inf")
	y_max = 0
	z_min = float("inf")
	z_max = 0
	for i in xrange(len(models)):
		data = scipy.io.loadmat(models[i])
		if(len(models) > 1):
			c = numpy.random.rand(3,1)
		else:
			c = 'k'
		for j in xrange(len(data['Rad'])):
			R = data['Rad'][j][0]
			p0 = numpy.array([data['Sta'][j][0],data['Sta'][j][1],data['Sta'][j][2]])
			p1 = numpy.array([data['Sta'][j][0]+data['Len'][j][0]*data['Axe'][j][0],data['Sta'][j][1]+data['Len'][j][0]*data['Axe'][j][1],data['Sta'][j][2]+data['Len'][j][0]*data['Axe'][j][2]])
			if(p0[0] != p1[0]):
				X,Y,Z = cylinder(p0,p1,R)
				cols=data['BoC'][j][1]
				ax.plot_surface(X,Y,Z,alpha=0.33,linewidth=0,color=colors[cols])
				for k in xrange(len(X)):
					for l in xrange(len(X[k])):
						if(X[k][l] < x_min):
							x_min = X[k][l]
						if(X[k][l] > x_max):
							x_max = X[k][l]
						if(Y[k][l] < y_min):
							y_min = Y[k][l]
						if(Y[k][l] > y_max):
							y_max = Y[k][l]
						if(Z[k][l] < z_min):
							z_min = Z[k][l]
						if(Z[k][l] > z_max):
							z_max = Z[k][l]
	max_range = numpy.array([x_max-x_min,y_max-y_min,z_max-z_min]).max()
	mean_x = (x_max+x_min)/2
	mean_y = (y_max+y_min)/2
	mean_z = (z_max+z_min)/2
	ax.set_xlim(roundd(mean_x-max_range/2,5,False),roundd(mean_x+max_range/2,5,True))
	ax.set_ylim(roundd(mean_y-max_range/2,5,False),roundd(mean_y+max_range/2,5,True))
	ax.set_zlim(roundd(mean_z-max_range/2,5,False),roundd(mean_z+max_range/2,5,True))
	if(pdf == False):
		matplotlib.pyplot.savefig(out+'_models.png',dpi=2400,bbox_inches='tight')
	else:
		matplotlib.pyplot.savefig(out+'_models.pdf')

def plotCloudsModels(clouds,models,azimuth,elevation,pdf,out):

# 	fig_width_pt = 418.82394
# 	inches_per_pt = 1.0/72.27
# 	golden_mean = (numpy.sqrt(5.0)-1.0)/2.0
# 	fig_width = fig_width_pt*inches_per_pt*1
# 	fig_height = fig_width*golden_mean
# 	fig_size = [fig_width,fig_height]
# 	params = {
# 			"backend": "pdflatex",
# 			"text.usetex": True,
# 			"font.family": "serif",
# 			"font.serif": [],
# 			"font.sans-serif": [],
# 			"font.monospace": [],
# 			"axes.labelsize": 10,
# 			"text.fontsize": 10,
# 			"legend.fontsize": 8,
# 			"xtick.labelsize": 8,
# 			"ytick.labelsize": 8,
# 			"figure.figsize": fig_size
# 		}
# 	matplotlib.rcParams.update(params)
	fig = matplotlib.pyplot.figure(frameon=False)
	ax = fig.add_subplot(111, projection='3d')
	ax.view_init(azim=azimuth,elev=elevation)
	matplotlib.pyplot.axis('off')
	x_min = float("inf")
	x_max = 0
	y_min = float("inf")
	y_max = 0
	z_min = float("inf")
	z_max = 0
	for i in xrange(len(clouds)):
		data = numpy.loadtxt(clouds[i])
		if(len(clouds) > 1):
			c = numpy.random.rand(3,1)
		else:
			c = 'k'
		ax.scatter(data[:,0],data[:,1],data[:,2],alpha=0.33,s=0.01,lw=0,color=c)
		for j in xrange(len(data)):
			if(data[j][0] < x_min):
				x_min = data[j][0]
			if(data[j][0] > x_max):
				x_max = data[j][0]
			if(data[j][1] < y_min):
				y_min = data[j][1]
			if(data[j][1] > y_max):
				y_max = data[j][1]
			if(data[j][2] < z_min):
				z_min = data[j][2]
			if(data[j][2] > z_max):
				z_max = data[j][2]
	colors = "bgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykw"
	for i in xrange(len(models)):
		data = scipy.io.loadmat(models[i])
		if(len(models) > 1):
			c = numpy.random.rand(3,1)
		else:
			c = 'r'
		for j in xrange(len(data['Rad'])):
			R = data['Rad'][j][0]
			p0 = numpy.array([data['Sta'][j][0],data['Sta'][j][1],data['Sta'][j][2]])
			p1 = numpy.array([data['Sta'][j][0]+data['Len'][j][0]*data['Axe'][j][0],data['Sta'][j][1]+data['Len'][j][0]*data['Axe'][j][1],data['Sta'][j][2]+data['Len'][j][0]*data['Axe'][j][2]])
			if(p0[0] != p1[0]):
				X,Y,Z = cylinder(p0,p1,R)
				cols=data['BoC'][j][1]
				ax.plot_surface(X,Y,Z,alpha=0.66,linewidth=0,color=colors[cols])
	max_range = numpy.array([x_max-x_min,y_max-y_min,z_max-z_min]).max()
	mean_x = (x_max+x_min)/2
	mean_y = (y_max+y_min)/2
	mean_z = (z_max+z_min)/2
	ax.set_xlim(roundd(mean_x-max_range/2,5,False),roundd(mean_x+max_range/2,5,True))
	ax.set_ylim(roundd(mean_y-max_range/2,5,False),roundd(mean_y+max_range/2,5,True))
	ax.set_zlim(roundd(mean_z-max_range/2,5,False),roundd(mean_z+max_range/2,5,True))
	if(pdf == False):
		matplotlib.pyplot.savefig(out+'_cloudsmodels.png',dpi=2400,bbox_inches='tight')
	else:
		matplotlib.pyplot.savefig(out+'cloudsmodels.pdf')

if __name__ == "__main__":

    	parser = argparse.ArgumentParser()
    	parser.add_argument('-c','--clouds',nargs='*',default=False,help='ASCII xyx cloud')
    	parser.add_argument('-w','--wood',nargs='*',default=False,help='ASCII wood xyx cloud')
    	parser.add_argument('-l','--leaf',nargs='*',default=False,help='ASCII leaf xyx cloud')
    	parser.add_argument('-m','--models',nargs='*',default=False,help='.mat structural models')
	parser.add_argument('-ol','--overlay',default=False,action='store_true',help='cloud overlain model')
	parser.add_argument('-az','--azimuth',type=int,default=0,help='azimuth angle')
	parser.add_argument('-el','--elevation',type=int,default=0,help='elevation angle')
	parser.add_argument('-out','--outname',type=str,default="test",help='output name')
	parser.add_argument('-pdf','--pdf',default=False,action='store_true',help='print as vector')
	args = parser.parse_args()
	if(args.clouds != False and args.overlay == False):
		plotClouds(args.clouds,args.azimuth,args.elevation,args.pdf,args.outname)
	if(args.models != False and args.overlay == False):
		plotModels(args.models,args.azimuth,args.elevation,args.pdf,args.outname)
	if(args.overlay == True and args.clouds != False and args.models != False):
		plotCloudsModels(args.clouds,args.models,args.azimuth,args.elevation,args.pdf,args.outname)
	if(args.wood != False and args.leaf != False and args.clouds != False):
		plot2Clouds(args.clouds,args.wood,args.leaf,args.azimuth,args.elevation,args.pdf,args.outname)
