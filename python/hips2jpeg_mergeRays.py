#!/usr/bin/env python
#THIS WILL PUT DIFFERENT HIPS IMAGES (e.g. 1 rpp) TOGETHER TO CREATE A HIGHER RPP. IMAGE

from hips import *
from hips import read_header
import optparse
import sys
import glob
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np


def convertDN(arr,stretch):
    #convert 0-1 range to 0-255
#     print "min value: %f, max value: %f, mean value: %f" %(arr.min(), arr.max(), arr.mean())
#     plt.hist(arr)
#     plt.show()
    arr = (arr*1/stretch)#.astype(np.uint8)
    arr[arr == 0] = 1  #0 is ~ NA in these reflectance images
    arr[arr > 1] = 1  #0 is ~ NA in these reflectance images
#     print "min value: %f, max value: %f, mean value: %f \n" %(arr.min(), arr.max(), arr.mean())
#     plt.hist(arr)
#     plt.show()
    return arr

def rpp_merge(rpp, order=[2,1,0], imsave=True):
    
    rpp=int(rpp)
    file_list = glob.glob("*hips")
    firstFile=True
  
    for file in file_list[0:rpp]:
    	print file
    	img, bands, res_x, res_y, fmt = read_hips(file)

    	if firstFile:
        	res_xx=res_x
        	res_yy=res_y
        	arr_r = img[:, :, order[0]]
        	arr_g = img[:, :, order[1]]
        	arr_b = img[:, :, order[2]]
        	firstFile=False
    	elif res_x==res_xx and res_y==res_yy:
        	arr_r = arr_r + img[:, :, order[0]]
        	arr_g = arr_g + img[:, :, order[1]]
        	arr_b = arr_b + img[:, :, order[2]]
    	else:
        	print "not all image dimensions matching"
        	sys.exit()

	#you don't have to average as it is eventually stretched to 0-1 anyway        
	# arr_r=arr_r/rpp
	# arr_g=arr_g/rpp
	# arr_b=arr_b/rpp

    print "out"
    stretch_value=np.max((np.percentile(arr_r, 99),np.percentile(arr_g, 99),np.percentile(arr_b, 99)))
	#rescale all 3 bands the same using the highest 99th percentile value:
	# print "Conversion of RED band"
    arr_r=convertDN(arr_r,stretch_value)
	# print "Conversion of GREEN band"
    arr_g=convertDN(arr_g,stretch_value)
	# print "Conversion of BLUE band"
    arr_b=convertDN(arr_b,stretch_value)


    rgbArray = np.zeros((res_x,res_y,3), 'uint8') #uint8 has a dynamic range of [0 255]
    rgbArray[..., 0] = arr_r*255 #r
    rgbArray[..., 1] = arr_g*255 #g
    rgbArray[..., 2] = arr_b*255 #b
#     img_rgb=np.dstack((arr_r,arr_g,arr_b))
    print "helloe"
    img_rgb=Image.fromarray(rgbArray, 'RGB')

# 	fig, ax = plt.subplots(figsize=(10, 10))
# 	ax.imshow(img_rgb, cmap='spectral')
# 	plt.show()
    if imsave:  
	    img_rgb.save(file_list[0].split("_rpp")[0] +'_rpp_' + str(rpp) + '.jpeg')
# 	    plt.imsave(file_list[0].split("_rpp")[0] +'_rpp_' + str(rpp) + '.2.jpeg',img_rgb2,cmap='spectral')


# MAIN PROGRAM
def main(cmdargs):
    rpp_merge(cmdargs.rpp, order=[6,3,1], imsave=cmdargs.save)


# Command arguments
class CmdArgs:
  def __init__(self):
    p = optparse.OptionParser()
#     p.add_option("-i","--inFile", dest="inFile", default=None, help="Input hips file")
    p.add_option("-c","--save", dest="save", default=True, help="Save as jpeg. Default True")
    p.add_option("-r","--rpp", dest="rpp", default=1, help="rays per pixel to include")
    (options, args) = p.parse_args()
    self.__dict__.update(options.__dict__)
    
#     if self.inFile is None:
#         p.print_help()
#         print "Input hips must be specified."
#         sys.exit()


# Run the script
if __name__ == "__main__":
    cmdargs = CmdArgs()
    main(cmdargs)
