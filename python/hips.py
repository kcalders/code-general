import os
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import scipy.misc



def read_header(fname):

    """
    :param fname:
    :return:

    header_length, bands, res_x, res_y, fmt
    """

    # grab header info
    header_length = open(fname).read().find('\n.')
    header = open(fname).read()[:header_length].split()
    bands = int(header[1])
    res_x = int(header[2])
    res_y = int(header[3])
    fmt = int(header[4])

    return header_length, bands, res_x, res_y, fmt

def read_hips(fname):

    header_length, bands, res_x, res_y, fmt = read_header(fname)

    # extract image from hips
    hips = np.fromfile(fname, np.float32)
    hips_length = len(hips)
    img = hips[hips_length - (bands * res_x * res_y):].reshape(bands, res_x, res_y)
    img = np.rot90(img.T, 3)

    return img, bands, res_x, res_y, fmt


def hipstats(fname, unique=False):

    img, bands, res_x, res_y, fmt = read_hips(fname)
    out = [img.min(), img.max(), img.mean(), img.std()]
    if unique:
        out.append(np.unique(img))

    return out


def hips2img_g(fname, order=[0], stretch=True, imshow=True,
             imsave=False, ax=None):
    
    img, bands, res_x, res_y, fmt = read_hips(fname)
    
    #stretch
    for b in range(bands):
        arr_b = img[:, :, b]
        if stretch:
            arr_b = ((arr_b - np.percentile(arr_b, 2.5)) / np.percentile(arr_b, 97.5))
            arr_b[arr_b < 0] = 0
            arr_b[arr_b > 1] = 1
        img[:, :, b] = arr_b

    # display image
    if not ax:
        fig, ax = plt.subplots(figsize=(10, 10))

    if len(order) == 1 or bands == 1:
        order = order[0]
        cmap = 'gray'
    else:
        cmap = 'spectral'

    ax.imshow(img[:, :, order], cmap=cmap, interpolation='none')
    ax.axis('off')
    
    # save image
    if imsave:
        plt.imsave(os.path.splitext(fname)[0] + '.png', img[:, :, order],
                   cmap=cmap)

    # plot image to screen
    if imshow:
        plt.show()

    return ax
    
def hips2img_g(fname, order=[0], imshow=False, imsave=True, ax=None):
    
    img, bands, res_x, res_y, fmt = read_hips(fname)
    
    for b in range(bands):
    	arr_b = img[:, :, b]
    	#    arr_b[arr_b == 0] = 1  #0 is ~ NA in these reflectance images
    	print "band %i, min value: %f, max value: %f, mean value: %f" %(b, arr_b.min(), arr_b.max(), arr_b.mean())
#     	print "converted to:"
#     	img[:, :, b] = (arr_b*255).astype(np.uint8)
#     	print "band %i, min value: %f, max value: %f, mean value: %f \n" %(b, img[:, :, b].min(), img[:, :, b].max(), img[:, :, b].mean())

    # display image
    if not ax:
        fig, ax = plt.subplots(figsize=(10, 10))


    cmap = 'gray' #colormap

    ax.imshow(img[:, :, order[0]], cmap=cmap, interpolation='none')
    ax.axis('off')
    
    # save image
    if imsave:
        plt.imsave(os.path.splitext(fname)[0] + '.jpeg', img[:, :, order[0]],
                   cmap=cmap)

    # plot image to screen
    if imshow:
        plt.show()

    return ax


def hips2img_rgb(fname, order=[2,1,0], imshow=False, imsave=True):
    
    img, bands, res_x, res_y, fmt = read_hips(fname)

    #convert 0-1 range to 0-255
    for b in range(bands):
    	arr_b = img[:, :, b]
    	arr_b[arr_b == 0] = 1  #0 is ~ NA in these reflectance images
    	print "band %i, min value: %f, max value: %f, mean value: %f" %(b, arr_b.min(), arr_b.max(), arr_b.mean())
    	print "converted to:"
    	img[:, :, b] = (arr_b*255).astype(np.uint8)
    	print "band %i, min value: %f, max value: %f, mean value: %f \n" %(b, img[:, :, b].min(), img[:, :, b].max(), img[:, :, b].mean())

    # http://stackoverflow.com/questions/10443295/combine-3-separate-numpy-arrays-to-an-rgb-image-in-pytho
	# http://stackoverflow.com/questions/9193603/applying-a-coloured-overlay-to-an-image-in-either-pil-or-imagemagik
	rgbArray = np.zeros((res_x,res_y,3), 'uint8') #uint8 has a dynamic range of [0 255]
    rgbArray[..., 0] = img[:, :, order[0]] #r
    rgbArray[..., 1] = img[:, :, order[1]] #g
    rgbArray[..., 2] = img[:, :, order[2]] #b
    img_rgb=Image.fromarray(rgbArray, 'RGB')
    
    if imsave:  #$save image
    	img_rgb.save(os.path.splitext(fname)[0] + '.jpeg')
	
     # plot image to screen
    if imshow:
         img_rgb.show()

