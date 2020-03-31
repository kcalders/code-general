#!/usr/bin/env python

from hips import *
from hips import read_header
import optparse
import sys

# MAIN PROGRAM
def main(cmdargs):
   print read_header(cmdargs.inFile)
   if cmdargs.mode == "gray":
   	hips2img_g(cmdargs.inFile, order=[0], imsave=cmdargs.save,imshow=cmdargs.show)
   elif cmdargs.mode == "rgb":
	hips2img_rgb(cmdargs.inFile, order=[2,1,0], imsave=cmdargs.save,imshow=cmdargs.show)
   else:
	print "mode must be specified as gray or rgb"

# Command arguments
class CmdArgs:
  def __init__(self):
    p = optparse.OptionParser()
    p.add_option("-i","--inFile", dest="inFile", default=None, help="Input hips file")
    p.add_option("-m","--mode", dest="mode", default="gray", help="gray or rgb")
    p.add_option("-v","--show", dest="show", default=False, help="Stretch. Default False")
    p.add_option("-c","--save", dest="save", default=True, help="Save as jpeb. Default True")
    (options, args) = p.parse_args()
    self.__dict__.update(options.__dict__)
    
    if self.inFile is None:
        p.print_help()
        print "Input hips must be specified."
        sys.exit()


# Run the script
if __name__ == "__main__":
    cmdargs = CmdArgs()
    main(cmdargs)
