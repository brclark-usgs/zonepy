# test zonal stats on points

import os, sys
import numpy as np
sys.path.append('./')
import zonepy as zp

gdb = os.path.join('examples', 'zone10km.gpkg')
ras = os.path.join('examples','zone10kmCellNum.tif')
lyrName = 'modelgrid' 
inputfile = os.path.join('examples', 'outputfilepoly.csv')
outTiff = os.path.join('examples', 'statRaster.tif')

zs = zp.ZoneClass(gdb, ras, lyrName=lyrName, fldname=None , 
                 buffDist=0, fact=1, 
                 outND=np.nan, nd_thresh=100, filenm=None, output='csv')

zs.RasCreate(outTiff=outTiff, inputfile=inputfile)