# test zonal stats on points

import os, sys
import numpy as np
sys.path.append('./')
import zonepy as zp

gdb = os.path.join('examples', 'zonepy_pts.gpkg')
ras = os.path.join('examples','zonepy_test.tif')
lyrName = 'zonepy_pts'
filenm = os.path.join('examples', 'outputfile')

zc = zp.zonal.ZoneCat(gdb, ras, lyrName=lyrName, fldname=None , 
	projIn=None, projOut=None, buffDist=100, cmap=None, fact=30, 
	filenm=filenm, csvout=True)
zc.compute_category()

