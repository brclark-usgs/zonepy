# test extract razter value to points

import os, sys
import numpy as np
sys.path.append('./')
import zonepy as zp

gdb = os.path.join('examples', 'zonepy_pts.gpkg')
ras = os.path.join('examples','zone10kmCellNum.tif')
lyrName = 'zonepy_pts'
filenm = os.path.join('examples', 'outputfile')

zs = zp.ZoneClass(gdb, ras, lyrName=lyrName, fldname=None , 
                 buffDist=100, fact=30, 
                 outND=np.nan, nd_thresh=100, filenm=filenm, output='csv')

zs.extractByPoint()