# test extract razter value to points

import os, sys
import numpy as np
sys.path.append('./')
import zonepy as zp

gdb = os.path.join('examples', 'data', 'zonepy_pts.gpkg')
ras = os.path.join('examples', 'data', 'zone10kmCellNum.tif')
lyrName = 'zonepy_pts'
filenm = os.path.join('examples', 'data', 'outputfile')

zs = zp.ZoneClass(gdb, ras, lyrName=lyrName, fldname=None , 
                 buffDist=0, fact=30, 
                 outND=np.nan, nd_thresh=100, filenm=filenm, output='csv')

zs.extractByPoint()
zs.writeSHP()
zs.writeCSV()
zs.writePKL()