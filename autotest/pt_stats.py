# test zonal stats on points

import os, sys
import numpy as np
sys.path.append('./')
import zonepy as zp

gdb = os.path.join('examples', 'data', 'zonepy_pts.gpkg')
ras = os.path.join('examples', 'data', 'zone10kmCellNum.tif')
lyrName = 'zonepy_pts'
filenm = os.path.join('examples', 'data', 'outputfilepoly')

zs = zp.ZoneClass(gdb, ras, lyrName=lyrName, fldname=None , 
                 buffDist=100, fact=30, 
                 outND=np.nan, nd_thresh=100, filenm=filenm)

zs.compute_stats()
zs.writeSHP()
zs.writeCSV()
zs.writePKL()