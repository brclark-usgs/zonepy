# test zonal stats on points

import os, sys
import numpy as np
sys.path.append('./')
import zonepy as zp

gdb = os.path.join('examples', 'data', 'zone10km.gpkg')
ras = os.path.join('examples', 'data', 'zone10kmCellNum.tif')
lyrName = 'modelgrid'
filenm = os.path.join('examples', 'data', 'outputfilepoly')

zs = zp.ZoneClass(gdb, ras, lyrName=lyrName, fldname=None , 
                 buffDist=100, fact=30, 
                 outND=np.nan, nd_thresh=100, filenm=filenm, output='csv')

zs.compute_stats()