# test zonal stats on points

import os, sys
import numpy as np
sys.path.append('./')
import zonepy as zp

gdb = os.path.join('..', '..', '..', 'scratch', 'testProjection.shp')
ras = os.path.join('examples','zonepy_test.tif')
lyrName = None 
filenm = os.path.join('examples', 'outputfile')

zs = zp.zonal.ZoneClass(gdb, ras, lyrName=lyrName, fldname=None , 
                 projIn=None, projOut=None, buffDist=100, fact=30, 
                 outND=np.nan, nd_thresh=100, filenm=filenm, csvout=True)

zs.compute_stats()