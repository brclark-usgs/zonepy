# test zonal stats on points

import os, sys
import numpy as np
sys.path.append('./')
import zonepy as zp

gdb = os.path.join('examples', 'data', 'zonepy_pts.gpkg')
ras = os.path.join('examples', 'data', 'zone10kmCellNum.tif')
lyrName = 'zonepy_pts'
filenm = os.path.join('examples', 'data', 'outputfile')

zc = zp.ZoneClass(gdb, ras, lyrName=lyrName, fldname=None , 
                  buffDist=100, cmap=None, fact=30, 
                  filenm=filenm)

zc.compute_category()

zc.writeSHP()
zc.writeCSV()
zc.writePKL()
