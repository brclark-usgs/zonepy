import os,sys
import numpy as np
import pandas as pd
from osgeo import gdal, ogr, osr
from osgeo.gdalconst import *


'''
Export pickled pandas object to raster, appropriate for continuous 
vector datasets
'''

def RasCreate(pklfile, vec, lyrName=None, stat='mean', 
    nodata=-9999,  outTiff='outras.tif'):
    '''
    Define numpy array of raster shape and write raster

    Parameters
    ----------
    pklfile: str, Filepath of pickle file
    vec: str, Filepath of vector zones used to create pickle file (must be 
            continuous vector)
    lyrName: str, default None, Name of feature class of points. Default 
            will select first layer in geodatabase. Also use default for 
            shapefiles.
    stat: str, default = 'mean', Summary statistic to use for raster creation.
            Options are 'min' (minimum), 'mean', 'max' (maximum), 'std' 
            (standard deviation), 'sum', 'count', 'median'
    nodata: float, default=-9999., NoData value to use in raster
    outTiff: str, Filepath and name of created raster
    '''

    # Open feature class
    vds = ogr.Open(vec)  
    if lyrName != None:
        vlyr = vds.GetLayerByName(lyrName)
    else:
        vlyr = vds.GetLayer(0)
    # print(vlyr)

    # Get extent of single polygon
    feat = vlyr.GetNextFeature()
    featext = feat.GetGeometryRef().GetEnvelope()
    print(featext)
    xsz = featext[1] - featext[0] #resolution in x
    ysz = featext[3] - featext[2] #resolution in y
    sz = (xsz, ysz)
    print(sz)

    # Get projection and extent
    vsr = vlyr.GetSpatialRef()
    # print(vsr)
    vext = vlyr.GetExtent()
    cols = int((vext[1] - vext[0]) / xsz)
    rows = int((vext[3] - vext[2]) / ysz)
    print(rows,cols)

    # Read pickle file of summary stats
    df = pd.read_pickle(pklfile)
    # print(df.shape)
    a = df[stat].values

    # Reshape to zone shape
    a = a.reshape((rows,cols))

    # Write raster
    writeRaster(a, vext, sz, vsr, outTiff=outTiff, nodata=nodata)


def writeRaster(grid, ext, sz, prj, outTiff='outras.tif', nodata=-9999):
    '''
    Create raster from numpy array and extent
    
    Parameters
    ----------
    grid: numpy array, numpy array of raster to be created with correct shape 
            (rows and columns)
    ext: tuple, tuple of floats for raster extent, ordered as
            x minimum, x maximum, y minimum, y maximum
    sz: tuple, tuple of floats for raster call resoluation, ordered
            as x size, y size
    prj: well-known text object, WKT projection info from vector (zone)
    outTiff: str, default 'outras.tif', Filepath of raster to be created.
        Default creates file in working directory.
    nodata: float, default -9999., no data value to write out to raster
    '''

    print('writing tiff...')

    # Calculate angle and coordinates
    angle = 0
    angle = angle * -1. * np.pi / 180.
    geotransform = (vext[0], sz[0] * np.cos(angle), -sz[0] * np.sin(angle),
                    vext[3], -sz[1] * np.sin(angle), -sz[1] * np.cos(angle)) 

    # Create raster to hold data
    nrow = grid.shape[0]
    ncol = grid.shape[1]
    drv = gdal.GetDriverByName('GTiff')
    ds = drv.Create(outTiff, ncol, nrow, 1 ,gdal.GDT_Float32)

    # Set raster no data value
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)

    # Set coordinates and projection
    ds.SetGeoTransform(geotransform)  
    ds.SetProjection(str(vsr))

    # Write array to raster
    band.WriteArray(grid)


if __name__ == '__main__':

    dws = r'D:\Projects\2015_NAQWA_METX' #main workspace
    znpath = os.path.join(dws, 'METX_GIS', 'ZonalTests') #folder location
    flname = os.path.join(znpath, 'zonetest_farmN.pkl') #pickle file
    vec =  os.path.join(znpath, 'natgrid_gpck', 'natGrid1km2.gpkg') #vector file (continuous)
    outTiff = os.path.join(znpath, 'Nfarm_ras.tif')

    RasCreate(flname, vec, lyrName='modelgrid', stat='mean', 
        nodata=-9999,  outTiff=outTiff)


   
