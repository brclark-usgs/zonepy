from __future__ import division
from __future__ import print_function
'''
Zonal Analysis

Developed by Brian Clark, Katherine Knierim, and Leslie Duncan. Portions of 
this code were modified from Copyright 2013 Matthew Perry, which were 
licensed under BSD-3 and included in this repo.


'''
import os, sys, cmath
from osgeo import gdal, ogr, osr
from osgeo.gdalconst import *
import numpy as np
import pandas as pd
from math import radians, degrees, atan
import warnings
gdal.PushErrorHandler('CPLQuietErrorHandler')
 
class ZoneStat(object):
    """
    Parameters
    ----------
    gdb: str, Filepath of vector to use for zonal analysis; 
            shapefile or geodatabase containing feature class, either points 
            or polygons 
    ras: str, Filepath of raster to use for zonal analysis;
        geotif format
    lyrName: str, default None, Name of feature class of points. Default 
            will select first layer in geodatabase. Also use default for 
            shapefiles.
    fldname: str, default None, Unique identifier field of vector data. 
            Default will select first column in feature class
    projIn: proj4 string, default None, Defines projection of vector data
    projOut: proj4 string, default None, Used to transform the vector data
             to match the projection of the raster for correct
             intersection
    buffDist: int, default 0, Linear distance to buffer vector data, 
            in same units as projOut. Default will return the raster cell 
            value for point data.
    fact: int, default 30, Ratio of vector area to raster cell size area. 
            Used to resample the raster to a smaller cell size
            when vector area is much smaller than raster cell size
    outND: float, default numpy NAN (np.nan), no data value to use in output
    nd_thresh: float, default 100, threshold percent of no data within vector 
            to return no data value (outND); for example, 60 means that 
            if there is greater than 60% no data within vector area, 
            all statistics will return as no data
    filenm: str, default 'outputfile', Filepath and name of output file 
            without an extension specified. Default will create a file in the 
            current working directory.
    csvout: boolean, default False, False = pickled dataframe will be created 
            as output file, True = csv file will be created as output file
    cmap: dict, default None, Dictionary of raster values (keys, as int) and 
            category names (values, as str), E.g. {1:'X', 2:'Y'}.
            Default will create a dictionary using unique raster Values.
    """


    def __init__(self, gdb, ras, lyrName=None, fldname=None, projIn=None,
                 projOut=None, buffDist=0, fact=30, outND=np.nan, 
                 nd_thresh=100, filenm='outputfile', csvout=False, cmap=None):

        self.gdb = gdb
        self.ras = ras
        self.lyrName = lyrName
        self.fldname = fldname
        self.projIn = projIn
        self.projOut = projOut
        self.buffDist = buffDist
        self.fact = fact
        self.outND = outND
        self.nd_thresh = nd_thresh
        self.filenm = filenm
        self.csvout = csvout
        self.cmap = cmap

        self.__theta = 0
        # res: float, uhhhh Leslie??
        self.__res = 1


    def extractByPoint(self):
        """
        warning - does not work on rotated rasters
        """

        src_ds = gdal.Open(self.ras, GA_ReadOnly)
        assert(rds)
        rb = rds.GetRasterBand(1)
        self.__gt = rds.GetGeoTransform()

        gdal.UseExceptions() #so it doesn't print to screen everytime point is outside grid
        # Convert from map to pixel coordinates.
        # Only works for geotransforms with no rotation.
        # If raster is rotated, see http://code.google.com/p/metageta/source/browse/trunk/metageta/geometry.py#493

        # Open feature class
        vds = ogr.Open(self.gdb, GA_ReadOnly)  
        if lyrName != None:
            vlyr = vds.GetLayerByName(self.lyrName)
        else:
            vlyr = vds.GetLayer(0)


        ptval = {}
        for feat in vlyr:   
            fldid = feat.GetFID()
            geom = feat.GetGeometryRef()
            mx = geom.GetX()
            my = geom.GetY()

            px = int((mx - self.__gt[0]) / self.__gt[1]) #x pixel
            py = int((my - self.__gt[3]) / self.__gt[5]) #y pixel
            try: #in case raster isnt full extent
                structval = rb.ReadRaster(px,py,1,1,buf_type=gdal.GDT_Float32) #Assumes 32 bit int aka 'float'
                intval = struct.unpack('f' , structval) #use the 'float' format code (8 bytes) not int (4 bytes)
                val = intval[0]
                if intval[0] < -9999:
                    val = -9999
                ptval[fldid] = val
            except:
                ptval[fldid] = self.outND
                # pass
        src_ds = None
        vlyr = None

        return(ptval)
