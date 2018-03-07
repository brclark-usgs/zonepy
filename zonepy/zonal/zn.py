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
from scipy.ndimage import zoom
from math import radians, degrees, atan
import warnings
# gdal.PushErrorHandler('CPLQuietErrorHandler')
 
class ZoneClass(object):
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


    def __init__(self, gdb='', ras='', lyrName=None, fldname=None, projIn=None,
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

        return

    def bbox_to_pixel_offsets(self, bbox):
        """
        Define raster cells within polygon zone

        Parameters
        ----------
        gt: geotransform of raster
        bbox: bounding box from vector (polygon or point buffer)

        rsize: size of raster (total rows and columns)

        """

        originX = self.__gt[0]
        originY = self.__gt[3]
        # pixel_width = self.__gt[1]
        pixel_height = self.__gt[5]

        if abs(self.__theta) > 0: # if raster is rotated
            self.__theta = self.__theta * -1 
            vx1,vy1 = rotatePt(bbox[0][0], bbox[0][1], self.__gt, self.__theta)
            vx2,vy2 = rotatePt(bbox[1][0], bbox[1][1], self.__gt, self.__theta)
            vx3,vy3 = rotatePt(bbox[2][0], bbox[2][1], self.__gt, self.__theta)
            vx4,vy4 = rotatePt(bbox[3][0], bbox[3][1], self.__gt, self.__theta)
            wkt = 'LINESTRING ({} {}, {} {}, {} {}, {} {})'.format(vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4)
            geom = ogr.CreateGeometryFromWkt(wkt)
            env = geom.GetEnvelope()
            x1 = int((env[0] - originX) / self.__res)
            x2 = int((env[1] - originX) / self.__res) + 1
            y1 = int((env[3] - originY) / -self.__res)
            y2 = int((env[2] - originY) / -self.__res) + 1
        else: 
            x1 = int((bbox[0] - originX) / pixel_width)
            x2 = int((bbox[1] - originX) / pixel_width) + 1
            y1 = int((bbox[3] - originY) / pixel_height)
            y2 = int((bbox[2] - originY) / pixel_height) + 1
        
        # "Clip" the geometry bounds to the overall raster bounding box
        # This should avoid any rasterIO errors for partially overlapping polys
        if x1 < 0:
            x1 = 0
        if x2 > self.__rsize[0]:
            x2 = self.__rsize[0]
        if y1 < 0:
            y1 = 0
        if y2 > self.__rsize[1]:
            y2 = self.__rsize[1]

        xsize = x2 - x1
        ysize = y2 - y1

        # print('src_offset = {},{},{},{}'.format(x1, y1, xsize, ysize))
        return(x1,y1,xsize,ysize)


    def rotatePt(self, dX, dY, gt):
        """
        Rotate coordinates to user specified angle
        
        Parameters
        ----------
        dX: float, original x coordinate
        dY: float, original y coordinate
        gt: geotransform object
        theta: degree, angle to rotate coordinates
        
        """
        pol = cmath.polar(complex(dX-gt[0],dY-gt[3]))
        newTheta = radians(self.__theta) + pol[1]
        newxy = cmath.rect(pol[0],newTheta)
        newx = newxy.real + gt[0]
        newy = newxy.imag + gt[3]

        return(newx,newy)
     
    def openRaster(self):
        # Open raster
        rds = gdal.Open(self.ras, GA_ReadOnly)
        assert(rds)
        rb = rds.GetRasterBand(1)
        self.__gt = rds.GetGeoTransform()
        self.__rsize = (rds.RasterXSize, rds.RasterYSize)
        
        # Get raster cell size
        dx1 = self.__gt[1] # pixel width
        dy1 = self.__gt[2] 
        dx2 = self.__gt[4] 
        dy2 = self.__gt[5] # pixel height
        len1 = np.sqrt(dx1**2 + dy1 **2)
        len2 = np.sqrt(dx2**2 + dy2 **2)
        ras_area = len1 * len2

        #if needed, angle of rotation (rotated rasters)
        xorig = self.__gt[0] 
        yorig = self.__gt[3]
        self.__theta = degrees(atan(dy1/dx1)) # atan returns values in the interval [-pi/2,pi/2] 
                                       # (or values in the range[-90, 90])
        if yorig < xorig: 
            self.__theta = 90 - (self.__theta * -1)
        rads = self.__theta * -1 * np.pi/180.
        if rads == 0.:
            rads = 1
        res = self.__gt[2]/-np.sin(rads)

        #Get NoData value
        self.__orig_nodata = rb.GetNoDataValue()

    def openVector(self):
        # Open feature class
        vds = ogr.Open(self.gdb, GA_ReadOnly)  
        if self.lyrName != None:
            self.__vlyr = vds.GetLayerByName(self.lyrName)
        else:
            self.__vlyr = vds.GetLayer(0)

        # Create memory drivers to hold arrays
        mem_drv = ogr.GetDriverByName('Memory')
        driver = gdal.GetDriverByName('MEM')

        # Deal with projections
        if self.projIn == None:
            self.projIn = self.__vlyr.GetSpatialRef().ExportToProj4()

        srcproj = osr.SpatialReference()
        srcproj.ImportFromProj4(self.projIn)

        if self.projOut != None:
            targproj = osr.SpatialReference()
            targproj.ImportFromProj4(self.projOut)
            self.__transform = osr.CoordinateTransformation(srcproj,targproj)
     
    def bufferGeom(self):
        #Buffer well points, using buffDist input (in meters)
        geom = feat.GetGeometryRef()
        if self.projOut != None:
            geom.Transform(self.__transform)
        if isinstance(self.buffDist, str):
            self.buffDist = feat.GetField(self.buffDist)

        buff = geom.Buffer(self.buffDist) 
        vec_area = buff.GetArea()
        # print('Ratio Buffer to Raster: ', vec_area/ras_area)
        verts = buff.GetGeometryRef(0)
        verts = verts.GetPoints()

    def getField(self):
        if self.fldname is None:
            self.__fldid = feat.GetFID()
        else:
            self.__fldid = feat.GetField(self.fldname)


    def computeMask(self):
            if vec_area / ras_area >= self.fact:
                self.__zooms = 1 #no resampling needed
            else:
                self.__zooms = int(np.sqrt(ras_area / (vec_area / self.fact)))

            # src_array = rb.ReadAsArray(*src_offset)
            src_array = rb.ReadAsArray(src_offset[0],src_offset[1],src_offset[2],src_offset[3])

            #Calculate new geotransform of the feature subset
            if abs(self.__theta) > 0: # if raster is rotated
                self.__new_gt = ((self.__gt[0] + (src_offset[0] * res)), # top left x coordinate
                           res/ zooms,
                           0.,
                          (self.__gt[3] + (src_offset[1] * -res)), # top left y coordinate
                           0.,
                           -res/zooms)
            else:
                self.__new_gt = ((self.__gt[0] + (src_offset[0] * self.__gt[1])),
                           self.__gt[1]/zooms,
                           0.0,
                          (self.__gt[3] + (src_offset[1] * self.__gt[5])),
                           0.0,
                           self.__gt[5]/zooms)
            
                            # Create a temporary vector layer in memory
                mem_ds = mem_drv.CreateDataSource('out')
                mem_layer = mem_ds.CreateLayer('poly', None, ogr.wkbPolygon)
                mem_poly = ogr.Feature(mem_layer.GetLayerDefn())
                if abs(self.__theta) > 0: # if raster is rotated
                    ring = ogr.Geometry(ogr.wkbLinearRing)
                    for v in verts:
                        x,y = rotatePt(v[0],v[1],rgt,-theta)
                        ring.AddPoint(x, y)
                    ply = ogr.Geometry(ogr.wkbPolygon25D)
                    ply.AddGeometry(ring)
                    mem_poly.SetGeometryDirectly(ply)
                else:
                    mem_poly.SetGeometryDirectly(buff)
                mem_layer.CreateFeature(mem_poly)

                # Rasterize the polygon
                rvds = driver.Create('', src_offset[2]*zooms, src_offset[3]*zooms, 1, gdal.GDT_Byte)
                rvds.SetGeoTransform(new_gt)
                gdal.RasterizeLayer(rvds, [1], mem_layer, None, None, [1], ['ALL_TOUCHED=True']) #burn_values=[1])
                rv_array = rvds.ReadAsArray()

                # Resample the raster (only changes when zooms not 1)
                src_re = zoom(src_array, zooms, order = 0)

    def computeArrays(self):
        # Create a temporary vector layer in memory
        mem_ds = mem_drv.CreateDataSource('out')
        mem_layer = mem_ds.CreateLayer('poly', None, ogr.wkbPolygon)
        mem_poly = ogr.Feature(mem_layer.GetLayerDefn())
        mem_poly.SetGeometryDirectly(buff)
        mem_layer.CreateFeature(mem_poly)
        
        # Rasterize the polygon
        rvds = driver.Create('', src_offset[2]*zooms, src_offset[3]*zooms, 1, gdal.GDT_Byte)
        rvds.SetGeoTransform(new_gt)
        gdal.RasterizeLayer(rvds, [1], mem_layer, burn_values=[1])
        rv_array = rvds.ReadAsArray()

        # Resample the raster (only changes when zooms not 1)
        src_re = zoom(src_array, zooms, order = 0)

    def createOutput(self):
        # Create dataframe from dictionary, transpose
        df = pd.DataFrame(self.__statDict)
        df = df.T
        df = df.replace(np.nan,0) # NAN values are true zeros
        df = df.reset_index()
        cols = df.columns.tolist()
        if self.fldname is None:
            cols[0] = 'uniqueID'
        else:
            cols[0] = self.fldname
        df.columns = cols

        ## OUTPUT options
        if self.filenm == 'outputfile':
            cwd = os.getcwd()
            self.filenm = os.path.join(cwd, 'outputfile')
        if self.csvout == True:
            print(self.filenm)
            df.to_csv('{}.csv'.format(self.filenm), index=False)
        else:
            df.to_pickle('{}.pkl'.format(self.filenm))

    def compute_category(self):

        """
        Compute percent categories of raster cells within vector zone
        """
        self.openRaster()

        # Create cmap if none provided
        cmapFlag = True
        if self.cmap is None:
            self.cmap = {}
            cmapFlag = False

        print('Raster NODATA Value: ', orig_nodata)
        print('Category keys and values:', self.cmap)

        self.openVector()

        # Loop through vectors
        stats = []
        statDict = {}
        lenid = len(vlyr)
        for i, feat in enumerate(vlyr):
            self.getField()

            sys.stdout.write('\r{} of {}, staid: {}'.format(i+1, lenid, fldid))
            sys.stdout.flush()

            self.bufferGeom()

            src_offset = self.bbox_to_pixel_offsets(buff.GetEnvelope())

            if src_offset[2] <= 0 or src_offset[3] <=0:
                #if point falls outside raster grid, include nodata as zone analysis
                self.__masked = None
            else: 
                self.computeMask()

                self.computeArrays()

                # Mask the source data array with our current feature
                # we take the logical_not to flip 0<->1 to get the correct mask effect
                masked = np.ma.MaskedArray(src_re, mask = np.logical_not(rv_array))

            #If mask is empty, use NoData column
            if masked is None:
                pixel_count = {}
                pixel_count['NoData'] = 100
            else:
                keys, counts = np.unique(masked.compressed(), return_counts=True)
                if cmapFlag == False:
                    for k in keys:
                        self.cmap[k] = k
                else:
                    self.cmap[orig_nodata] = 'NoData'

                # Create Dictionary of cmaps and counts
                pixel_count = dict(zip([self.cmap[k] for k in keys],
                                  [np.asscalar(c) for c in counts]))
                pixtot = float(sum(pixel_count.values()))
                for k, v in pixel_count.items():
                    pixel_count[k] = v / pixtot * 100

            # Create dictionary of station ids with pixel counts
            statDict[fldid] = pixel_count

        #clearing memory
            rvds = None
            mem_ds = None
     
        vds = None
        rds = None

        # Create dataframe from dictionary, transpose
        self.createOutput()


    def compute_stats(self):
        """
        Compute summary statistics of raster cells within vector zone
        """

        # Open raster
        self.openRaster() # check nodata value

        # Open feature class
        self.openVector()

        # Loop through vectors
        statDict = {}
        lenid = len(self.__vlyr)
        for i, feat in enumerate(self.__vlyr):
            getField()
            sys.stdout.write('\r{} of {}, id: {}\n'.format(i+1, lenid, fldid))
            sys.stdout.flush()

            #Buffer well points, using buffDist input
            self.bufferGeom()
            
            src_offset = self.bbox_to_pixel_offsets(buff.GetEnvelope())

            if src_offset[2] <= 0 or src_offset[3] <=0:
                #if point falls outside raster grid, include nodata as zone analysis
                self.__masked = None
            else: 
                self.computeMask()

                self.computeArrays()
                
                # Mask the source data array with our current feature
                # we take the logical_not to flip 0<->1 to get the correct mask effect
                # we also mask out nodata values explictly
                masked = np.ma.MaskedArray(src_re,
                    mask=np.logical_or(
                        src_re == nodata_value,
                        np.logical_not(rv_array))) 

                # Calculate the percent of No Data in masked array
                nd = 0
                masked_nd = np.ma.MaskedArray(src_re, mask = np.logical_not(rv_array))
                keys, counts = np.unique(masked_nd.compressed(), return_counts=True)
                mDict = dict(zip(keys,counts))
                if nodata_value in keys:
                    nd = mDict[nodata_value] / (masked_nd.shape[0] * masked_nd.shape[1]) * 100


                feature_stats = {
                    'min': float(np.ma.masked_invalid(masked).min()),
                    'mean': float(np.ma.masked_invalid(masked).mean()),
                    'max': float(np.ma.masked_invalid(masked).max()),
                    'std': float(np.ma.masked_invalid(masked).std()),
                    'sum': float(np.ma.masked_invalid(masked).sum()),
                    'count': int(np.ma.masked_invalid(masked).count()),
                    'median': float(np.ma.median(np.ma.masked_invalid(masked)))}

            no_stats = {
                'min': -9999.,
                'mean': -9999.,
                'max': -9999.,
                'std': -9999.,
                'sum': -9999.,
                'count': -9999.,
                'median': -9999.}

            # print('no data percent: {}, no data threshold: {}\n'.format(nd,nd_thresh))
            if masked is not None:
                if np.isnan(float(np.ma.masked_invalid(masked).mean())):
                    self.__statDict[fldid] = no_stats # if all NAN, return -9999
                else:
                    if nd >= self.nd_thresh: # insufficient data, return -9999
                        self.__statDict[fldid] = no_stats
                    else: # sufficient data, return stats
                        self.__statDict[fldid] = feature_stats
            else:
                self.__statDict[fldid] = no_stats # if outside of raster extent, return -9999

        #clearing memory
            rvds = None
            mem_ds = None
        vds = None
        rds = None

        ##OUTPUT
        self.createOutput()