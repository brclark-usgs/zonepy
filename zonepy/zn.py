from __future__ import division
from __future__ import print_function
'''
Zonal Analysis

Developed by Brian Clark, Katherine Knierim, and Leslie Duncan. Portions of 
this code were modified from Copyright 2013 Matthew Perry, which were 
licensed under BSD-3 and included in this repo.


'''
import os, sys, cmath, copy
from osgeo import gdal, ogr, osr
from osgeo.gdalconst import *
import numpy as np
import pandas as pd
from scipy.ndimage import zoom
from math import radians, degrees, atan, isclose
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
    buffDist: int, default 0, Linear distance to buffer vector data, 
            in same units as vector data. Default will return the raster cell 
            value for point data.
    fact: int, default 30, Ratio of vector area to raster cell size area. 
            Used to resample the raster to a smaller cell size
            when vector area is much smaller than raster cell size
    outND: float, default numpy NAN (np.nan), no data value to use in output
            when using compute_stats() method.
    nd_thresh: float, default 100, threshold percent of no data within vector 
            to return no data value (outND); for example, 60 means that 
            if there is greater than 60% no data within vector area, 
            all statistics will return as no data
    filenm: str, default 'outputfile', Filepath and name of output file 
            without an extension specified. Default will create a file in the 
            current working directory.
    cmap: dict, default None, Dictionary of raster values (keys, as int) and 
            category names (values, as str), E.g. {1:'X', 2:'Y'} used in
            compute_category() method.
            Default will create a dictionary using unique raster Values.
    """


    def __init__(self, gdb='', ras='', lyrName=None, fldname=None,buffDist=0, 
				fact=30, outND=np.nan, nd_thresh=100, filenm='outputfile', cmap=None,
                 extractVal=None):

        self.gdb = gdb
        self.ras = ras
        self.lyrName = lyrName
        self.fldname = fldname
        self.buffDist = buffDist
        self.fact = fact
        self.outND = outND
        self.nd_thresh = nd_thresh
        self.filenm = filenm
        self.cmap = cmap
        self.extractVal = extractVal

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
        pixel_width = self.__gt[1]
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
        self.__rds = gdal.Open(self.ras, GA_ReadOnly)
        assert(self.__rds)
        #get projection info
        self.__rproj4 = self.__rds.GetProjection()
        if self.__rproj4 is None:
            print('raster has no projection, need to define projection prior to running zonepy')
        rsr = osr.SpatialReference()
        rsr.ImportFromWkt(self.__rproj4)
        self.__rproj4 = rsr.ExportToProj4()
        self.__rb = self.__rds.GetRasterBand(1)
        self.__gt = self.__rds.GetGeoTransform()
        self.__rsize = (self.__rds.RasterXSize, self.__rds.RasterYSize)
        
        # Get raster cell size
        dx1 = self.__gt[1] # pixel width
        dy1 = self.__gt[2] 
        dx2 = self.__gt[4] 
        dy2 = self.__gt[5] # pixel height
        len1 = np.sqrt(dx1**2 + dy1 **2)
        len2 = np.sqrt(dx2**2 + dy2 **2)
        self.__ras_area = len1 * len2

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
        self.__orig_nodata = self.__rb.GetNoDataValue()
        print('raster no data value', self.__orig_nodata)

    def openVector(self, extractVal=None):
        # Open feature class
        if extractVal == None:
            self.__vds = ogr.Open(self.gdb, GA_ReadOnly)  
        else:
            self.__vds = ogr.Open(self.gdb, 1)  

        if self.lyrName != None:
            self.__vlyr = self.__vds.GetLayerByName(self.lyrName)
        else:
            self.__vlyr = self.__vds.GetLayer(0)

        # get type of feature class
        self.__geomType = self.__vlyr.GetGeomType()
        # -2147483647 Point25D
        # -2147483646 LineString25D
        # -2147483645 Polygon25D
        # -2147483644 MultiPoint25D
        # -2147483643 MultiLineString25D
        # -2147483642 MultiPolygon25D
        #  0 Geometry
        #  1 Point
        #  2 Line
        #  3 Polygon
        #  4 MultiPoint
        #  5 MultiLineString
        #  6 MultiPolygon
        # 100 No Geometry

        # Deal with projections      
        vrs = self.__vlyr.GetSpatialRef()
        if vrs is None:
            print('vector has no projection, need to define projection prior to running zonepy')
            sys.exit()
        self.__vproj4 = vrs.ExportToProj4()
        self.__srcproj = osr.SpatialReference()
        self.__srcproj.ImportFromProj4(self.__vproj4)


    def checkProj(self):

        # if raster and vector projections differ,
        # create transform object
        if self.__rproj4 != self.__vproj4:
            self.__targproj = osr.SpatialReference()
            self.__targproj.ImportFromProj4(self.__rproj4)
            self.__transform = osr.CoordinateTransformation(self.__srcproj,self.__targproj)
        else:
            self.__targproj = self.__srcproj

    def vectorTest(self):
        # test if buffer is zero and geometry is point 
        # advise user to implement extractByPoint

        if self.__geomType == 1 and self.buffDist <= 0:
            print('Cannot calculate value for point with zero buffer')
            print('Consider using extractByPoint() method')
            sys.exit()

    def bufferGeom(self, feat):
        #Buffer well points, using buffDist input (in meters)
        geom = feat.GetGeometryRef()
        if self.__rproj4 != self.__vproj4:
            geom.Transform(self.__transform)
        if isinstance(self.buffDist, str):
            self.buffDist = feat.GetField(self.buffDist)

        if self.buffDist == 0 and self.__geomType == 1: 
            return(geom, 0.)
        elif self.buffDist == 0 and self.__geomType != 1:
            self.buffDist = -0.0001
            # return(geom, geom.GetArea())

        buff = geom.Buffer(self.buffDist) 

        vec_area = buff.GetArea()
        # print('Ratio Buffer to Raster: ', vec_area/ras_area)
        # verts = buff.GetGeometryRef(0)
        # verts = verts.GetPoints()

        return(buff, vec_area)

    def getField(self, feat):
        if self.fldname is None:
            self.__fldid = feat.GetFID()
        else:
            self.__fldid = feat.GetField(self.fldname)


    def computeMask(self, vec_area, src_offset, buff):
            if vec_area / self.__ras_area >= self.fact:
                zooms = 1 #no resampling needed
            else:
                zooms = int(np.sqrt(self.__ras_area / (vec_area / self.fact)))

            # src_array = rb.ReadAsArray(*src_offset)
            src_array = self.__rb.ReadAsArray(src_offset[0],src_offset[1],src_offset[2],src_offset[3])

            #Calculate new geotransform of the feature subset
            if abs(self.__theta) > 0: # if raster is rotated
                new_gt = ((self.__gt[0] + (src_offset[0] * res)), # top left x coordinate
                           res/ zooms,
                           0.,
                          (self.__gt[3] + (src_offset[1] * -res)), # top left y coordinate
                           0.,
                           -res/zooms)
            else:
                new_gt = ((self.__gt[0] + (src_offset[0] * self.__gt[1])),
                           self.__gt[1]/zooms,
                           0.0,
                          (self.__gt[3] + (src_offset[1] * self.__gt[5])),
                           0.0,
                           self.__gt[5]/zooms)
            
                # Create a temporary vector layer in memory
                mem_drv = ogr.GetDriverByName('Memory')
                mem_ds = mem_drv.CreateDataSource('out')
                mem_layer = mem_ds.CreateLayer('poly', srs=self.__srcproj, geom_type=ogr.wkbPolygon)
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
                driver = gdal.GetDriverByName('MEM')
                rvds = driver.Create('', src_offset[2]*zooms, src_offset[3]*zooms, 1, gdal.GDT_Byte)
                rvds.SetGeoTransform(new_gt)
                rvds.SetProjection(self.__srcproj.ExportToWkt())

                gdal.RasterizeLayer(rvds, [1], mem_layer, None, None, [1], ['ALL_TOUCHED=True']) #burn_values=[1])
                rv_array = rvds.ReadAsArray()

                rvds = None
                mem_ds = None
                # Resample the raster (only changes when zooms not 1)
                src_re = zoom(src_array, zooms, order = 0)
                
            return(src_re, rv_array)


    def createDF(self, flag=''):

        if self.filenm == 'outputfile':
            cwd = os.getcwd()
            self.filenm = os.path.join(cwd, 'outputfile')
        # Create dataframe from dictionary, transpose
        df = pd.DataFrame(self.__statDict)
        df = df.T
        df = df.rename(columns={self.__orig_nodata:'NoData'})
        if flag == 'zonecat':
                df = df.replace(np.nan,0) # NAN values are true zeros
        df = df.reset_index()
        cols = df.columns.tolist()
        if self.fldname is None:
            cols[0] = 'uniqueID'
        else:
            cols[0] = self.fldname
        if flag == 'zonepoint':
            cols[1] = 'value'
            cols[2] = 'long'
            cols[3] = 'lat'
        cols = [str(x) for x in cols]
        df.columns = cols
        self.df = df
        
        ## OUTPUT options
        print('\n{}'.format(self.filenm))

        return

    def compute_category(self):

        """
        Compute percent categories of raster cells within vector zone
        """
        # be sure to start with empty dict
        self.__statDict = {}

        self.openRaster()

        # Create cmap if none provided
        cmapFlag = True
        if self.cmap is None:
            self.cmap = {}
            cmapFlag = False

        print('Raster NODATA Value: ', self.__orig_nodata)
        print('Category keys and values:', self.cmap)

        self.openVector()
        self.vectorTest()
        self.checkProj()

        # Loop through vectors
        lenid = len(self.__vlyr)
        for i, feat in enumerate(self.__vlyr):
            self.getField(feat)
            sys.stdout.write('\r{} of {}, staid: {}'.format(i+1, lenid, self.__fldid))
            sys.stdout.flush()

            buff, vec_area = self.bufferGeom(feat)
            cengeom = buff.Centroid()
            mx = cengeom.GetX()
            my = cengeom.GetY()

            src_offset = self.bbox_to_pixel_offsets(buff.GetEnvelope())

            if src_offset[2] <= 0 or src_offset[3] <=0:
                #if point falls outside raster grid, include nodata as zone analysis
                masked = None
            else: 
                src_re, rv_array = self.computeMask(vec_area, src_offset, buff)


                # Mask the source data array with our current feature
                # we take the logical_not to flip 0<->1 to get the correct mask effect
                masked = np.ma.MaskedArray(src_re, mask = np.logical_not(rv_array))

            #If mask is empty, use NoData column
            if masked is None:
                pixel_count = {}
                pixel_count['NoData'] = 100
            else:
                keys, counts = np.unique(masked.compressed(), return_counts=True)
                #print(keys)
                if cmapFlag == False:
                    for k in keys:
                        self.cmap[k] = k
                else:
                    self.cmap[self.__orig_nodata] = 'NoData'

                # Create Dictionary of cmaps and counts
                pixel_count = dict(zip([self.cmap[k] for k in keys],
                                  [np.asscalar(c) for c in counts]))
                pixtot = float(sum(pixel_count.values()))
                for k, v in pixel_count.items():
                    pixel_count[k] = v / pixtot * 100
            pixel_count['long'] = mx
            pixel_count['lat'] = my

            # Create dictionary of station ids with pixel counts
            self.__statDict[self.__fldid] = pixel_count

     
        self.__vds = None
        self.__rds = None

        # Create dataframe from dictionary, transpose
        self.__flag = 'zonecat'
        self.createDF(self.__flag)


    def compute_stats(self):
        """
        Compute summary statistics of raster cells within vector zone
        """
        # be sure to start with empty dict
        self.__statDict = {}

        # Open raster
        self.openRaster() # check nodata value

        # Open feature class
        self.openVector()
        self.vectorTest()
        self.checkProj()

        # Loop through vectors
        statDict = {}
        lenid = len(self.__vlyr)
        for i, feat in enumerate(self.__vlyr):
            self.getField(feat)
            sys.stdout.write('\r{} of {}, id: {}'.format(i+1, lenid, self.__fldid))
            sys.stdout.flush()

            #Buffer well points, using buffDist input
            buff, vec_area = self.bufferGeom(feat)
            cengeom = buff.Centroid()
            mx = cengeom.GetX()
            my = cengeom.GetY()

            src_offset = self.bbox_to_pixel_offsets(buff.GetEnvelope())
            
            if src_offset[2] <= 0 or src_offset[3] <=0:
                #if point falls outside raster grid, include nodata as zone analysis
                masked = None
            else: 
                src_re, rv_array = self.computeMask(vec_area, src_offset, buff)
                # self.computeArrays(zooms) # duplicated?
                

                # Calculate the percent of No Data in masked array
                masked_nd = np.ma.MaskedArray(src_re, mask = np.logical_not(rv_array))
                # print(masked_nd)
                if self.__orig_nodata is not None:
                    if self.__orig_nodata < -3.4e38: # arc no data value
                        masked_nd[masked_nd < -3.4e38] = -9999
                    if self.__orig_nodata > 1.7e38: # QGIS no data value
                        masked_nd[masked_nd > 1.7e38] = -9999
                # print(masked_nd)

                keys, counts = np.unique(masked_nd.compressed(), return_counts=True)
                mDict = dict(zip(keys,counts))
                # print('\t')
                # print(mDict)
                # print('\t')
                # print('{:.20f}, {:.20f}'.format(keys[0], self.__orig_nodata))

                # calculate nd value in buffer
                if self.__orig_nodata in keys:
                    nd = mDict[self.__orig_nodata] / (masked_nd.shape[0] * masked_nd.shape[1]) * 100
                    # print('calculated no data: ', np.round(nd,2))
                elif -9999 in keys:
                    nd = mDict[-9999] / (masked_nd.shape[0] * masked_nd.shape[1]) * 100
                    # print('calculated no data: ', np.round(nd,2))
                else:
                    nd = 0
                                     
                # Mask the source data array with our current feature
                # we take the logical_not to flip 0<->1 to get the correct mask effect
                # we also mask out nodata values explictly
                if -9999 in keys: # to deal with ArcGIS/QGIS no data value, force no data to -9999
                    masked = np.ma.MaskedArray(src_re,
                        mask=np.logical_or(
                            src_re == -9999, # nodata_value,
                            np.logical_not(rv_array))) 
                else:
                    masked = np.ma.MaskedArray(src_re,
                        mask=np.logical_or(
                            src_re == self.__orig_nodata, # nodata_value,
                            np.logical_not(rv_array))) 
                

                feature_stats = {
                    'min': float(np.ma.masked_invalid(masked).min()),
                    'mean': float(np.ma.masked_invalid(masked).mean()),
                    'max': float(np.ma.masked_invalid(masked).max()),
                    'std': float(np.ma.masked_invalid(masked).std()),
                    'sum': float(np.ma.masked_invalid(masked).sum()),
                    'count': int(np.ma.masked_invalid(masked).count()),
                    'median': float(np.ma.median(np.ma.masked_invalid(masked))),
                    'long': mx,
                    'lat': my}


            no_stats = {
                'min': self.outND,
                'mean': self.outND,
                'max': self.outND,
                'std': self.outND,
                'sum': self.outND,
                'count': self.outND,
                'median': self.outND,
                'long': mx,
                'lat': my}

            # print('no data percent: {}, no data threshold: {}\n'.format(nd,nd_thresh))
            if masked is not None:
                # print('\t')
                # print('calculating stats...')
                if np.isnan(float(np.ma.masked_invalid(masked).mean())):
                    # print('all nan')
                    self.__statDict[self.__fldid] = no_stats # if all NAN, return nodata value
                else:
                    # print('no data percent: {}, no data threshold: {}\n'.format(np.round(nd,2),self.nd_thresh))
                    if nd >= self.nd_thresh: # insufficient data, return nodata value
                        self.__statDict[self.__fldid] = no_stats
                        # print('insufficient data')
                    else: # sufficient data, return stats
                        self.__statDict[self.__fldid] = feature_stats
                        # print('stats calculated')
                        # print(self.__statDict[self.__fldid])
            else:
                # print('\t')
                # print('point outside of raster extent')
                self.__statDict[self.__fldid] = no_stats # if outside of raster extent, return nodata value

        #clearing memory
        self.__vds = None
        self.__rds = None

        ##OUTPUT
        self.__flag = 'zonestat'
        self.createDF(self.__flag)

        

    def extractByPoint(self, extractVal='extractVal'):
        """
        warning - does not work on rotated rasters

        extractVal: string, default None, column name for new field if extractByPoint
                  method is used
        """
        import struct

        gdal.UseExceptions() #so it doesn't print to screen everytime point is outside grid
        # Convert from map to pixel coordinates.
        # Only works for geotransforms with no rotation.
        # If raster is rotated, see http://code.google.com/p/metageta/source/browse/trunk/metageta/geometry.py#493

        self.openRaster()

        # Open feature class
        # self.openVector(extractVal=extractVal)
        self.openVector()
        self.checkProj()

        lyrDef = self.__vlyr.GetLayerDefn()
        fieldExists = False
        for i in range(lyrDef.GetFieldCount()):
            if lyrDef.GetFieldDefn(i).GetName() == extractVal:
                fieldExists = True
                break
        if not fieldExists:
            fieldDef = ogr.FieldDefn(extractVal, ogr.OFTReal)
            self.__vlyr.CreateField(fieldDef)

        self.__statDict = {}
        for feat in self.__vlyr: 
            self.getField(feat)  

            buff, vec_area = self.bufferGeom(feat)
            mx = buff.GetX()
            my = buff.GetY()

            px = int((mx - self.__gt[0]) / self.__gt[1]) #x pixel
            py = int((my - self.__gt[3]) / self.__gt[5]) #y pixel
            try: #in case raster isnt full extent
                structval = self.__rb.ReadRaster(px,py,1,1,buf_type=gdal.GDT_Float32) #Assumes 32 bit int aka 'float'
                intval = struct.unpack('f' , structval) #use the 'float' format code (8 bytes) not int (4 bytes)
                val = intval[0]
                if intval[0] < -9999:
                    val = self.outND
            except:
                val = self.outND
            # print(self.__fldid, val)
            if self.__orig_nodata is not None:
                if val < -3.4e38: # arc no data value
                    val = self.outND
                if val > 1.7e38: # QGIS no data value
                    val = self.outND
                    # print('new val!', val)
                if val == self.__orig_nodata:
                    val = self.outND
                    # print('new val!', val)

            self.__statDict[self.__fldid] = (val, mx, my)
        # print('no data',self.__orig_nodata)  
        self.__vds = None
        self.__rds = None

        #Output
        self.createDF('zonepoint')


    def writeRaster(self, stat='mean', outTiff='outras.tif',
                  inputfile=None):
        '''
        Define numpy array of raster shape and write raster

        Parameters
        ----------
        stat: str, default = 'mean', Summary statistic to use for raster creation.
                Options are 'min' (minimum), 'mean', 'max' (maximum), 'std' 
                (standard deviation), 'sum', 'count', 'median'
        outTiff: str, Filepath and name of created raster
        inputfile: str, filename of pickle or csv from compute_stats. if None, it
                 uses self.__statDict
        '''

        # Open feature class
        self.openVector()

        # Get extent of single polygon
        feat = self.__vlyr.GetNextFeature()
        featext = feat.GetGeometryRef().GetEnvelope()
        # print(featext)
        xsz = featext[1] - featext[0] #resolution in x
        ysz = featext[3] - featext[2] #resolution in y
        sz = (xsz, ysz)
        # print(sz)

        # Get projection and extent
        vsr = self.__vlyr.GetSpatialRef()
        # print(vsr)
        vext = self.__vlyr.GetExtent()
        cols = int(round(((vext[1] - vext[0]) / xsz),0))
        rows = int(round(((vext[3] - vext[2]) / ysz),0))
        #print(rows,cols)

        if inputfile == None:
            df = self.df
        elif inputfile.endswith('pkl'):
            # Read pickle file of summary stats
            df = pd.read_pickle(pklfile)
        elif inputfile.endswith('csv'):
            df = pd.read_csv(inputfile)

        if self.__flag == 'zonestat':
            a = df[stat].values

            # Reshape to zone shape
            a = a.reshape((rows,cols))

            # Write raster
            self.array2Raster(a, vext, sz, outTiff=outTiff)
        else: #category
            dfcols = df.columns.tolist()
            dfcols = dfcols[1:-2]
            for c in dfcols:
                a = df[c].values
                a = a.reshape((rows,cols))
                outname = os.path.splitext(outTiff)[0] + '{}.tif'.format(c)
                self.array2Raster(a, vext, sz, outTiff=outname)


        # Clear memory
        self.__vds = None



    def array2Raster(self, grid, ext, sz, outTiff='outras.tif'):
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
        outTiff: str, default 'outras.tif', Filepath of raster to be created.
            Default creates file in working directory.
        '''

        print('writing tiff...')

        # Calculate angle and coordinates
        angle = 0
        angle = angle * -1. * np.pi / 180.
        geotransform = (ext[0], sz[0] * np.cos(angle), -sz[0] * np.sin(angle),
                        ext[3], -sz[1] * np.sin(angle), -sz[1] * np.cos(angle)) 

        # Create raster to hold data
        nrow = grid.shape[0]
        ncol = grid.shape[1]
        drv = gdal.GetDriverByName('GTiff')
        ds = drv.Create(outTiff, ncol, nrow, 1 ,gdal.GDT_Float32)

        # Set raster no data value
        band = ds.GetRasterBand(1)
        band.SetNoDataValue(self.outND)

        # Set coordinates and projection
        ds.SetGeoTransform(geotransform)  
        ds.SetProjection(self.__srcproj.ExportToWkt())

        # Write array to raster
        band.WriteArray(grid)

        ds = None

    def writeCSV(self):

        self.df.to_csv('{}.csv'.format(self.filenm), index=False)

    def writePKL(self):

        self.df.to_pickle('{}.pkl'.format(self.filenm))

    def writeSHP(self):
        shp = os.path.join('{}.shp'.format(self.filenm))
        #make new shpfile here
        drv = ogr.GetDriverByName('ESRI Shapefile')
        if os.path.exists(shp):
            drv.DeleteDataSource(shp)
        ds = drv.CreateDataSource(shp) 
        # lyr = ds.CreateLayer('lyr', geom_type=ogr.wkbPoint)
        self.df2features(ds)

    def writeGPKG(self):
        shp = os.path.join('{}.gpkg'.format(self.filenm))
        #make new shpfile here
        drv = ogr.GetDriverByName('GPKG')
        # if os.path.exists(shp):
            # drv.DeleteDataSource(shp)
        ds = drv.CreateDataSource(shp) 
        self.df2features(ds)

    def df2features(self, ds=''):
        lyr = ds.CreateLayer('lyr', geom_type=ogr.wkbPoint, srs=self.__srcproj)
        newtrns = osr.CoordinateTransformation(self.__targproj, self.__srcproj)
        cols = self.df.columns.tolist()
        cols = [str(x) for x in cols]

        for c in cols:
            # print(c)
            if self.df[c].dtype == np.float64:
                ogrt = ogr.OFTReal
                # print('real')
            elif self.df[c].dtype == np.int:
                ogrt = ogr.OFTInteger64
                self.df[c] = self.df.astype({c:np.int64})
            else:
                ogrt = ogr.OFTString
                # print('str')
            fieldDef = ogr.FieldDefn(str(c), ogrt)
            lyr.CreateField(fieldDef)

        for row in self.df.itertuples():
            # print(row)
            pt = ogr.Geometry(ogr.wkbPoint)
            pt.AddPoint(row.long, row.lat)
            pt.Transform(newtrns)
            feat = ogr.Feature(lyr.GetLayerDefn())
            for i, v in enumerate(row[1:]):
                feat.SetField(cols[i], v)
            feat.SetGeometry(pt)
            lyr.CreateFeature(feat)

        lyr = None
        ds = None
        drv = None