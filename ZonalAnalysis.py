from __future__ import division
from __future__ import print_function
'''
Zonal Analysis

Developed by Brian Clark, Katherine Knierim, and Leslie Duncan. Portions of 
this code were modified from Copyright 2013 Matthew Perry, which were 
licensed under BSD-3 and included in this repo.



'''
import os,sys,cmath
from osgeo import gdal, ogr, osr
from osgeo.gdalconst import *
import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import zoom
from math import radians, degrees, atan
import matplotlib.pyplot as plt
import warnings
gdal.PushErrorHandler('CPLQuietErrorHandler')
 

def bbox_to_pixel_offsets(gt, bbox, rsize, theta=0, res=1):
	"""
	Define raster cells within polygon zone

	Parameters
	----------
	gt: geotransform of raster
	bbox: bounding box from vector (polygon or point buffer)
	rsize: size of raster (total rows and columns)
	theta: degree, angle to rotate coordinates
	res: float, uhhhh Leslie??

    """

	originX = gt[0]
	originY = gt[3]
	pixel_width = gt[1]
	pixel_height = gt[5]

	if abs(theta) > 0: # if raster is rotated
		theta = theta * -1 
		vx1,vy1 = rotatePt(bbox[0][0],bbox[0][1],gt,theta)
		vx2,vy2 = rotatePt(bbox[1][0],bbox[1][1],gt,theta)
		vx3,vy3 = rotatePt(bbox[2][0],bbox[2][1],gt,theta)
		vx4,vy4 = rotatePt(bbox[3][0],bbox[3][1],gt,theta)
		wkt = 'LINESTRING ({} {}, {} {}, {} {}, {} {})'.format(vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4)
		geom = ogr.CreateGeometryFromWkt(wkt)
		env = geom.GetEnvelope()
		x1 = int((env[0] - originX) / res)
		x2 = int((env[1] - originX) / res) + 1
		y1 = int((env[3] - originY) / -res)
		y2 = int((env[2] - originY) / -res) + 1
	else: 
		x1 = int((bbox[0] - originX) / pixel_width)
		x2 = int((bbox[1] - originX) / pixel_width) + 1
		y1 = int((bbox[3] - originY) / pixel_height)
		y2 = int((bbox[2] - originY) / pixel_height) + 1
	
	# "Clip" the geometry bounds to the overall raster bounding box
	# This should avoid any rasterIO errors for partially overlapping polys
	if x1 < 0:
		x1 = 0
	if x2 > rsize[0]:
		x2 = rsize[0]
	if y1 < 0:
		y1 = 0
	if y2 > rsize[1]:
		y2 = rsize[1]

	xsize = x2 - x1
	ysize = y2 - y1

	# print('src_offset = {},{},{},{}'.format(x1, y1, xsize, ysize))
	return(x1,y1,xsize,ysize)

def rotatePt(dX,dY,gt,theta):
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
	newTheta = radians(theta) + pol[1]
	newxy = cmath.rect(pol[0],newTheta)
	newx = newxy.real + gt[0]
	newy = newxy.imag + gt[3]

	return(newx,newy)
 

def zonal_stats(gdb, ras, lyrName=None, fldname=None , 
	projIn=None, projOut=None, buffDist=0, fact=30, 
	outND=np.nan, nd_thresh=100, filenm='outputfile', csvout=False):
	"""
	Compute summary statistics of raster cells within vector zone

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
	"""
	# Open raster
	rds = gdal.Open(ras, GA_ReadOnly)
	assert(rds)
	rb = rds.GetRasterBand(1)
	rgt = rds.GetGeoTransform()
	rsize = (rds.RasterXSize, rds.RasterYSize)

	# Get raster cell size
	dx1 = rgt[1] # pixel width
	dy1 = rgt[2] 
	dx2 = rgt[4] 
	dy2 = rgt[5] # pixel height
	len1 = np.sqrt(dx1**2 + dy1 **2)
	len2 = np.sqrt(dx2**2 + dy2 **2)
	ras_area = len1 * len2

	#if needed, angle of rotation (rotated rasters)
	xorig = rgt[0] 
	yorig = rgt[3]
	theta = degrees(atan(dy1/dx1)) # atan returns values in the interval [-pi/2,pi/2] 
								   # (or values in the range[-90, 90])
	if yorig < xorig: 
		theta = 90 - (theta * -1)
	rads = theta * -1 * np.pi/180.
	res = rgt[2]/-np.sin(rads)

   # Get NoData value
	nodata_value = rb.GetNoDataValue()
	print('Raster NODATA Value: ', nodata_value)

	# Open feature class
	vds = ogr.Open(gdb, GA_ReadOnly)  
	if lyrName != None:
		vlyr = vds.GetLayerByName(lyrName)
	else:
		vlyr = vds.GetLayer(0)

	#Create memory drivers to hold arrays
	mem_drv = ogr.GetDriverByName('Memory')
	driver = gdal.GetDriverByName('MEM')
 
	# Deal with projection
	if projIn == None:
		projIn = vlyr.GetSpatialRef().ExportToProj4()

	srcproj = osr.SpatialReference()
	srcproj.ImportFromProj4(projIn)

	if projOut != None:
		targproj = osr.SpatialReference()
		targproj.ImportFromProj4(projOut)
		transform = osr.CoordinateTransformation(srcproj,targproj)

	# Loop through vectors
	statDict = {}
	lenid = len(vlyr)
	for i, feat in enumerate(vlyr):
		if fldname is None:
			fldid = feat.GetFID()
			lyrdef = vlyr.GetLayerDefn()
			fldname = lyrdef.GetFieldDefn(0).GetName()
		else:
			fldid = feat.GetField(fldname)

		fldid = feat.GetField(fldname)

		sys.stdout.write('\r{} of {}, staid: {}\n'.format(i+1, lenid, fldid))
		sys.stdout.flush()

		#Buffer well points, using buffDist input
		geom = feat.GetGeometryRef()
		if projOut != None:
			geom.Transform(transform)
		buff = geom.Buffer(buffDist) 
		# geom = geom.Buffer(buffDist) 
		# buff = geom
		vec_area = buff.GetArea()
		# verts = geom.GetGeometryRef(0)
		verts = buff.GetGeometryRef(0)
		verts = verts.GetPoints()
		# print('Ratio Buffer to Raster: ', vec_area/ras_area)
		
		# if abs(theta) > 0:
		# 	src_offset = bbox_to_pixel_offsets(rgt,verts,rsize,theta,res)
		# else:
		src_offset = bbox_to_pixel_offsets(rgt, buff.GetEnvelope(), rsize, theta, res)

		if src_offset[2] <= 0 or src_offset[3] <=0:
			#if point falls outside raster grid, include nodata as zone analysis
			masked = None
		else:    
			if vec_area/ras_area >= fact:
				zooms = 1 #no resampling needed
			else:
				zooms = int(np.sqrt(ras_area/(vec_area/fact)))

			# src_array = rb.ReadAsArray(*src_offset)
			src_array = rb.ReadAsArray(src_offset[0],src_offset[1],src_offset[2],src_offset[3])
			# print('src_array = \n{}'.format(src_array))
			# plt.imshow(src_array)
			# plt.colorbar()
			# plt.show()

			#Calculate new geotransform of the feature subset
			if abs(theta) > 0: # if raster is rotated
				new_gt = ((rgt[0] + (src_offset[0] * res)), # top left x coordinate
						   res/ zooms,
						   0.,
						  (rgt[3] + (src_offset[1] * -res)), # top left y coordinate
						   0.,
						   -res/zooms)
			else:
				new_gt = ((rgt[0] + (src_offset[0] * rgt[1])),
						   rgt[1]/zooms,
						   0.0,
						  (rgt[3] + (src_offset[1] * rgt[5])),
						   0.0,
						   rgt[5]/zooms)

			# Create a temporary vector layer in memory
			mem_ds = mem_drv.CreateDataSource('out')
			mem_layer = mem_ds.CreateLayer('poly', None, ogr.wkbPolygon)
			mem_poly = ogr.Feature(mem_layer.GetLayerDefn())
			if abs(theta) > 0: # if raster is rotated
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
			# print('rv_array = \n{}'.format(rv_array)) 
			# plt.imshow(rv_array)
			# plt.colorbar()
			# plt.show()
			# print('\n{}'.format(rv_array.shape))

			# Resample the raster (only changes when zooms not 1)
			src_re = zoom(src_array, zooms, order = 0)
			# print('rc_re = \n{}'.format(src_re))
			# plt.imshow(src_re)
			# plt.show()
			# print('\n{}'.format(src_re.shape))
			
			# Mask the source data array with our current feature
			# we take the logical_not to flip 0<->1 to get the correct mask effect
			# we also mask out nodata values explictly
			masked = np.ma.MaskedArray(src_re,
				mask=np.logical_or(
					src_re == nodata_value,
					np.logical_not(rv_array))) 
			# print('masked = \n{}'.format(masked))
			# plt.imshow(masked)
			# plt.show()
			# print('\n{}'.format(masked.shape))

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
				statDict[feat.GetField(fldname)] = no_stats # if all NAN, return -9999
			else:
				if nd >= nd_thresh: # insufficient data, return -9999
					statDict[feat.GetField(fldname)] = no_stats
				else: # sufficient data, return stats
					statDict[feat.GetField(fldname)] = feature_stats
		else:
			statDict[feat.GetField(fldname)] = no_stats # if outside of raster extent, return -9999

	#clearing memory
		rvds = None
		mem_ds = None
	vds = None
	rds = None

	##OUTPUT
	df = pd.DataFrame(statDict)
	df = df.T
	df = df.reset_index()
	cols = df.columns.tolist()
	cols[0] = fldname
	df.columns = cols
	# print(df)

	## OUTPUT options
	if csvout == True:
		df.to_csv('{}.csv'.format(filenm), index=False)
	else:
		df.to_pickle('{}.pkl'.format(filenm))




def zonal_category(gdb, ras, lyrName=None, fldname=None, projIn=None, 
				projOut=None, buffDist=0, cmap=None, fact=30, filenm='outputfile',
				csvout=False):
	"""
	Compute percent categories of raster cells within vector zone

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
	cmap: dict, default None, Dictionary of raster values (keys, as int) and 
			category names (values, as str), E.g. {1:'X', 2:'Y'}.
			Default will create a dictionary using unique raster Values.
	fact: int, default 30, Ratio of vector area to raster cell size area. 
			Used to resample the raster to a smaller cell size
			when vector area is much smaller than raster cell size
	outND: float, default , No Data value to use as output
	filenm: str, default 'outputfile', Filepath and name of output file 
			without an extension specified. Default will create a file in the 
			current working directory.
	csvout: boolean, default False, False = pickled dataframe will be created 
			as output file, True = csv file will be created as output file
	"""
	# Open raster
	rds = gdal.Open(ras, GA_ReadOnly)
	assert(rds)
	rb = rds.GetRasterBand(1)
	rgt = rds.GetGeoTransform()
	rsize = (rds.RasterXSize, rds.RasterYSize)
	
	# Get raster cell size
	dx1 = rgt[1] # pixel width
	dy1 = rgt[2] 
	dx2 = rgt[4] 
	dy2 = rgt[5] # pixel height
	len1 = np.sqrt(dx1**2 + dy1 **2)
	len2 = np.sqrt(dx2**2 + dy2 **2)
	ras_area = len1 * len2

	#if needed, angle of rotation (rotated rasters)
	xorig = rgt[0] 
	yorig = rgt[3]
	theta = degrees(atan(dy1/dx1)) # atan returns values in the interval [-pi/2,pi/2] 
								   # (or values in the range[-90, 90])
	if yorig < xorig: 
		theta = 90 - (theta * -1)
	rads = theta * -1 * np.pi/180.
	res = rgt[2]/-np.sin(rads)

	#Get NoData value
	orig_nodata = rb.GetNoDataValue()

	# Create cmap if none provided
	if cmap is None:
		vals = np.unique(rb.ReadAsArray())
		cmap = {v:v for v in vals}
	else:
		cmap[orig_nodata] = 'NoData'
	print('Raster NODATA Value: ', orig_nodata)
	print('Category keys and values:', cmap)

	# Open feature class
	vds = ogr.Open(gdb, GA_ReadOnly)  
	if lyrName != None:
		vlyr = vds.GetLayerByName(lyrName)
	else:
		vlyr = vds.GetLayer(0)

	# Create memory drivers to hold arrays
	mem_drv = ogr.GetDriverByName('Memory')
	driver = gdal.GetDriverByName('MEM')

	# Deal with projections
	if projIn == None:
		projIn = vlyr.GetSpatialRef().ExportToProj4()

	srcproj = osr.SpatialReference()
	srcproj.ImportFromProj4(projIn)

	if projOut != None:
		targproj = osr.SpatialReference()
		targproj.ImportFromProj4(projOut)
		transform = osr.CoordinateTransformation(srcproj,targproj)

	# Loop through vectors
	stats = []
	statDict = {}
	lenid = len(vlyr)
	for i, feat in enumerate(vlyr):
		if fldname is None:
			fldid = feat.GetFID()
			lyrdef = vlyr.GetLayerDefn()
			fldname = lyrdef.GetFieldDefn(0).GetName()
		else:
			fldid = feat.GetField(fldname)
		sys.stdout.write('\r{} of {}, staid: {}'.format(i+1, lenid, fldid))
		sys.stdout.flush()

		#Buffer well points, using buffDist input (in meters)
		geom = feat.GetGeometryRef()
		if projOut != None:
			geom.Transform(transform)
		if isinstance(buffDist, str):
			buffDist = feat.GetField(buffDist)

		buff = geom.Buffer(buffDist) 
		vec_area = buff.GetArea()
		# print('Ratio Buffer to Raster: ', vec_area/ras_area)

		src_offset = bbox_to_pixel_offsets(rgt, buff.GetEnvelope(),rsize)

		if src_offset[2] <= 0 or src_offset[3] <=0:
			#if point falls outside raster grid, masked = None
			masked = None
		else:    
			if vec_area/ras_area >= fact:
				zooms = 1
			else:
				zooms = int(np.sqrt(ras_area/(vec_area/fact)))

			src_array = rb.ReadAsArray(*src_offset)

			# Calculate new geotransform of the feature subset
			if abs(theta) > 0: # if raster is rotated
				new_gt = ((rgt[0] + (src_offset[0] * res)), # top left x coordinate
						   res/ zooms,
						   0.,
						  (rgt[3] + (src_offset[1] * -res)), # top left y coordinate
						   0.,
						   -res/zooms)
			else:
				new_gt = ((rgt[0] + (src_offset[0] * rgt[1])),
						   rgt[1]/zooms,
						   0.0,
						  (rgt[3] + (src_offset[1] * rgt[5])),
						   0.0,
						   rgt[5]/zooms)
			
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

			# Mask the source data array with our current feature
			# we take the logical_not to flip 0<->1 to get the correct mask effect
			masked = np.ma.MaskedArray(src_re, mask = np.logical_not(rv_array))

		#If mask is empty, use NoData column
		if masked is None:
			pixel_count = {}
			pixel_count['NoData'] = 100
		else:
			keys, counts = np.unique(masked.compressed(), return_counts=True)

			# Create Dictionary of cmaps and counts
			pixel_count = dict(zip([cmap[k] for k in keys],
							  [np.asscalar(c) for c in counts]))
			pixtot = float(sum(pixel_count.values()))
			for k, v in pixel_count.items():
				pixel_count[k] = v / pixtot * 100

		# Create dictionary of station ids with pixel counts
		statDict[feat.GetField(fldname)] = pixel_count

	#clearing memory
		rvds = None
		mem_ds = None
 
	vds = None
	rds = None

	# Create dataframe from dictionary, transpose
	df = pd.DataFrame(statDict)
	df = df.T
	df = df.replace(np.nan,0) # NAN values are true zeros
	df = df.reset_index()
	cols = df.columns.tolist()
	cols[0] = fldname
	df.columns = cols

	## OUTPUT options
	if filenm == 'outputfile':
		cwd = os.getcwd()
		filenm = os.path.join(cwd, 'outputfile')
	if csvout == True:
		df.to_csv('{}.csv'.format(filenm), index=False)
	else:
		df.to_pickle('{}.pkl'.format(filenm))
