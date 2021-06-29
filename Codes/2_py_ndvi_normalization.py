#!/usr/bin/python

import glob

import os
import osr
import sys
import gdal
import time
import numpy as np

def create_output_file(base_filepath, out_filepath, dataType = gdal.GDT_Float32, imageFormat = 'GTiff'):
    
  driver = gdal.GetDriverByName(imageFormat)
  base_ds = gdal.Open(base_filepath)

  x_start, pixel_width, _, y_start, _, pixel_height = base_ds.GetGeoTransform()
  x_size = base_ds.RasterXSize 
  y_size = base_ds.RasterYSize
  
  out_srs = osr.SpatialReference()
  out_srs.ImportFromWkt(base_ds.GetProjectionRef())

  output_img_ds = driver.Create(out_filepath, x_size, y_size, 1, dataType, ['COMPRESS=LZW', 'TILED=True', 'BIGTIFF=YES'])
  output_img_ds.SetGeoTransform((x_start, pixel_width, 0, y_start, 0, pixel_height))
  output_img_ds.SetProjection(out_srs.ExportToWkt())

  return output_img_ds

def merge_unique_values(result, uniq_vals, count_vals):

	for i in range(0, len(uniq_vals)):
		uniq_val = uniq_vals[i]
		if uniq_val not in result:
			result[uniq_val] = 0
		result[uniq_val] = result[uniq_val] + count_vals[i]

def calcNdviMin(uniqValues, nPixels1pct):
	valueCount = 0
	values = []
	
	for i in range(0,len(uniqValues)):
		uniqValue = uniqValues[i]
		valueCount = valueCount + ndviFreq[uniqValue]
		if valueCount > nPixels1pct:
			break
		else:
			values.append(uniqValue * ndviFreq[uniqValue])

	ndvi_min = np.sum(values) / valueCount
	return ndvi_min
	print(ndvi_max)

def calcNdviMax(uniqValues, nPixels1pct):
	
	valueCount = 0
	values = []

	for i in range(len(uniqValues)-1,-1,-1):
		uniqValue = uniqValues[i]
		valueCount = valueCount + ndviFreq[uniqValue]
		if valueCount > nPixels1pct:
			break
		else:
			values.append(uniqValue * ndviFreq[uniqValue])

	ndvi_max = np.sum(values) / valueCount
	return ndvi_max
	print(ndvi_max)


def calcHistogram(ndviDs, ndviFreq, totalPixels, XChunk):
	Xsize = ndviDs.RasterXSize
	YSize = ndviDs.RasterYSize

	ndviBand = ndviDs.GetRasterBand(1)

	for x in range(0,Xsize,XChunk):
		if (x+XChunk) > Xsize:
			XChunk = Xsize - x

		ndvi = ndviBand.ReadAsArray(x,0,XChunk,YSize);
		ndvi[ndvi==nodataValue] = np.nan
		
		#pasture = ndimage.binary_dilation(pasture, structure=struct2).astype(pasture.dtype)
		# Binary dilation https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.ndimage.morphology.binary_dilation.html#scipy.ndimage.morphology.binary_dilation

		maskFilter = np.logical_not(np.isnan(ndvi))

		ndvi_vals = ndvi[maskFilter]
		
		if (len(ndvi_vals) > 0):
			
			uniq_vals, count_vals = np.unique(ndvi_vals, return_counts=True)
			totalPixels = totalPixels + np.sum(count_vals)

			merge_unique_values(ndviFreq, uniq_vals, count_vals)
			print("Reading " + str(float(x)/float(Xsize)*100.0) + "%")

	return totalPixels

	return totalPixels

start = time.time()

input_dir = sys.argv[1]
output_dir = sys.argv[2]
nodataValue = int(sys.argv[3])

os.makedirs(output_dir,  exist_ok = True)

driver = gdal.GetDriverByName('GTiff')
input_files = glob.glob(input_dir + "/*.tif")
print(input_files) ###

totalSum = 0
XChunk = 1024

ndvi001 = []
ndvi099 = []
ndviFreq = {}
totalPixels = 0

for ndviFilepath in input_files:
	ndviDs = gdal.Open(ndviFilepath, gdal.GA_ReadOnly)
	totalPixels = calcHistogram(ndviDs, ndviFreq, totalPixels, XChunk)

nPixels1pct = totalPixels*0.01
uniqValues = sorted(ndviFreq.keys())

minNdviVal = np.min(uniqValues)
maxNdviVal = np.max(uniqValues)

ndviMin = float(calcNdviMin(uniqValues, nPixels1pct))
ndviMax = float(calcNdviMax(uniqValues, nPixels1pct))

print('ndvi001:', ndviMin, minNdviVal, nPixels1pct)
print('ndvi099:', ndviMax, maxNdviVal, nPixels1pct)

#ndviMin = 0.37099849587271955 # -0.5696202516555786 3554811.74
#ndviMax = 0.7466436368063499 # 0.9097744226455688 3554811.74

for ndviFilepath in input_files:
	ndviDs = gdal.Open(ndviFilepath, gdal.GA_ReadOnly)
	outputFile = os.path.join(output_dir, os.path.basename(ndviFilepath))

	ndviBand = ndviDs.GetRasterBand(1)
	
	Xsize = ndviDs.RasterXSize
	YSize = ndviDs.RasterYSize

	resultDs = create_output_file(ndviFilepath, outputFile)
	resultBand = resultDs.GetRasterBand(1)

	for x in range(0,Xsize,XChunk):
		if (x+XChunk) > Xsize:
			XChunk = Xsize - x

		ndvi = ndviBand.ReadAsArray(x,0,XChunk,YSize);

		ndvi[ndvi==nodataValue] = np.nan

		maskFiler = np.logical_not(np.isnan(ndvi))

		if (len(ndvi[maskFiler]) > 0):
			print('NDVI:', np.mean(ndvi[maskFiler]), np.min(ndvi[maskFiler]), np.max(ndvi[maskFiler]))

			ndvi[np.logical_and(maskFiler, ndvi < ndviMin)] = ndviMin
			ndvi[np.logical_and(maskFiler, ndvi > ndviMax)] = ndviMax
			ndvi[maskFiler] = (ndvi[maskFiler] - ndviMin) / (ndviMax - ndviMin)

			print('NDVIq:', np.mean(ndvi[maskFiler]), np.min(ndvi[maskFiler]), np.max(ndvi[maskFiler]))
			
		resultBand.WriteArray(ndvi, x, 0)

		print("Writing " + str(float(x)/float(Xsize)*100.0) + "%")

	end = time.time()
	print('Time processing', end - start)
	resultBand.FlushCache()

	ndviDs = None
	pastureDs = None
	resultDs = None

	end = time.time()
	print('All time', end - start)

