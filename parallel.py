"""
* Divide and Conquour script to divide, process in parallel and merge again without having
*	to worry about edge effects.
* This version is intended to be used with call gvi.py
"""

from gvi import f
from math import ceil
from os import makedirs
from os.path import exists
from datetime import datetime
from multiprocessing import Pool
from fiona import open as fi_open
from subprocess import call, DEVNULL
from rasterio import open as rio_open
from rasterio.mask import mask as rio_mask
from rasterio.merge import merge as rio_merge
from shapely.geometry import Polygon, mapping


def extractPolygon(src, poly):
	"""
	* extract a subset from a raster according to the specified bounds
	* returns the dataset (numpy array) and metadata object
	"""

	# extract the bit of src that intersects poly
	out_data, out_transform = rio_mask(src, poly, crop=True, all_touched=True)

	# create the metadata for the dataset
	out_meta = src.meta
	out_meta.update({
		"height": out_data.shape[1],
		"width": out_data.shape[2],
		"transform": out_transform
	})

	# return the mask and the metadata
	return out_data[0], out_meta


def outputFiles(crs, masks):
	"""
	* Output illustrative GIS (Shapefile and GeiTiff) files of the extracted
	*  zones and the aoi polygons
	"""

	# create the file in write mode, with the required CRS (Proj4 string) and an empty schema
	with fi_open('./tmp/aoi.shp', 'w',
		driver = 'ESRI Shapefile',
		crs = crs,
		schema = {'geometry': 'Polygon', 'properties': {}}
		) as o:

		# loop through resulting masks
		for m in range(len(masks)):

			# write the splitter to shapefile
			o.write({'geometry': mapping(masks[m]["aoi"]),'properties': {}})

			# create a raster file
			with rio_open(f'./tmp/{m}.tif', 'w',
				driver =	masks[m]["meta"]["driver"],
				height = 	masks[m]["meta"]["height"],
				width =		masks[m]["meta"]["width"],
				count =		masks[m]["meta"]["count"],
				dtype =		masks[m]["meta"]["dtype"],
				crs =		masks[m]["meta"]["crs"],
				transform =	masks[m]["meta"]["transform"],
				) as dst:

				# write the data to the raster
				dst.write(masks[m]["dtm"], 1)


# def f(mask):
# 	"""
# 	* just for testing...
# 	"""
# 	return mask["meta"]['crs']


'''
SETTINGS
'''
''' NB: We assume that the three datasets are identical in terms of bounds and resolution '''
out_path = './out/gvi.tif'
dtm_path = './data2/DTM_testArea.tif'
dsm_path = './data2/DSM_testArea.tif'
green_path = './data2/Green_noGreen_testArea.tif'
boundary_path = './data2/TestArea2.shp'
padding = 100		# the required padding for each mask in m (e.g. radius of a viewshed)
sections = [2, 2]	# x, y
'''---'''

# log start time
time = datetime.now()

print(f"a {sections[0]}x{sections[1]} grid, {sections[0]*sections[1]} processes, started at {time}.")

# initialise masks array for the results
masks = []

# get aoi bounds from shapefile
with fi_open(boundary_path) as boundary:
	bounds = boundary.bounds

# read in raster data
with rio_open(dtm_path) as dtm_input:
	with rio_open(dsm_path) as dsm_input:
		with rio_open(green_path) as green_input:

			# verify raster dimensions and resolution
			if (dsm_input.width == dtm_input.width == green_input.width) == False or \
				(dsm_input.height == dtm_input.height == green_input.height) == False or \
				(dsm_input.res[0] == dtm_input.res[0] == green_input.res[0]) == False:
				print("rasters do not match!")
				print("width: \t\t", dsm_input.width == dtm_input.width == green_input.width)
				print("height: \t", dsm_input.height == dtm_input.height == green_input.height)
				print("resolution: \t", dsm_input.res[0] == dtm_input.res[0] == green_input.res[0])
				exit()

			# store metadata from dtm
			meta = dtm_input.meta

			# read data bands
			dtm_data = dtm_input.read(1)
			dsm_data = dsm_input.read(1)
			green_data = green_input.read(1)

			# get clip dimensions
			aoiWidth =  ceil(((bounds[2]) - (bounds[0])) / sections[0])
			aoiHeight = ceil(((bounds[3]) - (bounds[1])) / sections[1])
			print(f"aoi dimensions (m): {aoiWidth}, {aoiHeight}")

			# pre-calculate origin location polygon extraction
			xOrigin = bounds[0] - padding
			yOrigin = bounds[1] - padding

			# loop through mask matrix
			for x in range(sections[0]):
				for y in range(sections[1]):

					# construct a Shapely polygon for use in processing
					polygon = Polygon([
						(bounds[0] + (aoiWidth * x) - padding, bounds[1] + (aoiHeight * y) 	- padding), # bl
						(bounds[0] + (aoiWidth * x) - padding, bounds[1] + (aoiHeight * (y+1)) + padding), # tl
						(bounds[0] + (aoiWidth * (x+1)) + padding, bounds[1] + (aoiHeight * (y+1)) + padding), # tr
						(bounds[0] + (aoiWidth * (x+1)) + padding, bounds[1] + (aoiHeight * y) - padding)  # br
					])

					# construct a shapely polygon for use in analysis and trimming the results
					aoi = Polygon([
						(bounds[0] + (aoiWidth * x), bounds[1] + (aoiHeight * y)), # bl
						(bounds[0] + (aoiWidth * x), bounds[1] + (aoiHeight * (y+1))), # tl
						(bounds[0] + (aoiWidth * (x+1)), bounds[1] + (aoiHeight * (y+1))), # tr
						(bounds[0] + (aoiWidth * (x+1)), bounds[1] + (aoiHeight * y)) # br
					])

					# extract the polygon and append to masks list
					dtm, meta = extractPolygon(dtm_input, [polygon])
					dsm, _ = extractPolygon(dsm_input, [polygon])
					green, _ = extractPolygon(green_input, [polygon])

					# make result object and append to masks list
					masks.append({
						"dtm": dtm,
						"dsm": dsm,
						"green": green,
						"meta": meta,
						"aoi": aoi,
						"options": {
							"radius":	100,	# viewshed radius
							"o_height": 2,		# observer height
							"t_height": 0		# target height
						}
					})

			# output files for debugging purposes
			# outputFiles(dtm_input.crs, masks)


# make as many processes as are required and launch them
with Pool(sections[0] * sections[1]) as p:
	results = p.map(f, masks)

print(results)

# open all result files in read mode
files = []
for filepath in results:
	files.append(rio_open(filepath, 'r'))

# merge result files
merged, out_transform = rio_merge(files)

# update the metadata
out_meta = files[0].meta.copy()
out_meta.update({
		"height": merged.shape[1],
		"width": merged.shape[2],
		"transform": out_transform,
		"dtype": 'float64'
		})

# check that tmp folder exists
if not exists('./out/'):
	makedirs('out')

# create a raster file
with rio_open(out_path, 'w', **out_meta) as dst:
	dst.write(merged[0], 1)

# use GDAL binary to calculate histogram and statistics
call(["gdalinfo", "-stats", out_path], stdout=DEVNULL)

# print how long it took
print(datetime.now() - time)
print("done!")
