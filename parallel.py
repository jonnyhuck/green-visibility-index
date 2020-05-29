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
from argparse import ArgumentParser
from subprocess import call, DEVNULL
from rasterio import open as rio_open
from rasterio.transform import rowcol, xy
from rasterio.mask import mask as rio_mask
from rasterio.merge import merge as rio_merge
from shapely.geometry import Polygon, mapping


def coords2Array(a, coords):
	"""
	* convert between coords and array position
	*  returns row,col (y,x) as expected by rasterio
	"""
	x, y = coords
	r, c = rowcol(a, x, y)
	return int(r), int(c)


def array2Coords(a, index):
	"""
	* convert between array position and coords
	*  params are row,col (y,x) as expected by rasterio
	*  returns coords at the CENTRE of the cell
	"""
	row, col = index
	x, y = xy(a, row, col)
	return int(x), int(y)


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
	* Output illustrative GIS (Shapefile and GeoTiff) files of the extracted
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


# get settings from args
parser = ArgumentParser(description="Labib's Greenspace Visibility Tool")
parser.add_argument('--dtm', help='DTM layer for analysis', required=True)
parser.add_argument('--dsm', help='DSM layer for analysis', required=True)
parser.add_argument('--green', help='Binary GreenSpace layer for analysis', required=True)
parser.add_argument('--aoi', help='Boundary of Area of Interest', required=True)
parser.add_argument('--padding', help='Required padding for each mask in CRS units (e.g. radius of a viewshed)', required=True)
parser.add_argument('--sections', nargs=2, help='x and y divisions for masking (e.g. 16 divisions would be `4 4`)', required=True)
parser.add_argument('--out', help='Path for result file', required=True)
args = vars(parser.parse_args())

# get args
dtm_path = args['dtm']
dsm_path = args['dsm']
green_path = args['green']
boundary_path = args['aoi']
padding = int(args['padding'])
sections = [ int(x) for x in args['sections'] ]
out_path = args['out']

# log start time and log start
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

			# adjust bounds to the raster by transforming to image space and back again
			minx, miny = array2Coords(dtm_input.transform, coords2Array(dtm_input.transform, [bounds[0], bounds[1]]))
			maxx, maxy = array2Coords(dtm_input.transform, coords2Array(dtm_input.transform, [bounds[2] + dtm_input.res[0], bounds[3] + dtm_input.res[0]]))

			# get width and height of study area
			w, h = (maxx - minx), (maxy - miny)

			# get width and height of all aois (absorb offcut into first entry)
			aoiWidth =  [int(w / sections[0])] * sections[0]
			aoiWidth[0] += w - sum(aoiWidth)
			aoiHeight = [int(h / sections[1])] * sections[1]
			aoiHeight[0] += h - sum(aoiHeight)

			# pre-calculate origin location polygon extraction
			xOrigin = bounds[0] - padding
			yOrigin = bounds[1] - padding

			# loop through mask matrix
			for x in range(sections[0]):
				for y in range(sections[1]):

					# construct a Shapely polygon for use in processing
					polygon = Polygon([
						(bounds[0] + sum(aoiWidth[:x]) - padding, bounds[1] + sum(aoiHeight[:y]) - padding), # bl
						(bounds[0] + sum(aoiWidth[:x]) - padding, bounds[1] + sum(aoiHeight[:y+1]) + padding), # tl
						(bounds[0] + sum(aoiWidth[:x+2]) + padding, bounds[1] + sum(aoiHeight[:y+1]) + padding), # tr
						(bounds[0] + sum(aoiWidth[:x+1]) + padding, bounds[1] + sum(aoiHeight[:y]) - padding)  # br
					])

					# construct a shapely polygon for use in analysis and trimming the results
					aoi = Polygon([
						(bounds[0] + sum(aoiWidth[:x]), bounds[1] + sum(aoiHeight[:y])), # bl
						(bounds[0] + sum(aoiWidth[:x]), bounds[1] + sum(aoiHeight[:y+1])), # tl
						(bounds[0] + sum(aoiWidth[:x+1]), bounds[1] + sum(aoiHeight[:y+1])), # tr
						(bounds[0] + sum(aoiWidth[:x+1]), bounds[1] + sum(aoiHeight[:y]))  # br
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
							"radius": padding,	# viewshed radius
							"o_height": 1.7,		# observer height
							"t_height": 0		# target height
						}
					})

					# print(masks[0])
					# exit()

			# output files for debugging purposes
			# outputFiles(dtm_input.crs, masks)
			# exit()

# make as many processes as are required and launch them
with Pool(sections[0] * sections[1]) as p:
	results = p.map(f, masks)

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

# create a raster file
with rio_open(out_path, 'w', **out_meta) as dst:
	dst.write(merged[0], 1)

# use GDAL binary to calculate histogram and statistics
# call(["gdalinfo", "-stats", out_path], stdout=DEVNULL)

# close all of the files
for file in files:
	file.close()

# print how long it took
print(datetime.now() - time)
print("done!")
