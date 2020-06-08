"""
* Green Visibility Index Script
"""
from uuid import uuid1
from os import makedirs
from affine import Affine
from os.path import exists
from math import exp, hypot
# from time import perf_counter
from rasterio import open as rio_open
from rasterio.transform import rowcol, xy
from skimage.draw import line, disk, circle_perimeter
from numpy import zeros, unique, multiply, array, column_stack


def coords2Array(a, x, y):
	"""
	* convert between coords and array position
	*  returns row,col (y,x) as expected by rasterio
	"""
	r, c = rowcol(a, x, y)
	return int(r), int(c)


def array2Coords(a, row, col):
	"""
	* convert between array position and coords
	*  params are row,col (y,x) as expected by rasterio
	*  returns coords at the CENTRE of the cell
	"""
	x, y = xy(a, row, col)
	return int(x), int(y)


def lineOfSight(r0, c0, r1, c1, observer_height, resolution, target_height, dsm_data, dtm_data, output):
	"""
	 * Runs a single ray-trace from one point to another point, returning a list of visible cells
	"""

	# init variables for loop
	cur_dydx = 0 		  	# current dydx (base of object)
	max_dydx = 0 	  		# biggest dydx so far
	# top_dydx = 0 		    # current dydx (top of object)
	distance_travelled = 0  # how far we have travelled along the ray

	# get the viewer height
	height0 = dtm_data[(r0, c0)] + observer_height

	# get the pixels in the line (excluding the first one	)
	pixels = column_stack(line(r0, c0, r1, c1))[1:]

	# loop along the pixels in the line
	for r, c in pixels:

		# distance travelled so far
		distance_travelled = hypot(c0 - c, r0 - r)

		''' comment this out as long as we use 0 as target offset '''
		## set cell as visible if the height of the top of the object from the DTM > previous max
		# top_dydx = (dsm_data[(r, c)] - height0 + target_height) / distance_travelled
		# if (top_dydx >= max_dydx):
		# 	output[(r, c)] = 1
		#
		## update max dydx the height of the base of the object on the DSM > previous max
		# cur_dydx = (dsm_data[(r, c)] - height0) / distance_travelled
		# if (cur_dydx > max_dydx):
		# 	max_dydx = cur_dydx

		# update max dydx the height of the base of the object on the DSM > previous max
		cur_dydx = (dsm_data[(r, c)] - height0) / (distance_travelled * resolution)
		if (cur_dydx > max_dydx):
			max_dydx = cur_dydx
			output[(r, c)] = 1

	# return updated output surface
	return output


def viewshed(r0, c0, radius_px, resolution, observerHeight, targetHeight, dsm_data, dtm_data, a):
	"""
	* Use Bresenham's Circle / Midpoint algorithm to determine endpoints for viewshed
	"""

	# create output array at the same dimensions as data for viewshed
	output = zeros(dtm_data.shape)

	# set the start location as visible automatically
	output[(r0, c0)] = 1

	# get pixels in the circle
	for r, c in column_stack(circle_perimeter(r0, c0, radius_px)):

		# calculate line of sight to each pixel
		output = lineOfSight(r0, c0, r, c, resolution, observerHeight, targetHeight, dsm_data, dtm_data, output)

	# return the resulting viewshed
	return output


def f(mask):
	"""
	* main function for running with parallel.py
	"""

	# create an output array at the same dimensions as data for output
	gvi = zeros((mask["meta"]["height"], mask["meta"]["width"]))

	# radius in pixels
	radius_px = int(mask["options"]["radius"] // mask['meta']['transform'][0])

	# build weighting mask
	weighting_mask = zeros((radius_px*2, radius_px*2))
	for r, c in column_stack(disk((radius_px, radius_px), radius_px, shape=weighting_mask.shape)):
		weighting_mask[(r, c)] = exp(-0.0003 * (hypot(radius_px - c, radius_px - r) * mask['meta']['transform'][0]))

	# get pixel references for aoi extents
	min_r, min_c  = coords2Array(mask["meta"]["transform"], mask["aoi"].bounds[0], mask["aoi"].bounds[3])
	max_r, max_c  = coords2Array(mask["meta"]["transform"], mask["aoi"].bounds[2], mask["aoi"].bounds[1])

	# loop through dataset rows and columns
	for r in range(min_r, max_r+1):
		for c in range(min_c, max_c+1):

			# print(r, c)
			# t1_start = perf_counter()

			# call (weighted) viewshed
			output = viewshed(r, c, radius_px, 		# coords and radius in pixels
				mask['meta']['transform'][0],		# resolution of datasets
				mask["options"]["o_height"], 		# observer height
				mask["options"]["t_height"],		# target height
				mask["dsm"], 						# dsm dataset
				mask["dtm"],						# dtm dataset
				mask["meta"]["transform"])			# affine transform

			# extract the viewshed data from the output surface and apply weighting mask
			visible = output[r-radius_px:r+radius_px, c-radius_px:c+radius_px] * weighting_mask

			# print(f"\tviewshed took {perf_counter() - t1_start}s", visible.shape, visible.sum())
			# t1_start = perf_counter()

			# multiply extract of (weighted) viewshed with extract of (weighted) green dataset
			visible_green = visible * (mask["green"][r-radius_px:r+radius_px, c-radius_px:c+radius_px] * weighting_mask)

			# print(f"\tvisible green {perf_counter() - t1_start}s", visible_green.shape, visible_green.sum())
			# t1_start = perf_counter()

			# get the ratio for greenness in the view
			gvi[(r,c)] = visible_green.sum() / visible.sum()

			# print(f"\tgvi took {perf_counter() - t1_start}s", gvi[(r,c)])
			# print()

	# clip gvi to aoi bounds
	gvi = gvi[min_r:max_r+1, min_c:max_c+1]

	# check that tmp folder exists
	if not exists('./tmp/'):
	    makedirs('tmp')

	# make unique filename
	filepath = f'./tmp/{str(uuid1())[:8]}.tif'

	# output file with updated dimensions and transform
	with rio_open(filepath, 'w',
		driver =	mask["meta"]["driver"],
		height = 	gvi.shape[0],
		width =		gvi.shape[1],
		count =		mask["meta"]["count"],
		dtype =		'float64',
		crs =		mask["meta"]["crs"],
		transform =	Affine(
			mask['meta']['transform'][0],
			mask['meta']['transform'][1],
			mask["aoi"].bounds[0],
			mask['meta']['transform'][3],
			mask['meta']['transform'][4],
			mask["aoi"].bounds[3]),
		) as dst:
		dst.write(gvi, 1)

	# return the filepath to the result
	return filepath


"""
* Do not call this script directly
"""
if __name__ == '__main__':
	print("please call this script using parallel.py")
