"""
* Green Visibility Index Script
* Implements Huck and Gullick's Lightweight Viewshed algorithm
"""

from uuid import uuid1
from os import makedirs
from affine import Affine
from os.path import exists
from subprocess import call, DEVNULL
from rasterio import open as rio_open
from rasterio.transform import rowcol, xy
from numpy import zeros, unique, multiply

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


def distance(x1, y1, x2, y2):
	"""
	* Cartesian distance using Pythagoras
	"""
	return ((x1 - x2)**2 + (y1 - y2)**2)**0.5


# def biliniearInterp(x,  y,  br,  bl,  tl,  tr,  xr,  xl,  yt,  yb):
# 	"""
# 	* Bilinear Interpolation
# 	*
# 	* x - the x coordinate of the current value
# 	* y - the y coordinate of the current value
# 	* br - value to the bottom right of the coordinate
# 	* bl - value to the bottom left of the coordinate
# 	* tl - value to the top left of the coordinate
# 	* tr - value to the top right of the coordinate
# 	* xr - the x coordinate for the values to the right of the coordinate
# 	* xl - the x coordinate for the values to the left of the coordinate
# 	* yt - the y coordinate for the values above the coordinate
# 	* yb - the y coordinate for the values below the coordinate
# 	"""
# 	# bottom two x values
# 	R1 = ((xr - x) / (xr - xl)) * bl + ((x - xl) / (xr - xl)) * br
#
# 	# top two x values
# 	R2 = ((xr - x) / (xr - xl)) * tl + ((x - xl) / (xr - xl)) * tr
#
# 	# combine with y values
# 	return ((yt - y) / (yt - yb)) * R1 + ((y - yb) / (yt - yb)) * R2
#
#
# def getbiliniearInterpHeight(d, x, y, resolution):
# 	"""
# 	* Turn a set of coordinates into an interpolated height
# 	"""
# 	# work out the position in the cell (used to determine which cells to interpolate)
# 	xposincell = (x - d.bounds.left) % resolution
# 	yposincell = (y - d.bounds.bottom) % resolution
#
# 	# work out cells for the interpolation window
# 	(xl, xr) = doLeftRight(xposincell, x, y, resolution)
# 	(yb, yt) = doTopBottom(yposincell, x, y, resolution)
#
# 	return biliniearInterp(x, y, d[coords2Array(xr, yb)], d[coords2Array(xl, yb)],
# 		d[coords2Array(xl, yt)], d[coords2Array(xr, yt)], xr, xl, yt, yb)
#
#
# def doLeftRight(xposincell, x, y, resolution):
# 	"""
# 	* Calculate x values for the left and right of the bilinear interp window
# 	"""
# 	# work out (left/ right)
# 	if xposincell > (resolution/2):        # right
# 		(xr, na) = array2Coords(coords2Array(x + resolution, y))
# 		(xl, na) = array2Coords(coords2Array(x, y))
#
# 	elif xposincell < (resolution/2):      # left
# 		(xl, na) = array2Coords(coords2Array(x - resolution, y))
# 		(xr, na) = array2Coords(coords2Array(x, y))
#
# 	else:
# 		(xl, na) = array2Coords(coords2Array(x, y))
# 		xl = xr
#
# 	return (xl, xr)
#
#
# def doTopBottom(yposincell, x, y, resolution):
# 	"""
# 	* Calculate y values for the top and bottom of the bilinear interp window
# 	"""
# 	# work out (top/bottom)
# 	if yposincell > (resolution/2): # top
# 		(na, yt) = array2Coords(coords2Array(x, y + resolution))
# 		(na, yb) = array2Coords(coords2Array(x, y))
#
# 	elif yposincell < (resolution/2): # bottom
# 		(na, yb) = array2Coords(coords2Array(x, y - resolution))
# 		(na, yt) = array2Coords(coords2Array(x, y))
#
# 	else: # centre
# 		(na, yb) = array2Coords(coords2Array(x, y))
# 		yt = yb
#
# 	return (yb, yt)


def lineOfSight(x1, y1, x2, y2, resolution, observerHeight, targetHeight, dsm_data, dtm_data, a, output):
	"""
	 * Runs a single ray-trace from one point to another point, returning a list of visible cells
	"""

	# init variables
	deltax = abs(x2 - x1)	  # the difference in x axis
	deltay = abs(y2 - y1)	  # the difference in y axis
	biggestDYDXSoFar = 0 	  # biggest peak so far
	currentDYDX = 0 		  # current peak
	visible = 0   			  # used a char here, bool doesnt exist!
	tempHeight = 0 		      # temp height used for offset comparisons.
	distanceTravelled = 0     # how far we have travelled along the ray
	x = x1         		      # Start x off at the first pixel
	y = y1         		      # Start y off at the first pixel

	# get the direction on x axis
	if (x2 >= x1):      # The x-values are increasing
		xinc1 = resolution
		xinc2 = resolution
	else:             	# The x-values are decreasing
		xinc1 = -resolution
		xinc2 = -resolution

	# get the direction on y axis
	if (y2 >= y1):       # The y-values are increasing
		yinc1 = resolution
		yinc2 = resolution

	else:             	 # The y-values are decreasing
		yinc1 = -resolution
		yinc2 = -resolution

	# there is at least one x-value for every y-value
	if (deltax >= deltay):
		xinc1 = 0
		yinc2 = 0
		den = deltax
		num = deltax / 2
		numadd = deltay
		numpixels = deltax

	# there is at least one y-value for every x-value
	else:
		xinc2 = 0
		yinc1 = 0
		den = deltay
		num = deltay / 2
		numadd = deltax
		numpixels = deltay

	# loop along the line in increments of resolution
	for curpixel in range(0, numpixels, resolution):

		# pixel location for end point
		x1pixel, y1pixel = coords2Array(a, x, y)

		# distance travelled so far
		distanceTravelled = distance(x1, y1, x, y)

		# if we are on the first pixel (center of the circle)
		if (curpixel == 0):

			# set the initial height
			initialHeight = dtm_data[x1pixel, y1pixel] + observerHeight
			#initialHeight = getbiliniearInterpHeight(dtm_data, x1pixel, y1pixel, resolution)

			# we of course can see ourselves
			output[x1pixel, y1pixel] = 1

		# we are on the second pixel
		elif (curpixel == resolution):

			# first step definitely visible, just record the DY/DX and move on
			biggestDYDXSoFar = (dsm_data[x1pixel, y1pixel] - initialHeight) / distanceTravelled
			#biggestDYDXSoFar = getbiliniearInterpHeight(dsm_data, x1pixel, y1pixel, resolution)

			# again, obviously visible
			output[x1pixel, y1pixel] = 1

		# we are past the second pixel
		else:

			# the height of the top of the object in the landscape
			tempHeight = (dtm_data[x1pixel, y1pixel] - initialHeight + targetHeight) / distanceTravelled
			# tempHeight = (getbiliniearInterpHeight(dtm_data, x1pixel, y1pixel, resolution) - initialHeight + targetHeight) / distanceTravelled

			# the height of the base of the object in the landscape
			currentDYDX = (dsm_data[x1pixel, y1pixel] - initialHeight) / distanceTravelled
			# currentDYDX = (getbiliniearInterpHeight(dsm_data, x1pixel, y1pixel, resolution) - initialHeight) / distanceTravelled

			# is the angle bigger than we have seen?
			if ((tempHeight - biggestDYDXSoFar) >= 0):
				output[x1pixel, y1pixel] = 1

			# if this angle is greater than the biggest we have seen before, remember it.
			if (currentDYDX >= biggestDYDXSoFar):
				biggestDYDXSoFar = currentDYDX

		# update iterators
		num += numadd		# Increase the numerator by the top of the fraction
		if (num >= den):	# Check if numerator >= denominator
			num -= den		# Calculate the new numerator value
			x += xinc1		# Change the x as appropriate
			y += yinc1		# Change the y as appropriate
		x += xinc2       	# Change the x as appropriate
		y += yinc2       	# Change the y as appropriate

	# return list of visible cells
	return output


def viewshed(x0, y0, radius, resolution, observerHeight, targetHeight, dsm_data, dtm_data, a):
	"""
	* Use Bresenham's Circle / Midpoint algorithm to determine endpoints for viewshed
	"""

	# create output array at the same dimensions as data for viewshed
	output = zeros(dtm_data.shape)

	# initialise variables
	x = radius - resolution
	y = 0
	dx = resolution
	dy = resolution
	err = dx - (radius << 1)
	out = [];

	# loop around the 8 octant arcs
	while (x >= y):

		# calculate one ray in each octant, update output array each time
		output = lineOfSight(x0, y0, x0 + x, y0 + y, resolution, observerHeight, targetHeight, dsm_data, dtm_data, a, output)
		output = lineOfSight(x0, y0, x0 + y, y0 + x, resolution, observerHeight, targetHeight, dsm_data, dtm_data, a, output)
		output = lineOfSight(x0, y0, x0 - y, y0 + x, resolution, observerHeight, targetHeight, dsm_data, dtm_data, a, output)
		output = lineOfSight(x0, y0, x0 - x, y0 + y, resolution, observerHeight, targetHeight, dsm_data, dtm_data, a, output)
		output = lineOfSight(x0, y0, x0 - x, y0 - y, resolution, observerHeight, targetHeight, dsm_data, dtm_data, a, output)
		output = lineOfSight(x0, y0, x0 - y, y0 - x, resolution, observerHeight, targetHeight, dsm_data, dtm_data, a, output)
		output = lineOfSight(x0, y0, x0 + y, y0 - x, resolution, observerHeight, targetHeight, dsm_data, dtm_data, a, output)
		output = lineOfSight(x0, y0, x0 + x, y0 - y, resolution, observerHeight, targetHeight, dsm_data, dtm_data, a, output)

		# adjust for error
		if (err <= 0):
			y += resolution
			err += dy
			dy += (2 * resolution)
		else:
			x -= resolution
			dx += (2 * resolution)
			err += dx - (radius << 1)

	# return the resulting viewshed
	return output


def f(mask):
	"""
	* main function for running with parallel.py
	"""

	# create an array  at the same dimensions as data for output
	gvi = zeros((mask["meta"]["height"], mask["meta"]["width"]))

	# get pixel references for aoi extents
	min_r, min_c  = coords2Array(mask["meta"]["transform"], mask["aoi"].bounds[0], mask["aoi"].bounds[3])
	max_r, max_c  = coords2Array(mask["meta"]["transform"], mask["aoi"].bounds[2], mask["aoi"].bounds[1])

	# loop through dataset rows and columns
	for r in range(min_r, max_r):
		for c in range(min_c, max_c):

			# call viewshed
			x, y = array2Coords(mask["meta"]["transform"], r, c)
			output = viewshed(x, y,
				mask["options"]["radius"], 			# radius
				int(mask['meta']['transform'][0]),	# resolution
				mask["options"]["o_height"], 		# observor height
				mask["options"]["t_height"],		# target height
				mask["dsm"], 						# dsm dataset
				mask["dtm"],						# dtm dataset
				mask["meta"]["transform"])			# affine transform

			# get visible count (catch if nothing is visible)
			_, counts = unique(output, return_counts=True)
			visibleCellCount = counts[1] if counts.shape[0] == 2 else 0

			# multiply with the green dataset
			visibleGreen = multiply(output, mask["green"])

			# get value counts (catch if nothing is visible)
			_, counts = unique(visibleGreen, return_counts=True)
			visibleGreenCount = counts[1] if counts.shape[0] == 2 else 0

			# get the ratio for greenness in the view
			gvi[r,c] = visibleGreenCount / visibleCellCount

	# clip gvi to aoi bounds
	gvi = gvi[min_r:max_r, min_c:max_c]

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
* Dummy data for testing
"""
if __name__ == '__main__':
	from rasterio.crs import CRS
	from numpy.random import random
	from shapely.geometry import Polygon
	b = (387589.3499999987, 394741.04330000095, 388410.3499999987, 395401.04330000095)
	dsm = random((173, 206))
	dtm = random((173, 206))
	green = random((173, 206))
	mask = {
		'dtm': dtm,
		'dsm': dsm,
		'green': green,
		'meta': {
			'driver': 'GTiff',
			'dtype': 'float32',
			'nodata': None,
			'width': 206,
			'height': 173,
			'count': 1,
			'crs': CRS.from_dict(init='epsg:27700'),
			'transform': Affine(5.0, 0.0, 387484.9969,0.0, -5.0, 395504.9997)
			},
		'aoi': Polygon([ [b[0], b[1]], [b[0], b[3]], [b[2], b[3]], [b[2], b[1]] ]),
		'options': {
			'radius': 100,
			'o_height': 2,
			't_height': 0
			}
		}
	f(mask)
