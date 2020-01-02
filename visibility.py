"""
* Green Visibility Index Script
* Implements Huck and Gullick's Lightweight Viewshed algorithm
"""
from datetime import datetime
from fiona import open as fi_open
from subprocess import call, DEVNULL
from rasterio import open as rio_open
from numpy import zeros, unique, multiply

def coords2Array(input, x, y):
	"""
	* convert between coords and array position
	*  returns row,col (y,x) as expected by rasterio
	"""
	r, c = input.index(x, y)
	return int(r), int(c)


def array2Coords(input, row, col):
	"""
	* convert between array position and coords
	*  params are row,col (y,x) as expected by rasterio
	*  returns coords at the CENTRE of the cell
	"""
	x, y = input.xy(row, col)
	return int(x), int(y)


def distance(x1, y1, x2, y2):
	"""
	* Cartesian distance using Pythagoras
	"""
	return ((x1 - x2)**2 + (y1 - y2)**2)**0.5


def biliniearInterp(x,  y,  br,  bl,  tl,  tr,  xr,  xl,  yt,  yb):
	"""
	* Bilinear Interpolation
	*
	* x - the x coordinate of the current value
	* y - the y coordinate of the current value
	* br - value to the bottom right of the coordinate
	* bl - value to the bottom left of the coordinate
	* tl - value to the top left of the coordinate
	* tr - value to the top right of the coordinate
	* xr - the x coordinate for the values to the right of the coordinate
	* xl - the x coordinate for the values to the left of the coordinate
	* yt - the y coordinate for the values above the coordinate
	* yb - the y coordinate for the values below the coordinate
	"""
	# bottom two x values
	R1 = ((xr - x) / (xr - xl)) * bl + ((x - xl) / (xr - xl)) * br

	# top two x values
	R2 = ((xr - x) / (xr - xl)) * tl + ((x - xl) / (xr - xl)) * tr

	# combine with y values
	return ((yt - y) / (yt - yb)) * R1 + ((y - yb) / (yt - yb)) * R2


def getbiliniearInterpHeight(d, x, y, resolution):
	"""
	* Turn a set of coordinates into an interpolated height
	"""
	# work out the position in the cell (used to determine which cells to interpolate)
	xposincell = (x - d.bounds.left) % resolution
	yposincell = (y - d.bounds.bottom) % resolution

	# work out cells for the interpolation window
	(xl, xr) = doLeftRight(xposincell, x, y, resolution)
	(yb, yt) = doTopBottom(yposincell, x, y, resolution)

	return biliniearInterp(x, y, d[coords2Array(xr, yb)], d[coords2Array(xl, yb)],
		d[coords2Array(xl, yt)], d[coords2Array(xr, yt)], xr, xl, yt, yb)


def doLeftRight(xposincell, x, y, resolution):
	"""
	* Calculate x values for the left and right of the bilinear interp window
	"""
	# work out (left/ right)
	if xposincell > (resolution/2):        # right
		(xr, na) = array2Coords(coords2Array(x + resolution, y))
		(xl, na) = array2Coords(coords2Array(x, y))

	elif xposincell < (resolution/2):      # left
		(xl, na) = array2Coords(coords2Array(x - resolution, y))
		(xr, na) = array2Coords(coords2Array(x, y))

	else:
		(xl, na) = array2Coords(coords2Array(x, y))
		xl = xr

	return (xl, xr)


def doTopBottom(yposincell, x, y, resolution):
	"""
	* Calculate y values for the top and bottom of the bilinear interp window
	"""
	# work out (top/bottom)
	if yposincell > (resolution/2): # top
		(na, yt) = array2Coords(coords2Array(x, y + resolution))
		(na, yb) = array2Coords(coords2Array(x, y))

	elif yposincell < (resolution/2): # bottom
		(na, yb) = array2Coords(coords2Array(x, y - resolution))
		(na, yt) = array2Coords(coords2Array(x, y))

	else: # centre
		(na, yb) = array2Coords(coords2Array(x, y))
		yt = yb

	return (yb, yt)


def lineOfSight(x1, y1, x2, y2, observerHeight, targetHeight):
	"""
	 * Runs a single ray-trace from one point to another point,
	 *	set output data to 1 for each visible cell
	"""
	# init variables
	deltax = abs(x2 - x1)
	deltay = abs(y2 - y1)
	count = 0 				  # this is how many pixels we are in to our ray
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
		x1pixel, y1pixel = coords2Array(dsm_input, x, y)

		# distance travelled so far
		distanceTravelled = distance(x1, y1, x, y)

		# if we are on the first pixel (center of the circle)
		if (count == 0):

			# set the initial height
			initialHeight = dtm_data[x1pixel, y1pixel] + observerHeight
			#initialHeight = getbiliniearInterpHeight(dtm_data, x1pixel, y1pixel, resolution)

			# we of course can see ourselves
			output[x1pixel, y1pixel] = 1

		# we are on the second pixel
		elif (count == 1):

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

		# increment ourselves along the line
		count += 1

		# update iterators
		num += numadd		# Increase the numerator by the top of the fraction
		if (num >= den):	# Check if numerator >= denominator
			num -= den		# Calculate the new numerator value
			x += xinc1		# Change the x as appropriate
			y += yinc1		# Change the y as appropriate
		x += xinc2       	# Change the x as appropriate
		y += yinc2       	# Change the y as appropriate


def viewshed(x0, y0, radius, observerHeight, targetHeight):
	"""
	* Use Bresenham's Circle / Midpoint algorithm to determine endpoints for viewshed
	"""
	# initialise variables
	x = radius - resolution
	y = 0
	dx = resolution
	dy = resolution
	err = dx - (radius << 1)
	out = [];

	# loop around the 8 octant arcs
	while (x >= y):

		# calculate one ray in each octant
		lineOfSight(x0, y0, x0 + x, y0 + y, observerHeight, targetHeight)
		lineOfSight(x0, y0, x0 + y, y0 + x, observerHeight, targetHeight)
		lineOfSight(x0, y0, x0 - y, y0 + x, observerHeight, targetHeight)
		lineOfSight(x0, y0, x0 - x, y0 + y, observerHeight, targetHeight)
		lineOfSight(x0, y0, x0 - x, y0 - y, observerHeight, targetHeight)
		lineOfSight(x0, y0, x0 - y, y0 - x, observerHeight, targetHeight)
		lineOfSight(x0, y0, x0 + y, y0 - x, observerHeight, targetHeight)
		lineOfSight(x0, y0, x0 + x, y0 - y, observerHeight, targetHeight)

		# adjust for error
		if (err <= 0):
			y += resolution
			err += dy
			dy += (2 * resolution)
		else:
			x -= resolution
			dx += (2 * resolution)
			err += dx - (radius << 1)


''' --- SETTINGS --- '''
radius = 		100		# viewshed radius
checkrasters = 	True	# validate dimensions and resolution
# dtm_path = 		'/Users/jonnyhuck/Documents/_LABIB/DTMextendedFilledal1.tif'
# dsm_path = 		'/Users/jonnyhuck/Documents/_LABIB/DSMextendedFilledal1.tif'
# green_path = 	'/Users/jonnyhuck/Documents/_LABIB/GreenNoGreen_50km_5mFullal2.tif'
# boundary_path = '/Users/jonnyhuck/Documents/_LABIB/GM-Boundary/GMBoundary.shp'
dtm_path = 		'./data2/DTM_testArea.tif'
dsm_path = 		'./data2/DSM_testArea.tif'
green_path = 	'./data2/Green_noGreen_testArea.tif'
boundary_path = './data2/TestArea2.shp'
out_path = 		'./out/gvi-test-4.tif'
''' ---------------- '''

time = datetime.now()

# open bounbdary dataset and store bounds
with fi_open(boundary_path) as boundary:
	bounds = boundary.bounds

#  open the input dataset
with rio_open(dsm_path) as dsm_input:
	with rio_open(dtm_path) as dtm_input:

		# get the extents of the analysis in image space
		extent_b, extent_l = coords2Array(dtm_input, bounds[0], bounds[1])
		extent_t, extent_r = coords2Array(dtm_input, bounds[2], bounds[3])

		# get the data from band 1
		dsm_data = dsm_input.read(1)
		dtm_data = dtm_input.read(1)

		# get the resolution (assumes square pixels with int value)
		resolution = int(dtm_input.res[0])

		# create an array  at the same dimensions as data for output
		gvi = zeros((dtm_input.height, dtm_input.width))

		# get required info from files so we can close them
		width = dtm_input.width
		height = dtm_input.height
		crs = dtm_input.crs
		transform = dtm_input.transform

		# print(dtm_input.bounds)
		# print(bounds)
		# print(extent_l, extent_b, extent_r, extent_t)
		# print(dtm_input.height, dtm_input.width)
		# exit()

		#  open the greenspace
		with rio_open(green_path) as green_input:

			# verify that the input layers are all the same!
			if (checkrasters):
				if (dsm_input.width == dtm_input.width == green_input.width) == False or \
					(dsm_input.height == dtm_input.height == green_input.height) == False or \
					(dsm_input.res[0] == dtm_input.res[0] == green_input.res[0]) == False:
					print("rasters do not match!")
					print("width: \t\t", dsm_input.width == dtm_input.width == green_input.width)
					print("height: \t", dsm_input.height == dtm_input.height == green_input.height)
					print("resolution: \t", dsm_input.res[0] == dtm_input.res[0] == green_input.res[0])
					exit()

			# get the data from band 1
			green = green_input.read(1)

# loop through dataset rows
for r in range(height):

	# report progress
	print(f"{r / height * 100:.2f}% complete", flush=True)

	# enforce the buffer (so we don't go out of bounds)
	if r < extent_b and r > extent_t:

		# loop through dataset rows
		for c in range(width):

			# enforce the buffer (so we don't go out of bounds)
			if c > extent_l and c < extent_r:

				# create / reset output array at the same dimensions as data for output
				output = zeros((height, width))

				# call viewshed
				x, y = array2Coords(dtm_input, r, c)
				viewshed(x, y, radius, 2, 0)

				# get visible count (catch if nothing is visible)
				_, counts = unique(output, return_counts=True)
				visibleCellCount = counts[1] if counts.shape[0] == 2 else 0

				# multiply with the green dataset
				visibleGreen = multiply(output, green)

				# get value counts (catch if nothing is visible)
				_, counts = unique(visibleGreen, return_counts=True)
				visibleGreenCount = counts[1] if counts.shape[0] == 2 else 0

				# get the ratio for greenness in the view
				gvi[r,c] = visibleGreenCount / visibleCellCount

# output output dataset
with rio_open(out_path, 'w', driver='GTiff', height=height,
	width=width, count=1, dtype='float64', crs=crs,
	transform=transform
) as out:

	# write the output data to the tif
	out.write(gvi, 1)

# print how long it took
print(datetime.now() - time)

# use GDAL binary to calculate histogram and statistics
call(["gdalinfo", "-stats", out_path], stdout=DEVNULL)
