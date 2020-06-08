from glob import glob
from rasterio import open as rio_open
from rasterio.merge import merge as rio_merge

# get all tif images in directory
images = glob('../CEF-Result/out/*.tif')

# get files
files = []
for file in images:
    files.append(rio_open(file, 'r'))

# merge result files
merged, out_transform = rio_merge(files, method='max')

# output dataset to raster (hardcoded crs as was causing error)
with rio_open("../CEF-Result/800m.tif", 'w', driver='GTiff', height=merged.shape[1],
    width=merged.shape[2], count=1, dtype='float64', crs="EPSG:27700",
    transform=out_transform) as out:

    out.write(merged[0], 1)

# close all of the files
for file in files:
    file.close()

print("done")
