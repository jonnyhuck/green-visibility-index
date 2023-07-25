# Greenness Visibility Index <img src="logo.png" align="right" height="139"/>

This software uses a Python Viewshed implementation to calculate a GVI surface from a DSM, DTM and Binary Greenness Surface. It is part of the paper published at:

[Labib, S.M., Huck, J.J. and Lindley, S., 2021. Modelling and mapping eye-level greenness visibility exposure using multi-source data at high spatial resolutions. **Science of the Total Environment**, 755, p.143050](https://doi.org/10.1016/j.scitotenv.2020.143050). 

## Usage

This software is designed to run in a parallel computing environment, but can also be run on a single process. It is called using a simple command line interface (CLI). Note that you should call `parallel.py`, even if you only want to run the script on a single core. Here is an example call:

```bash
python parallel.py --dtm ./dtm.tif  --dsm ./dsm.tif --green ./greenspace.tif --aoi ./aoi.shp  --padding 800 --sections 1 1 --out ./mark_test.tif
```

Here are the arguments and explanations:

```
options:
  -h, --help            show this help message and exit
  --dtm DTM             File path for the DTM/DEM layer
  --dsm DSM             File path for the DSM layer for analysis
  --green GREEN         File path for the Binary (1/0) GreenSpace layer
  --aoi AOI             Boundary of Area of Interest
  --padding PADDING     Viewshed radius in CRS units (normally metres)
  --sections SECTION 		x and y divisions for running in parallel (e.g. 16 divisions would be `4 4`)
  --out OUT             Filepath for result file (.tif)
```

The inputs for the model comprise:

* *Raster datasets:* Digital Elevation/Terrain Model (DEM/DTM), Digital Surface Model (DSM) and a binary Greenspace layer (1=green, 0=not green)
  * All three of these layers should be 'aligned': i.e., identical CRS, bounds and resolution. The 'align rasters' tool in QGIS ius a good way to achieve this
* *Vector datasets:* An Area of Interest (AOI) polygon. This should be smaller than the raster datasets by at least the viewshed radius, to allow the cells near the edge to calculate a full viewshed
  * An easy way to calculate this is to use the 'Extract Layer Extent' tool in QGIS to get a polygon representing your raster data extent, then buffer it by the **negative** of the viewshed distance (e.g., -1000m for a 1000m viewshed).
* The viewshed radius (`padding`) in the CRS units (normally metres)
* Instructions for the area should be subdivided for parallel computing (`sections`). Each section is computed separately and then the results are stitched back together. For example, `2 4`would subdivide the are into 8 separate sections (2 wide by 4 tall), each of which would be processed using a separate process.
  * If you don't want to run the tool in parallel, then simply put `1 1`

## Other languages

An R port of this repository is available at: [![STBrinkmann R Library](https://badgen.net/badge/STBrinkmann/R%20Library/blue?icon=github)](https://github.com/STBrinkmann/GVI)
