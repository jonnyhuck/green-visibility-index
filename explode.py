from geopandas import read_file, GeoSeries

data = read_file("../gm_boundary/gm_64.shp")

for i, feature in data.iterrows():
    GeoSeries(feature.geometry, crs=data.crs).to_file(f"../gm_boundary/exploded/{i+1}.shp")

print("done")
