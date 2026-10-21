# Inputs drawn for this study

Both GeoPackages are in NAD83 / Washington North (EPSG:32148, metres) and carry a `polygon_id` column.

- `stable_areas.gpkg`: 8 polygons of terrain assumed unchanged between 2013 and 2014, drawn by the author over the
  hillshades away from the slide and the river. They hold about 230,000 one-metre cells and carry the whole error
  model (slope term, interpolation term, systematic error).
- `zones_source_runout.gpkg`: `runout` is the USGS approximate extent of the 22 March 2014 deposit
  (`NewLandslidePoly`, Haugerud 2014, mapped from news photographs and described by USGS as not precise), converted
  from US survey feet. `source` is the convex hull of the USGS headscarp line and of the significant erosion deeper
  than 2 m within 400 m of it, minus the runout zone.
