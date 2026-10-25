# Inputs drawn for this study

The GeoPackages are in NAD83 / Washington North (EPSG:32148, metres) and carry a `polygon_id` column.

- `stable_areas.gpkg`: 8 polygons of terrain assumed unchanged between 2013 and 2014, drawn by the author over the
  hillshades away from the slide and the river. They cover about 0.66 km² (66,500 cells of 3 m), redrawn on 2026-10-25, and carry the
  whole error model (slope term, interpolation term, systematic error).
- `zones_source_runout.gpkg`: `runout` is the USGS approximate extent of the 22 March 2014 deposit
  (`NewLandslidePoly`, Haugerud 2014, mapped from news photographs and described by USGS as not precise), converted
  from US survey feet. `source` is the convex hull of the USGS headscarp line and of the significant erosion deeper
  than 2 m within 400 m of it, minus the runout zone.
- `water_2013.gpkg`: the North Fork Stillaguamish channel in 2013 (one polygon, about 0.24 km²), drawn by the author
  over the 2013 DTM hillshade, because the 2013 delivery has no water class.
- `water_2014.gpkg`: the water in 2014 (133 polygons): the delivered class 9 mask, polygonised and edited by hand
  over the hillshade to add the impoundment lake, where the water returns are too sparse to be bridged.
