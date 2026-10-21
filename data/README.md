# Data

The lidar is not in this repository: it is several GB, and the 2014 survey carries no stated licence (see
`LICENSE-docs.md`). Both surveys are free to download from the Puget Sound LiDAR Consortium (PSLC), which needs a
free account for its data pages.

| Epoch | Survey | Acquisition | Delivered CRS | Terms |
|---|---|---|---|---|
| t0 = 2013 | PSLC **2013 Tulalip Partnership**, Watershed Sciences Inc. | 17 Apr - 31 Jul 2013 | NAD83(HARN) / Washington North, US survey feet; NAVD88 feet | Public domain; credit PSLC |
| t1 = 2014 | WSDOT **2014 Stillaguamish lidar, including Oso**, Quantum Spatial | 24 Mar - 9 Apr 2014 (4 flights) | NAD83(2011) / Washington North, US survey feet; NAVD88 (GEOID12A) feet | Freely available; no licence stated |

Both deliveries carry EPSG:32148 (metres) in their GeoKeys while the coordinates are in feet: step 0 of the workflow
overrides the label, so do not reproject the tiles yourself.

## Download

1. Create a PSLC account: https://pugetsoundlidar.ess.washington.edu/
2. 2013: LAS tiles of the 2013 Tulalip project,
   https://pugetsoundlidar.ess.washington.edu/lidardata/restricted/las/pslc2013/tulalip/delivery1and2
   (tile index: https://pugetsoundlidar.ess.washington.edu/lidardata/data_index/tiles/2013tulalip_tile_index.7z).
   Tiles used: 48121C7306, 48121C7307, 48121C7308, 48121C7309, 48121C7310, 48121C7311, 48121C7312, 48121C7313,
   48121C7314, 48121C7315, 48121C7316, 48121C7317, 48121C7318, 48121C7319, 48121C7320
3. 2014: LAS tiles of the Oso survey,
   https://pugetsoundlidar.ess.washington.edu/lidardata/restricted/las/oso_dot2014
   Tiles used: 9 tiles, 48121C7307 to 48121C7319 (the tiles covering the AOI above)
4. Put them in `data/raw/2013/` and `data/raw/2014/` (one merged file per epoch works too).

Area of interest, NAD83 / Washington North (US survey feet), as in the delivered files:
X 1,389,745 - 1,399,054, Y 462,123 - 475,846.

## Which flights are in the 2014 cloud?

The 2014 files store GPS week time. In the merged cloud it runs from 81,371 s to 410,925 s, that is from Sunday
15:36 to Thursday 11:09 local time (UTC-7). The only Sunday flight is 6 April and the only Thursday flight 3 April,
so the cloud is a mosaic of flights; `scripts/check_flight_dates.py` gives the share of points per day and hour:

| Day | Local hours | Share of points | Flight |
|---|---|---|---|
| Sunday | 15-18 h | 67.8 % | 6 April |
| Thursday | 10-11 h | 18.5 % | 3 April |
| Monday | 10 h, 12 h | 13.7 % | 7 April (assumed; 24 March cannot be excluded from week time) |

## Reference documents (not redistributed)

- Quantum Spatial (2014), technical data report of the 2014 survey: flights, datum, classes (16 = bridges).
  https://pugetsoundlidar.ess.washington.edu/lidardata/proj_reports/Oso_Landslide_Stillaguamish_River_LiDAR_Report.pdf
- 2013 Tulalip metadata: https://pugetsoundlidar.ess.washington.edu/lidardata/metadata/2013tulalip_deliv2_be.html
