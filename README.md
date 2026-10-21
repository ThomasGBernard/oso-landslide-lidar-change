# Oso (SR530) landslide: lidar change detection and sediment budget

[![tests](https://github.com/GITHUB_USER/oso-landslide-lidar-change/actions/workflows/ci.yml/badge.svg)](https://github.com/GITHUB_USER/oso-landslide-lidar-change/actions/workflows/ci.yml)
[![Code: MIT](https://img.shields.io/badge/code-MIT-blue.svg)](LICENSE)
[![Text and figures: CC BY 4.0](https://img.shields.io/badge/text%20%26%20figures-CC%20BY%204.0-lightgrey.svg)](LICENSE-docs.md)

How much material left the source of the 2014 Oso landslide, and how much stayed on the valley floor? This repository
answers with repeat airborne lidar (2013 against spring 2014), a DEM of difference thresholded by a spatially
varying Level of Detection, and a budget whose uncertainty is stated term by term. The result matches the published
estimates within a few percent.

![Significant elevation change 2013–2014](figures/dod_map.png)

## The event

On 22 March 2014 a slope of glacial deposits about 180–190 m high failed above the North Fork Stillaguamish River,
near Oso, Washington. The debris crossed the valley, buried the Steelhead Haven neighbourhood and State Route 530,
dammed the river, and killed 43 people (Iverson et al., 2015). It is one of the deadliest landslides in U.S. history,
and it has been studied in detail by the U.S. Geological Survey (Haugerud, 2014; Iverson et al., 2015; Anderson et
al., 2017).

## Data

| | Pre-event (t0) | Post-event (t1) |
|---|---|---|
| Survey | Puget Sound LiDAR Consortium, 2013 Tulalip Partnership | Washington State DOT, 2014 Stillaguamish lidar including Oso |
| Vendor | Watershed Sciences Inc. | Quantum Spatial |
| Dates | 17 April – 31 July 2013 | flights 24 March – 9 April 2014; this area at least 3 and 6 April |
| Ground points used | 10.8 M (1.2 pts/m² under the forest) | 15.7 M |

Neither survey is redistributed here; `data/README.md` says where to download them (free account) and under which
terms.

## The trap in the data: metres that are feet

Both deliveries say EPSG:32148 (NAD83 / Washington North, **metres**) in their GeoKeys. The coordinates are in **US
survey feet**: the vendor report says so, and the Z range of 235–1,070 only fits the valley floor at Oso (72 m) in
feet. Read as labelled, every vertical change would be 3.28 times too large and every volume 35.3 times too large,
and a CRS check on the EPSG code would pass. The workflow overrides the label with the true feet-based CRS and
converts to metres on the same datum before anything else (`src/oso_cd/prepare.py`).

The 2013 survey is on NAD83(HARN) and the 2014 survey on NAD83(2011). The difference is centimetric at Oso and is
absorbed by the registration on stable terrain (next section).

## Method

1. **Units and ground.** PDAL pipeline: CRS override (ftUS), reprojection to metres, ASPRS class 2 only.
2. **Registration.** py4dgeo stable-area ICP (Yang and Schwieger, 2022; Besl and McKay, 1992; Anders et al., 2026)
   on the ground points, thinned to one point per 1 m voxel, in a local frame. Residual offset found: a few
   centimetres (−0.04, −0.05, −0.09 m).
3. **DTMs.** 1 m inverse-distance gridding of the ground points within res/√2 of each cell centre (as PDAL
   `writers.gdal`), then a GDAL inverse-distance fill of the gaps inside the common footprint. Every cell is
   recorded as measured or interpolated, with its distance to the nearest measured cell. Under the forest, 36 % of
   the 2013 cells have no ground point.
4. **Error model** from 8 hand-drawn stable areas (about 380,000 cells), where the true change is zero (after
   Hugonnet et al., 2022, with robust NMAD statistics after Höhle and Höhle, 2009):
   - slope term σ_DoD(θ) = √(a² + (b·tanθ)²) fitted on measured cells: a = 0.061 m, b = 0.279 m;
   - interpolation term as a function of the distance to the nearest measured cell;
   - systematic error σ_sys = 0.028 m, the spread of the per-area medians;
   - cells steeper than 45° get no σ and never count as change.
5. **DoD and Level of Detection.** ΔZ = Z₂₀₁₄ − Z₂₀₁₃ is significant where |ΔZ| > 1.96·√(σ₀² + σ₁²), the 95 %
   threshold of Brasington et al. (2003), Lane et al. (2003) and Wheaton et al. (2010). Median LoD: 0.33 m.
6. **Budget by zone.** Source zone (behind the USGS headscarp) and runout zone (USGS approximate deposit extent,
   Haugerud, 2014), each split into measured and interpolated cells.

![Error model](figures/error_model.png)

On the stable areas the DoD has an NMAD of 0.094 m on measured cells and a median offset of +0.009 m, in line with
the 0.05–0.13 m that Anderson et al. (2017) report for the same surveys.

## Results

| Zone | Area (km²) | Erosion (Mm³) | Deposition (Mm³) | Net (Mm³) | Stated uncertainty (Mm³) |
|---|---|---|---|---|---|
| Source | 0.175 | 3.92 | 0.04 | **−3.89** | ± 0.08 |
| Runout | 0.900 | 0.35 | 3.98 | **+3.63** | ± 0.27 |

- **Mass balance.** Over source and runout together, erosion is 4.27 Mm³ and deposition 4.02 Mm³ (net −6 % of the
  erosion): the budget closes within bulking and uncertainty.
- **Uncertainty.** Stated as the systematic term (σ_sys × zone area: ± 0.005 and ± 0.025 Mm³) plus a conservative
  interpolation term (a 1 m bias of one sign on every interpolated cell of change: ± 0.07 and ± 0.25 Mm³). Survey
  noise is not what limits this budget: interpolation under the forest and the zone outlines are. 65 % of the source
  volume and 26 % of the runout volume come from cells interpolated in at least one epoch, at a median of 1 m from a
  ground point. The approximate USGS outline is not included in the figure.

![Budget against published values](figures/budget_vs_published.png)

## Comparison with published results

| Quantity | This study | Published |
|---|---|---|
| Volume that left the source | 3.89 Mm³ | ≈ 3.8 Mm³: 8.3 Mm³ mobilised minus ≈ 4.5 Mm³ stranded in the source (range 2.8–4.7), Iverson et al. (2015) |
| Valley-floor deposit | 3.63 Mm³ | ≈ 3.9 Mm³: 3.4 Mm³ still in place in July 2015 plus 0.506 Mm³ eroded since 24 March 2014, Anderson et al. (2017) |
| Survey noise on stable ground | NMAD 0.094 m, offset +0.009 m | SD 0.05–0.13 m, systematic −0.009 to +0.044 m, Anderson et al. (2017) |

Two points of interpretation:

- **A DoD measures what left, not what moved.** Iverson et al. (2015) estimate 8.3 Mm³ mobilised, but about 4.5 Mm³
  of it slid and stayed in the source area, where it barely changes the surface elevation. The mobilised volume needs
  a basal slip surface; differencing gives the volume that left.
- **The 2014 surface is about two weeks old.** The post-event cloud comes from flights of 3 and 6 April 2014, not
  from the 24 March emergency flight, so a little erosion of the deposit by the river had already happened.
  The runout gain is still within 8 % of the deposit volume of Anderson et al. (2017).

A point-based cross-check with M3C2 (Lague et al., 2013), run in the author's in-house toolbox, gives 3.14 Mm³ of
erosion over the whole site: a lower bound, because M3C2 finds no distance where the 2013 survey has no ground points
under the forest, whereas the DoD fills those cells.

## Limitations

- The zone outlines are approximate (the USGS runout extent was mapped from news photographs) and their effect on the
  volumes is not quantified.
- The interpolated share is large under the forest; the conservative term above bounds it but is not a
  statistical confidence interval.
- The 2014 cloud mixes at least two flights; any change between them inside the slide is part of the result.
- The impoundment lake upstream of the deposit is a water surface in 2014 (flat at 87.6 m): it lies outside both
  zones and is not part of the budget.

## Reproduce

```bash
conda env create -f environment.yml
conda activate oso-cd
# download the lidar into data/raw/2013 and data/raw/2014 (data/README.md)
python -m oso_cd.pipeline config/oso.toml          # or run notebooks/oso_sediment_budget.ipynb
pytest
```

The run writes `results/` (volumes, error model, validation against the reference run) and `figures/`. Large
intermediate products go to the git-ignored `data/work/`. A full run takes about 7 minutes on a laptop.

| Path | Content |
|---|---|
| `src/oso_cd/` | the workflow: `prepare`, `registration`, `dtm`, `uncertainty`, `dod`, `plots`, `pipeline` |
| `config/oso.toml` | every parameter of the run |
| `inputs/` | stable areas and zones (drawn for this study) |
| `notebooks/` | the executed notebook, step by step |
| `results/`, `figures/` | published results; `results/reference/` holds the toolbox reference run |
| `scripts/check_flight_dates.py` | GPS time to flight days |

This repository is a lean, readable rewrite of the Package 1 workflow of the author's in-house toolbox (v0.4.1). It
reproduces the toolbox's zonal volumes within 2 % (`results/validation.json`).

## References

Site and published results

- Anderson, S.W., Keith, M.K., Magirl, C.S., Wallick, J.R., Mastin, M.C., and Foreman, J.R., 2017, Geomorphic
  response of the North Fork Stillaguamish River to the State Route 530 landslide near Oso, Washington: U.S.
  Geological Survey Scientific Investigations Report 2017–5055, 85 p. https://doi.org/10.3133/sir20175055
- Haugerud, R.A., 2014, Preliminary interpretation of pre-2014 landslide deposits in the vicinity of Oso,
  Washington: U.S. Geological Survey Open-File Report 2014–1065, 4 p. https://doi.org/10.3133/ofr20141065
- Iverson, R.M., George, D.L., Allstadt, K., Reid, M.E., Collins, B.D., Vallance, J.W., Schilling, S.P., Godt, J.W.,
  Cannon, C.M., Magirl, C.S., Baum, R.L., Coe, J.A., Schulz, W.H., and Bower, J.B., 2015, Landslide mobility and
  hazards: implications of the 2014 Oso disaster: Earth and Planetary Science Letters, 412, 197–208.
  https://doi.org/10.1016/j.epsl.2014.12.020

Data

- Puget Sound LiDAR Consortium, 2014, 2013 Tulalip Partnership topographic lidar (Watershed Sciences Inc.), public
  domain. https://www.fisheries.noaa.gov/inport/item/50174
- Washington State Department of Transportation, 2014, 2014 Stillaguamish lidar including Oso (Quantum Spatial),
  distributed by PSLC. https://pugetsoundlidar.ess.washington.edu/lidardata/restricted/nonpslc/oso2014
- Quantum Spatial, 2014, Oso Landslide and Stillaguamish River LiDAR technical data report.
  https://pugetsoundlidar.ess.washington.edu/lidardata/proj_reports/Oso_Landslide_Stillaguamish_River_LiDAR_Report.pdf

Methods and software

- Anders, K., et al., 2026, py4dgeo: open-source scientific software for topographic change analysis in 3D/4D
  geographic point clouds: SoftwareX, 34, 102670. https://doi.org/10.1016/j.softx.2026.102670
- Besl, P.J., and McKay, N.D., 1992, A method for registration of 3-D shapes: IEEE Transactions on Pattern Analysis
  and Machine Intelligence, 14(2), 239–256. https://doi.org/10.1109/34.121791
- Brasington, J., Langham, J., and Rumsby, B., 2003, Methodological sensitivity of morphometric estimates of coarse
  fluvial sediment transport: Geomorphology, 53(3–4), 299–316. https://doi.org/10.1016/S0169-555X(02)00320-3
- Butler, H., et al., 2024, PDAL: Point Data Abstraction Library. https://doi.org/10.5281/zenodo.2616780
- Höhle, J., and Höhle, M., 2009, Accuracy assessment of digital elevation models by means of robust statistical
  methods: ISPRS Journal of Photogrammetry and Remote Sensing, 64(4), 398–406.
  https://doi.org/10.1016/j.isprsjprs.2009.02.003
- Hugonnet, R., et al., 2022, Uncertainty analysis of digital elevation models by spatial inference from stable
  terrain: IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 15, 6456–6472.
  https://doi.org/10.1109/JSTARS.2022.3188922
- Lague, D., Brodu, N., and Leroux, J., 2013, Accurate 3D comparison of complex topography with terrestrial laser
  scanner: application to the Rangitikei canyon (N-Z): ISPRS Journal of Photogrammetry and Remote Sensing, 82,
  10–26. https://doi.org/10.1016/j.isprsjprs.2013.04.009
- Lane, S.N., Westaway, R.M., and Hicks, D.M., 2003, Estimation of erosion and deposition volumes in a large,
  gravel-bed, braided river using synoptic remote sensing: Earth Surface Processes and Landforms, 28(3), 249–271.
  https://doi.org/10.1002/esp.483
- Wheaton, J.M., Brasington, J., Darby, S.E., and Sear, D.A., 2010, Accounting for uncertainty in DEMs from repeat
  topographic surveys: improved sediment budgets: Earth Surface Processes and Landforms, 35(2), 136–156.
  https://doi.org/10.1002/esp.1886
- Yang, Y., and Schwieger, V., 2022, Supervoxel-based targetless registration and identification of stable areas for
  deformed point clouds: Journal of Applied Geodesy. https://doi.org/10.1515/jag-2022-0031

BibTeX of all of them: `references.bib`. How to cite this work: `CITATION.cff`.

## Author

Thomas Bernard, PhD in geomorphology and remote sensing (lidar), river morphodynamics and sediment transport,
Trento, Italy. WEBSITE_URL
