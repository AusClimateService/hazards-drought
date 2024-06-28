# ACS Hazard Team on Drought and Changes in Aridity

## Description
GitHub repository for ACS Hazard Team on Drought and Changes in Aridity to store, track and develop code.
Last updated: 28 June 2024

## Indices considered by the hazard team:
**Currently available:**
### Standardised Precipitation Index (SPI)
The Standardized Precipitation Index (SPI) is a widely used index that measures the amount of precipitation over a specific period relative to the long-term average for that period. It standardizes precipitation as a z-score, making it possible to compare different locations and time scales. The SPI is typically used to identify and quantify the severity of droughts. Positive SPI values indicate wetter-than-average conditions, while negative values indicate drier-than-average conditions. We use a SPI-3, aggregated over three months and a value of -1 or lower as a drought metric, signifying moderate drought conditions. It is valuable in water resource management, agriculture, and climate studies for its simplicity and effectiveness in drought monitoring (McKee et al, 1993).

### Aridity Index (AI)
The Aridity Index (AI) is a numerical indicator used to quantify the dryness of a region. It is calculated as the ratio of annual precipitation to potential evapotranspiration. Lower values of AI indicate more arid conditions, while higher values suggest more humid conditions. The AI is commonly used in climatology, agriculture, and environmental studies to classify climates, assess water availability, and manage land and water resources (UNEP, 1992).
Aridity categories based on AI values are as follows:
- Hyper-Arid: AI < 0.05
- Arid: 0.05 ≤ AI < 0.2
- Semi-Arid: 0.2 ≤ AI < 0.5
- Dry Sub-Humid: 0.5 ≤ AI < 0.65
- Humid: AI ≥ 0.65

### Rainfall percentiles
Rainfall percentiles are statistical measures used to evaluate and interpret precipitation data. They indicate the relative ranking of a given rainfall amount within a historical reference period. The 15th percentile we use on three-month rainfall aggregation (smiliar to SPI) represents a value below which 15% of the observed rainfall amounts fall, meaning it is drier than 85% of the reference data (WMO, 2017). The 15th percentile approximates an SPI value of -1, but instead if an inde value it provides the expected rainfall amount. This means that if the precipitation for a given time frame (e.g., GWL) is lower than the amount that is exceeded 85% of the time in the historical record, it is considered indicative of drought. This index helps in understanding and quantifying the severity and frequency of drought events by focusing on the lower tail of the precipitation distribution.


**For future deliveries:**
- Standardised soil mositure index (SSMI)
- Standardised Precipitation and Evapotranspiration Index (SPEI)
- Standardised runoff index (SRI)
- Evaporative Demand Drought Index (EDDI)
- Evaporative Stress Index (ESI)

## Products:
Status of the NCRA deliverables. 

The three dots (in order from first/top/left to last/bottom/right) represent the datasets used to compute indices:
- Dot 1: Pre-processed BARPA/CCAM – downscaled but NOT bias-corrected, 5 km (deliverable for 30 June)
- Dot 2: Bias-corrected BARPA/CCAM – downscaled AND bias-corrected, 5 km (deliverable for 31 July)
- Dot 3: National Hydrological Projections (NHP1.0) based on CMIP5 – bias-corrected, 5km
Where only one dot is in the cell the format type does not apply to the metric, e.g. no time series for rainfall 15th percentile.
 
In terms of the colors:
- :green_circle: The data is available in its final official form
- :yellow_circle: The data creation is currently in progress and available soon
- :red_circle: The data processing has not yet started
- :white_circle: Not intended for delivery/not applicable

| Index/metric | time series (ts) | GWLs ts | GWLs 2D | MME 2D | MME 2D change | Scheduled<br>delivery date | Data<br>location | Last update
|-----         | :-:              |:-:      |:-:      |:-:     |:-:            |------------    |-----             |-----
| SPI3 |:green_circle:<br>:yellow_circle:<br>:white_circle:|:green_circle:<br>:yellow_circle:<br>:white_circle:|:green_circle:<br>:yellow_circle:<br>:white_circle:|:green_circle:<br>:yellow_circle:<br>:white_circle:|:green_circle:<br>:yellow_circle:<br>:white_circle:|<ul><li>30 June</li><li>31 July</li><li>N/A</li></ul>|/g/data/ia39/ncra/<br>drought_aridity/spi/|27/06/24
| Rainfall 15th prctl |:white_circle:|:white_circle:|:green_circle:<br>:green_circle:<br>:white_circle:|:green_circle:<br>:green_circle:<br>:white_circle:|:green_circle:<br>:green_circle:<br>:white_circle:|<ul><li>30 June</li><li>31 July</li><li>N/A</li></ul>|/g/data/ia39/ncra/<br>drought_aridity/<br>rainfall_percentiles/|27/06/24
| AI |:white_circle:<br>:white_circle:<br>:green_circle:|:white_circle:<br>:white_circle:<br>:green_circle:|:white_circle:<br>:white_circle:<br>:green_circle:|:white_circle:<br>:white_circle:<br>:green_circle:|:white_circle:<br>:white_circle:<br>:green_circle:|<ul><li>N/A</li><li>N/A</li><li>31 July</li></ul>|/g/data/ia39/ncra/<br>drought_aridity/ai/|28/06/24

**Figures:** Figures for each 2D metric (SPI <= -1, AI, rainfall percentiles) and GWL as well as for changes relative to GWL 1.2 are located in the index directories (see data location in table above) in the sub directory `/figures/`. Example figure for 'time spent in drought (SPI <= -1)' change for GWL 2.0 (rel. to GWL 1.2) for the 10th, 50th and 90th percentile of the multi-model ensemble:
![Time spent in drought (SPI3 <= -1) change GWL 2 relative to GWL 1.2](figures/SPI3_change_GWL20_rel_to_GWL12.png)

**Area statistics:** Figures for changes of each 2D metric (SPI <= -1, AI, rainfall percentiles) and GWL relative to GWL 1.2 are located in the index directories (see data location in table above) in the sub directory `/MME_change_spatial_summary_statistics/`.
Example summary statistic for NCRA regions:
![GWL2.0_to_GWL1.2_change_for_MME_50th_percentile](figures/area_stat_GWL2.0_to_GWL1.2_change_for_MME_50th_percentile.png)

## References
* Mckee, T.B., Doesken, N.J., Kleist, J. (1993). The relationship of drought frequency and duration to time scales. In: Proceedings of the 8th Conference on Applied Climatology. American Meteorological Society, Boston, MA, pp. 179–183.
* United Nations Environment Programme (UNEP) (1992). World Atlas of Desertification. Edward Arnold. https://wedocs.unep.org/20.500.11822/42137
* WMO. (2017). WMO Guidelines on the Calculation of Climate Normals. WMO-No. 1203, 1203, 29. https://library.wmo.int/doc_num.php?explnum_id=4166


## Contributing
Open to contributions. 

## Authors and acknowledgment
Hazard team:
- [ ] David Hoffmann (BOM, lead)
- [ ] Tess Parker (CSIRO, alternate lead)
- [ ] Jessica Bhardwaj (BOM, contributor)
