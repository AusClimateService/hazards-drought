# ACS Hazard Team on Drought and Changes in Aridity

## Description
GitHub repository for ACS Hazard Team on Drought and Changes in Aridity to store, track and develop code. 

## Indices considered by the hazard team:
- Standardised Precipitation Index (SPI)
- Standardised soil mositure index (SSMI)
- Standardised runoff index (SRI)
- Rainfall percentiles
- Aridity Index (AI)
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

| Index/metric | time series (ts) | GWLs ts | GWLs 2D | MME 2D | MME 2D change | Notes | Data<br>location | Last update
|-----         | :-:              |:-:      |:-:      |:-:     |:-:            |-----    |-----          |-----
| SPI |:yellow_circle:<br>:red_circle:<br>:white_circle:|:yellow_circle:<br>:red_circle:<br>:white_circle:|:yellow_circle:<br>:red_circle:<br>:white_circle:|:yellow_circle:<br>:red_circle:<br>:white_circle:|:yellow_circle:<br>:red_circle:<br>:white_circle:|<ul><li>deliverable for 30 June</li><li>deliverable for 31 July</li><li>*N/A</li></ul>|/g/data/ia39/ncra/<br>drought_aridity/spi/|20/06/24
| SPI 'time spent in drought' |:yellow_circle:<br>:red_circle:<br>:white_circle:|:yellow_circle:<br>:red_circle:<br>:white_circle:|:yellow_circle:<br>:red_circle:<br>:white_circle:|:yellow_circle:<br>:red_circle:<br>:white_circle:|:yellow_circle:<br>:red_circle:<br>:white_circle:|<ul><li>deliverable for 30 June</li><li>deliverable for 31 July</li><li>N/A</li></ul>|/g/data/ia39/ncra/<br>drought_aridity/spi/|20/06/24
| Rainfall 15th prctl |:white_circle:|:white_circle:|:yellow_circle:<br>:red_circle:<br>:white_circle:|:yellow_circle:<br>:red_circle:<br>:white_circle:|:yellow_circle:<br>:red_circle:<br>:white_circle:|<ul><li>deliverable for 30 June</li><li>deliverable for 31 July</li><li>N/A</li></ul>|/g/data/ia39/ncra/<br>drought_aridity/<br>rainfall_percentiles/|20/06/24
| AI |:white_circle:<br>:white_circle:<br>:green_circle:|:white_circle:<br>:white_circle:<br>:green_circle:|:white_circle:<br>:white_circle:<br>:green_circle:|:white_circle:<br>:white_circle:<br>:green_circle:|:white_circle:<br>:white_circle:<br>:green_circle:|<ul><li>N/A</li><li>N/A</li><li>deliverable for 31 July</li></ul>|/g/data/ia39/ncra/<br>drought_aridity/ai/|20/06/24



## Roadmap
Coming soon..

## Contributing
Open to contributions. 

## Authors and acknowledgment
Hazard team:
- [ ] David Hoffmann (BOM, lead)
- [ ] Tess Parker (CSIRO, alternate lead)
- [ ] Jessica Bhardwaj (BOM, contributor)
