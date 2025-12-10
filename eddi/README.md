# ACS_eddi_et0

Python package to calculate:

- FAO-56 Penman–Monteith reference evapotranspiration (ET₀), following methodology in Allen et al. (2005). FAO Irrigation and Drainage Paper 56: Crop Evapotranspiration – Guidelines for Computing Crop Water Requirements. [https://www.mesonet.org/images/site/ASCE_Evapotranspiration_Formula.pdf](https://www.mesonet.org/images/site/ASCE_Evapotranspiration_Formula.pdf)
- Evaporative Drought Demand Index (EDDI), following methodology in Hobbins, M. T., Wood, A., McEvoy, D. J., Huntington, J. L., Morton, C., Anderson, M., & Hain, C. (2016). The Evaporative Demand Drought Index. Part I: Linking Drought Evolution to Variations in Evaporative Demand. Journal of Hydrometeorology, 17(6), 1745-1761. [https://doi.org/10.1175/JHM-D-15-0121.1](https://doi.org/10.1175/JHM-D-15-0121.1)

Input requirements for ET₀:
- Temperature: max/min in [K] or [˚C]
- Humidity: specific humidity in [kg/kg], relative humidity [%]
- Radiation: [W/m-2]
- Sea level pressure: [Pa]
- Wind: 10m in [m/s]
- Elevation: [m]

Input requirements for EDDI:
- ET₀: [mm/day]

## Installation

Clone the repository:
```bash
git clone <repo>
```
In your python script/Jupyter notebook:
```python
import sys
sys.path.append("/path/to/eddi-repo")
```

## Usage
```import xarray as xr
from eddi.pet_eddi import compute_daily_FAO56_PET, compute_EDDI

# Load data
tasmax = xr.open_dataarray("tasmax.nc")
tasmin = xr.open_dataarray("tasmin.nc")
rsds = xr.open_dataarray("rsds.nc")
hurs = xr.open_dataarray("hurs.nc")
sfcWind = xr.open_dataarray("sfcWind.nc")
psl = xr.open_dataarray("psl.nc")
elev = xr.open_dataarray("elev.nc")

# Compute daily ET0
ET0 = compute_daily_FAO56_PET(tasmax, tasmin, rsds, hurs, sfcWind, psl, elev)

# Compute EDDI
EDDI = compute_EDDI(ET0, ndays=30)
```

