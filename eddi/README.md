# ACS_eddi_et0

Python package to calculate:

- FAO-56 Penman–Monteith reference evapotranspiration (ET₀)
- Evaporative Drought Demand Index (EDDI)

Input requirements:
- Temperature: max/min in [K] or [˚C]
- Humidity: specific humidity in [kg/kg], relative humidity [%]
- Radiation: [W/m-2]
- Sea level pressure: [Pa]
- Wind: 10m in [m/s]
- Elevation: [m]

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

