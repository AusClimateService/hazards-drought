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

Clone the repository and install in editable mode:

```bash
git clone <repo_url>
cd eddi
pip install -e .
```

## Usage
```import xarray as xr
from eddi.utils import meanVar, satVapourPressure, actVapourPressure
from eddi.eddi_utils import compute_EDDI
from eddi.utils import convert_sea_level_pressure_to_station_pressure, calculate_FAO56_pmpet

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

