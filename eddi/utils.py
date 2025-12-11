"""
pet_utils.py
============
Utility functions for computing FAO-56 Penman–Monteith reference ET₀ and
supporting meteorological transforms.

Author
------
David Hoffmann (david.hoffmann@bom.gov.au)

Description
-----------
This module contains helper functions used for the calculation of reference
evapotranspiration (FAO-56) and for downstream drought metrics such as EDDI.
The functions implement equations from:

Allen et al. (2005). FAO Irrigation and Drainage Paper 56:
Crop Evapotranspiration – Guidelines for Computing Crop Water Requirements.
https://www.mesonet.org/images/site/ASCE_Evapotranspiration_Formula.pdf
"""

import numpy as np
import xarray as xr
import git
import sys

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def get_git_hash():
    """Returns the git hash for the working repository"""
    git_repo = git.Repo(sys.argv[0], search_parent_directories=True)
    git_root = git_repo.git.rev_parse("--show-toplevel")
    git_hash = git.Repo(git_root).heads[0].commit
    git_text = " (Git hash: %s)" %(str(git_hash)[0:7])
    return git_text
    
def convertTemperature(da):
    """
    Convert Kelvin to Celsius if necessary.
    """
    if da.attrs.get("units", "").lower() in ["k", "kelvin"]:
        print("Converting temperature from K to C.")
        out = da - 273.15
        out.attrs["units"] = "C"
        return out
    else:
        return da

def meanVar(varmax, varmin):
    """
    Daily mean from min/max following FAO-56.
    """
    return (varmax + varmin) / 2.0


# ---------------------------------------------------------------------------
# Vapour pressure calculations
# ---------------------------------------------------------------------------

def satVapourPressure(tasmax=None, tasmin=None, tas=None):
    """
    Parameters
    ----------
    tasmax/tasmin/tas : xarray.DataArray or float [degC]

    Returns
    -------
    esat : same type as input
        Saturation Vapour Pressure in kilo Pascals (kPa)

    Notes
    -----
    Saturation vapour pressure [esat [kPa] = f(Tmx [degC], Tmn [degC]) using  FAO-56 Chapter 3 - Meteorological Data for formulae, Allen 2005, Equations 6 and 7
    Accepts either tasmax/tasmin or tas.
    """
    if tasmax is not None and tasmin is not None:
        e0_tmax = 0.6108 * np.exp(17.27 * tasmax / (tasmax + 237.3))
        e0_tmin = 0.6108 * np.exp(17.27 * tasmin / (tasmin + 237.3))
        return (e0_tmax + e0_tmin) / 2.0

    elif tas is not None:
        return 0.6108 * np.exp(17.27 * tas / (tas + 237.3))

    else:
        raise ValueError("Provide either (tasmax & tasmin) or tas.")


def actVapourPressure(esat, psl=None, huss=None, hurs=None,
                      hursmax=None, hursmin=None, tasmax=None, tasmin=None):
    """
    Parameters
    ----------
    esat/psl/huss/hurs/hursmax/hursmin/tasmax/tasmin : xarray.DataArray or float
 
    Returns
    -------
    ea : same type as input
        actual vapour pressure in kilo Pascals (kPa)

    Notes
    -----
    Formula from Allen et al., FAO-56 (2005), Eq. 3:

        P = 101.3 * ((293 - 0.0065 * z)/293)^5.26 [kPa]

    Input is in Pa, output is in Pa.
    Actual vapour pressure (ea) using one of the FAO-56 compliant pathways:

    - If hurs is provided: ea = RH * esat
    - If hursmax/hursmin provided: FAO-56 Eq. 17
    - If specific humidity + pressure provided: WCI eq.
    """
    if hurs is not None:
        hurs_aligned = hurs.reindex_like(esat)
        ea = (hurs_aligned / 100.0) * esat

    elif (hursmax is not None) and (hursmin is not None):
        hursmax_al = hursmax.reindex_like(tasmin)
        hursmin_al = hursmin.reindex_like(tasmax)
        e0_tmin = satVapourPressure(tas=tasmin)
        e0_tmax = satVapourPressure(tas=tasmax)
        ea = 0.5 * (e0_tmin * (hursmax_al / 100.0) +
                    e0_tmax * (hursmin_al / 100.0))

    elif (psl is not None) and (huss is not None):
        w = huss / (1.0 - huss)
        P_kPa = psl / 1000.0
        ea = (P_kPa * w) / (0.622 + w)

    else:
        raise ValueError("Provide hurs OR (hursmax+hursmin) OR (psl+huss)")

    # enforce ea ≤ esat for physical consistency
    return ea.clip(max=esat)


# ---------------------------------------------------------------------------
# Pressure and wind conversions
# ---------------------------------------------------------------------------

def convert_sea_level_pressure_to_station_pressure(mslp_pa, h_m):
    """
    Convert mean sea-level pressure (Pa) to station pressure (Pa) using FAO-56 (Allen 2005, Eq. 3).

    Parameters
    ----------
    mslp_pa : xarray.DataArray or float
        Mean sea-level pressure in Pascals (Pa)
    h_m : xarray.DataArray, float, or int
        Station elevation in meters

    Returns
    -------
    P_station : same type as input
        Atmospheric pressure at station elevation in Pascals (Pa)

    Notes
    -----
    Formula from Allen et al., FAO-56 (2005), Eq. 3:

        P = 101.3 * ((293 - 0.0065 * z)/293)^5.26 [kPa]

    Input is in Pa, output is in Pa.
    """
    # Convert Pa → kPa
    mslp_kpa = mslp_pa / 1000.0

    # FAO-56 equation - Adjust proportionally to actual sea-level pressure
    P_station_kpa = mslp_kpa * ((293.0 - 0.0065 * h_m)/293.0)**5.26

    # Return in Pa
    P_station_pa = P_station_kpa * 1000.0

    # Preserve attrs if input is xarray.DataArray
    if hasattr(mslp_pa, "attrs"):
        P_station_pa = P_station_pa.copy()
        P_station_pa.attrs = mslp_pa.attrs.copy()

    return P_station_pa



def wind2m(u10):
    """
    Convert 10 m wind speed to 2 m.

    Parameters
    ----------
    u10 : xarray.DataArray
        10 m wind speed.

    Returns
    -------
    u2 : DataArray
         2m based on the full logarithmic wind speed profile equation 33 in Allen et al. 2005
    """

    return u10 * (4.87 / np.log(67.8 * 10.0 - 5.42))


# ---------------------------------------------------------------------------
# FAO-56 PM final equation
# ---------------------------------------------------------------------------

def calculate_FAO56_pmpet(Rn, t, u2, esat, ea, ps, crop="short"):
    """
    Parameters
    ----------
    Rn : xarray.DataArray or float
        Net radiation (MJ/(m2.day))
    t: xarray.DataArray or float
        2m mean temperature (degC)
    u2: xarray.DataArray or float
        2m mean wind speed (m/s)
    esat: xarray.DataArray or float
        saturation vapour pressure (kPa)
    ea: xarray.DataArray or float
        actual vapour pressure (kPa)
    ps: xarray.DataArray or float
        Atmospheric pressure at station elevation (kPa)
        
    h_m : xarray.DataArray, float, or int
        Station elevation in meters

    Returns
    -------
    P_station : same type as input
        Atmospheric pressure at station elevation in Pascals (Pa)

    Notes
    -----
    Calculate FAO-56 PM reference ET₀ (mm/day).
    Based on Allen et al. (2005) Eq. 6.

    crop = "short" (0.12 m) or "tall" (0.50 m)
    """
    # G: ground heat flux  [MJ/(m2.day)], on daily time scales, G ≈ 0 [FAO56 Eq. 42]
    # Cn: numerator crop constant [K.m.s3/(kg.day)]
    # Cd: denominator crop constant [s/m]
    if crop == "short":
        G = 0.0
        Cn = 900.0
        Cd = 0.34
    elif crop == "tall":
        G = 0.0
        Cn = 1600.0
        Cd = 0.38
    else:
        raise ValueError("crop must be 'short' or 'tall'.")

    delta = 2503. * np.exp(17.27 * t / (t + 237.3)) / (t + 237.3)**2.
    gamma = 0.000665 * ps * 1e-3 #[kPa/°C]

    # ASCE-EWRI ET0 [mm/day] = f(0.408 [mm.m2/J], delta [kPa/K],
    # gamma [Pa/K], u2 [m/sec], Rn [J/m2/day], G [J/m2/day], T [degC],
    # esat [Pa], ea [Pa], Cn_alf [K.mm.sec3/Mg/day], Cd_alf [sec/m])
    ET0 = (0.408 * delta * (Rn - G) + gamma * (Cn / (t + 273.)) * u2 * (esat - ea)) / (delta + gamma * (1. + Cd * u2))

    # Attributes
    ET0.attrs.update({"long_name": f"FAO-56 {crop} reference crop ET [mm/day]",
                       "standard_name": "reference_ET",
                       "units": "mm day-1",
                       "cell_methods": "time: mean (interval: 1D)",
                       "grid_mapping": "crs"}
                    )
    
    # List of coordinates to drop    
    coords_to_drop = ['height', 'crs', 'level_height', 'model_level_number', 'sigma']
    # Drop coordinates only if they are present in the dataset
    ET0 = ET0.drop([coord for coord in coords_to_drop if coord in ET0.coords]).rename('ET0').to_dataset()

    return ET0
