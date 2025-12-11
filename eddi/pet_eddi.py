
"""
pet_eddi.py
===========
Compute FAO-56 Penman–Monteith reference ET₀ and the Evaporative
Demand Drought Index (EDDI).

Author
------
David Hoffmann (david.hoffmann@bom.gov.au)

Description
-----------
This module exposes two main public functions:

- compute_daily_FAO56_PET(...)
- compute_ref_eddi(...)

Both functions accept xarray objects and return xarray objects, suitable
for large gridded datasets (e.g., BARRA, ERA5, ...).

Dependencies: numpy, xarray, scipy (for rankdata), cmdline_provenance, sys, git
"""

import numpy as np
import xarray as xr
from scipy.stats import rankdata
import cmdline_provenance as cmdprov

import utils

# Global attributes from BARRA-R2/C2 that must be kept
barra_glob_attrs = ["axiom_version","axiom_schemas_version","axiom_schema","productive_version",
              "variable_version","Conventions","activity_id","title","source","project",
              "program","summary","publisher_institution","processing_level","version_realisation",
              "institution_id","source_id","driving_variant_label","driving_experiment_id",
              "driving_institution_id","driving_experiment","driving_source_id","nesting_source_id",
              "project_id","domain","domain_id","standard_name_vocabulary","tracking"]


# ---------------------------------------------------------------------------
# FAO-56 PET wrapper
# ---------------------------------------------------------------------------

def compute_daily_FAO56_PET(
    tasmax, tasmin, rsds, hurs, sfcWind, psl, elev,
    coords=None,
    crop="short",
    freq="day",
    global_attrs=None
):
    """
    Compute daily FAO-56 Penman–Monteith ET₀.
    See documentation inside for a full description.
    """

    print("----- Start calculation ... -----")

    # Convert temperature if needed
    tasmax = utils.convertTemperature(tasmax)
    tasmin = utils.convertTemperature(tasmin)

    # Temperature
    t = utils.meanVar(tasmax, tasmin)

    # Vapour pressures
    esat = utils.satVapourPressure(tasmax=tasmax, tasmin=tasmin)
    ea = utils.actVapourPressure(esat, hurs=hurs)

    # Pressure [kPa]
    ps = (utils.convert_sea_level_pressure_to_station_pressure(psl, elev) * 1000)

    # Slope of vapour pressure curve
    delta = 2503 * np.exp(17.27 * t / (t + 237.3)) / (t + 237.3)**2

    # Radiation: downwelling short wave radiation [J/(m2.day)]
    Rs = rsds * 86400 / 1e6 #[MJ/(m2.day)]

    u2 = utils.wind2m(sfcWind)

    # Day-of-year for radiation geometry
    jday = xr.DataArray(
        rsds.time.dt.dayofyear.values,
        coords=[('day', rsds.time.dt.dayofyear.values)],
        dims=["day"]
    )

    Gsc = 4.92 #[MJ/(m2.hr)]
    Lat_rad = np.deg2rad(tasmax["lat"])
    
    # decl: declination [radians] = f(pi[unitless], jday[unitless]) [FAO56 Eq. 42]
    decl = 0.409 * np.sin((2*np.pi*jday/365) - 1.39)
    arg = (-np.tan(Lat_rad) * np.tan(decl)).clip(-1, 1)
    # ws: sunset hour angle ws[radians] = f(latitude[rad], dec[rad]) [FAO56 Eq. 25]
    ws = np.arccos(arg)
    # dr: inverse relative distance of earth from sun [unitless] [FAO56 Eq. 23]
    dr = 1. + 0.033 * np.cos(2*np.pi*jday/365)

    # Ra: extra-terrestrial (TOA) SW radiation [J/(m2.day)] = 
        #f(Gsc [J/m2/hour], dr [unitless], ws [rad], Lat [rad], decl [rad]) 
    Ra = (
        (24/np.pi) * Gsc * dr *
        (ws * np.sin(Lat_rad) * np.sin(decl) +
         np.cos(Lat_rad) * np.cos(decl) * np.sin(ws))
    )

    if "day" in Ra.dims:
        Ra = Ra.rename({"day": "time"})
    Ra = Ra.assign_coords(time=tasmax["time"])

    elev_m = elev
    # Rso: clear-sky SW radiation at surface [J/m2/day] = f(Ra[J/(m2.day)], elev[m]) [FAO56 Eq. 37]
    Rso = Ra * (0.75 + 2e-5 * elev_m)
    Rso = Rso.reindex_like(Rs).where(Rso > 0)

    RsRso = (Rs / Rso).fillna(0).clip(0.3, 1.0)
    # fcd: cloudiness function [unitless] [FAO56 pg. 79]
    fcd = (1.35 * RsRso - 0.35).clip(0.05, 1.0)

    sigma = 4.901e-9 #[MJ/(K4.m2.day)]
    albedo = 0.23
    
    # Rns: Net shortwave radiation
    Rns = (1 - albedo) * Rs

    # Rnl: net LW radiation upward[J/(m2.day)] = f(sigma [J/(K4.m2.day)], eact[Pa], T[°C])   [FAO56 Eq. 39]
    Rnl = (
        sigma * fcd * (0.34 - 0.14 * (ea)**0.5) * 
        (0.5 * ((tasmax+273.16)**4. + (tasmin+273.16)**4.))
    ) #[MJ/(m2.day)]
    # Rn: net radiation [J/(m2.day)]
    Rn = Rns - Rnl

    # Final PM equation
    ET = utils.calculate_FAO56_pmpet(Rn, t, u2, esat, ea, ps, crop)
    ET = ET.astype("float32")

    # Monthly mean if needed
    if freq == "month":
        ET = ET.resample(time="ME").mean()

    # Metadata
    my_log = cmdprov.new_log()
    if global_attrs: 
        # Replace with the filtered attributes
        ET.attrs = global_attrs
        
    ET.attrs.update({
        "description": "FAO-56 crop reference evapotranspiration (ETo) according to Allen et al. 2005 (https://www.mesonet.org/images/site/ASCE_Evapotranspiration_Formula.pdf). Produced by the Australian Climate Service (ACS) Drought and Aridity Hazard Team."
    })
    ET.attrs.update({'history': f"{cmdprov.new_log()}"})
    # ET.attrs.update({'history': f"{cmdprov.new_log(extra_notes=[utils.get_git_hash()])} {ET.attrs['history']}"})

    return ET


# ---------------------------------------------------------------------------
# EDDI — Evaporative Demand Drought Index
# ---------------------------------------------------------------------------

def compute_EDDI(
    ds_ET0,
    ndays,
    compute_categories=True,
    outdir="EDDI_output",
    filename="EDDI",
    enc = None,
    global_attrs=None
):
    import os
    import numpy as np
    import xarray as xr
    from dask.diagnostics import ProgressBar

    # ---------------------------------------------------
    # Step 1 — Rechunk for efficient rolling
    # ---------------------------------------------------
    ds_ET0 = ds_ET0.chunk({"time": -1, "lat": 100, "lon": 100})

    # Remove Feb 29 completely.
    # Because:
    #     It keeps DOYs aligned
    #     It avoids DOY=60 having 25% more samples
    #     It eliminates shifts in all DOYs after Feb 29
    ds_ET0 = ds_ET0.where(~((ds_ET0['time.month'] == 2) & (ds_ET0['time.day'] == 29)), drop=True)

    # ---------------------------------------------------
    # Step 2 — Rolling n-day cumulative ET
    # ---------------------------------------------------
    ds_ET0_agg = (
        ds_ET0.rolling(time=ndays, center=False)
        .sum()
        .dropna(dim="time")
    )

    ds_ET0_agg = ds_ET0_agg.persist()

    # ---------------------------------------------------
    # Step 3 — Assign DOY and year
    # ---------------------------------------------------
    da = ds_ET0_agg

    da = da.assign_coords(
        year=da['time'].dt.year,
        doy=da['time'].dt.dayofyear
    )

    # ---------------------------------------------------
    # Step 4 — Rank by DOY (vectorised, fast version)
    # ---------------------------------------------------
    from scipy.stats import rankdata

    def rank_by_doy(da):
        """
        Rank EDDI values within each DOY across all years.
        Returns ranks 1..N where N = number of years in record.
        """
        # ---------------------------------------------------
        # Ensure correct chunking for rankdata(axis=0)
        # ---------------------------------------------------
        if da.chunks is None:
            # not yet chunked at all
            da = da.chunk({"time": -1, "lat": 100, "lon": 100})
        else:
            # find which axis corresponds to time
            time_dim = da.get_axis_num("time")
            time_chunks = da.chunks[time_dim]
        
            # if time axis has more than 1 chunk → must rechunk
            if len(time_chunks) > 1:
                da = da.chunk({"time": -1, "lat": 100, "lon": 100})
    
        # add year/doy if absent
        if "doy" not in da.coords:
            da = da.assign_coords(doy=da.time.dt.dayofyear)

        # group by day of year
        groups = da.groupby("doy")
    
        def rank_func(x):
            # x has shape (time, lat, lon) for each DOY group
            return xr.apply_ufunc(
                rankdata,
                x,
                input_core_dims=[["time"]],
                output_core_dims=[["time"]],
                vectorize=True,
                kwargs={"axis": 0, "method": "dense"},
                dask="parallelized",
                output_dtypes=[x.dtype],
            )
    
        ranked = groups.map(rank_func)
        ranked_sorted = ranked.sortby("time")
    
        return ranked_sorted
    da_ranked = rank_by_doy(da)

    # ---------------------------------------------------
    # Step 5 — EDDI computation
    # ---------------------------------------------------

    climstart = int(da.year.min())
    climend = int(da.year.max())
    N = climend - climstart + 1
    den = N + 0.33

    rank = da_ranked

    P = (rank - 0.33) / den

    # W transform
    W_lower = np.sqrt(-2 * np.log(P))
    W_upper = np.sqrt(-2 * np.log(1 - P))

    W = xr.where(P <= 0.5, W_lower, W_upper)

    # Constants
    C0, C1, C2 = 2.515517, 0.802853, 0.010328
    d1, d2, d3 = 1.432788, 0.189269, 0.001308

    # EDDI
    EDDI = xr.where(
        P <= 0.5,
        -1.0 * (W - (C0 + C1 * W + C2 * W**2) /
                (1 + d1 * W + d2 * W**2 + d3 * W**3)),
        (W - (C0 + C1 * W + C2 * W**2) /
         (1 + d1 * W + d2 * W**2 + d3 * W**3))
    )

    # Persist full EDDI time series before writing yearly files
    EDDI = EDDI.persist()

    # ---------------------------------------------------
    # Step 6 — Optional categories
    # ---------------------------------------------------
    if compute_categories:
        brks = [
            -2.05374891,
            -1.6448363,
            -1.28155157,
            -0.8416212,
            -0.52440051,
            0.52440051,
            0.8416212,
            1.28155157,
            1.6448363,
            2.05374891,
        ]

        vals = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]

        EDDI_cat = xr.full_like(EDDI, -9999)

        for low, high, val in zip(brks[:-1], brks[1:], vals[1:-1]):
            EDDI_cat = xr.where((EDDI > low) & (EDDI <= high), val, EDDI_cat)

        EDDI_cat = xr.where(EDDI <= brks[0], vals[0], EDDI_cat)
        EDDI_cat = xr.where(EDDI > brks[-1], vals[-1], EDDI_cat)

        EDDI_cat = EDDI_cat.persist()

    # ---------------------------------------------------
    # Step 7 — Write yearly output to avoid one massive file
    # ---------------------------------------------------
    os.makedirs(outdir, exist_ok=True)

    years = np.unique(da.year.values)

    for y in years:
        print(f"Writing {y}...")

        ed = EDDI.sel(time=str(y)).drop_vars(['year', 'doy'], errors="ignore")
        ed = ed.transpose("time", "lat", "lon")
        ed.attrs.update({"long_name": "Evapoarative demand drought index",
                       "standard_name": "EDDI",
                       "units": "",
                       "cell_methods": "time: mean (interval: 1D)",
                       "grid_mapping": "crs"}
                        )

        if compute_categories:
            edc = EDDI_cat.sel(time=str(y)).drop_vars(['year', 'doy'], errors="ignore")
            edc = edc.transpose("time", "lat", "lon")
            edc.attrs.update({"long_name": "USDM categorised evapoarative demand drought index",
                       "standard_name": "EDDI_cat",
                       "units": "",
                       "cell_methods": "time: mean (interval: 1D)",
                       "grid_mapping": "crs"}
                        )

            ds_out = xr.Dataset({"EDDI": ed, "EDDI_cat": edc})
        else:
            ds_out = xr.Dataset({"EDDI": ed})

        # Metadata
        my_log = cmdprov.new_log()
        if global_attrs: 
            # Replace with the filtered attributes
            ds_out.attrs = global_attrs
        
        ds_out.attrs.update({"description": "Evaporative demand drought index (EDDI, Hobbins et al. 2016, https://doi.org/10.1175/JHM-D-15-0121.1) based on FAO-56 crop reference evapotranspiration (ETo) according to Allen et al. 2005 (https://www.mesonet.org/images/site/ASCE_Evapotranspiration_Formula.pdf). Produced by the Australian Climate Service (ACS) Drought and Aridity Hazard Team."}
                           )
        ds_out.attrs.update({'history': f"{cmdprov.new_log()}"})
        # ds_out.attrs.update({'history': f"{cmdprov.new_log(extra_notes=[utils.get_git_hash()])} {ds_out.attrs['history']}"})

        ds_out.to_netcdf(f"{outdir}/{filename}_{y}.nc", encoding= enc or {})

    print("Done.")

    return ds_out

