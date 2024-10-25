#< import modules
import xarray as xr
import numpy as np
import cmdline_provenance as cmdprov
import git
import tempfile
import argparse
import os
import sys
import re
import xclim.indices
import multiprocessing
import matplotlib.pyplot as plt
import time

sys.path.append('/g/data/mn51/users/dh4185/hazards-drought/aridity/')
from lib_david import *

#< Functions
def filter_files_by_year(file_list, start_year=1995, end_year=2014):
    """
    Filter the list of files to include only those containing specified years
    in the format YYYYMMDD-YYYYMMDD.

    Parameters:
    - file_list: List of file paths to filter
    - start_year: Start year for filtering (default: 1995)
    - end_year: End year for filtering (default: 2014)

    Returns:
    - List of filtered file paths
    """
    # Create a regex pattern to match the date format
    year_pattern = '|'.join(str(year) for year in range(start_year, end_year + 1))
    regex = re.compile(rf'.*?({year_pattern})\d{{2}}\d{{2}}-\d{{4}}\d{{2}}\d{{2}}.*')

    # Filter the files based on the regex
    filtered_files = [file for file in file_list if regex.match(file)]
    
    return filtered_files
    
def _calc_kbdi(pr, tmax, pr_annual):
    return xclim.indices.fire._ffdi.keetch_byram_drought_index(pr, tmax, pr_annual, kbdi0=None).rename("KBDI")

def calc_kbdi(tmax_files, pr_files, var_tmax, var_pr, lats, lons, pr_annual):
    tmax = xr.open_mfdataset(tmax_files,combine='nested', concat_dim='time').sel(lat=lats,lon=lons)[var_tmax].compute()
    pr = xr.open_mfdataset(pr_files,combine='nested', concat_dim='time').sel(lat=lats,lon=lons)[var_pr].compute()
    
    pr_annual_local = pr_annual.sel(lat=lats,lon=lons).compute()
    kbdi = _calc_kbdi(pr,tmax,pr_annual_local)

    tmax.close()
    pr.close()
    pr_annual_local.close()

    tmax = None
    pr = None
    pr_annual_local = None
    
    return kbdi

def get_dim_and_length(file):
    ds = xr.open_dataset(file)
    lats = ds["lat"]
    lons = ds["lon"]
    nlat = len(lats)
    nlon = len(lons)

    return lats, lons, nlat, nlon

def get_batch(x, size):
    if size is None:
        size = len(x)
    return [x[i:i+size] for i in range(0, len(x), size)]

data_source = {
    'ERA5':{'var_p':'tp', 'var_sm':'swvl{}', 'var_pet':'evspsblpot','var_lat':'lat','var_lon':'lon'}, #swvl1+swvl2+swvl3 is volume of water in 1m soil column
    'AGCD':{'var_p':'precip', 'var_lat':'latitude','var_lon':'longitude'},
    'AWRA':{'var_sm':'s{}', 'var_pet':'e0','var_lat':'latitude','var_lon':'longitude'}, #s0 or ss
    'CMIP6':{'var_tmax':'tasmax','var_tmin':'tasmin','var_p':'pr', 'var_sm':'mrsos', 'var_et':'evspsbl', 'var_pet':'evspsblpot',
             'var_lat':'lat','var_lon':'lon',
             'CMCC-ESM2':{'variant-id':'r1i1p1f1','version':'v1'}, 
             'ACCESS-ESM1-5':{'variant-id':'r6i1p1f1','version':'v1'},
             'ACCESS-CM2':{'variant-id':'r4i1p1f1','version':'v1'},
             'EC-Earth3':{'variant-id':'r1i1p1f1','version':'v1'},
             'MPI-ESM1-2-HR':{'variant-id':'r1i1p1f1','version':'v1'}, 
             'CESM2':{'variant-id':'r11i1p1f1','version':'v1'},
             'NorESM2-MM':{'variant-id':'r1i1p1f1','version':'v1'},
             'CNRM-ESM2-1':{'variant-id':'r1i1p1f2','version':'v1'}
        }
    }
 

def main(inargs):
    
    #< Keep attributes from input datasets
    xr.set_options(keep_attrs = True)

    RCM = inargs.RCM
    #MPI-ESM1-2-HR // NorESM2-MM for BARPA only and CNRM-ESM2-1 for CCAM only
    model_list = ['ACCESS-CM2']#,'CMCC-ESM2', 'ACCESS-ESM1-5', 'EC-Earth3', 'CESM2'] 
    # model_list = model_list + ['CNRM-ESM2-1'] if RCM == 'CCAM-v2203-SN' else model_list + ['MPI-ESM1-2-HR', 'NorESM2-MM']


    for model in model_list:
        print('========= '+RCM+'_'+model+' =========')
        bc_string = '_ACS-{}-{}-{}-{}.nc'.format(inargs.bcMethod, inargs.bcSource, '1960' if inargs.bcSource == 'AGCD' else '1979', '2022') if inargs.bc == 'output' else '.nc'
        variant_id = data_source['CMIP6'][model]['variant-id'] # lib_david.data_source
        file_name = "{}/KBDI_{}_{}_{}_{}_{}_{}_{}{}".format(inargs.outputDir, inargs.index,'AGCD-05i',model,'ssp370',variant_id,'BOM' if RCM == 'BARPA-R' else 'CSIRO','v1-r1','_raw.nc' if inargs.bc == 'raw' else bc_string)

        ###############################  Compute KBDI ############################################

        if os.path.exists(file_name)==False:
            print(f"Computing annual time series. File: {file_name}...")
        
            # read input data for AI calculation
            syear = 1960
            eyear = 2000

            if inargs.index == 'pet_thornthwaite':
                pet_method = "Using Thornthwaite method to estimate potential evapotranspiration."
                input_array_p, var_pr = lib_david.load_target_files('var_p', RCM, model, bc=inargs.bc, bc_method=inargs.bcMethod, bc_source=inargs.bcSource)
                input_array_tmax, var_tmax = lib_david.load_target_files('var_tmax', RCM, model, bc=inargs.bc, bc_method=inargs.bcMethod, bc_source=inargs.bcSource)
      
                print("================ input_files_p =============")
                pr_files = filter_files_by_year(input_array_p, start_year=syear, end_year=eyear)
                print(var_pr)
                print(pr_files)

                print("================ input_array_tmax =============")
                tmax_files = filter_files_by_year(input_array_tmax, start_year=syear, end_year=eyear)
                print(var_tmax)
                print(tmax_files)            

            pr_annual = xr.open_mfdataset(pr_files).sel(time=slice('1971','2000')).resample(time='YE').sum('time').mean('time')[var_pr]
            print(" Annual climatology calculated.")
            
            lats, lons, nlat, nlon = get_dim_and_length(pr_files[0])
            time_len = ((eyear-syear)+1)*365
            print(time_len)

            batch_size = int(np.sqrt(inargs.batchSizeMB*1e6/(4*time_len))) #-> for 20MB chunk size
            print(f"batch size: {batch_size}")
            
            lat_batch = get_batch(lats.values,batch_size) #inargs.batchSize)
            lon_batch = get_batch(lons.values,batch_size) #inargs.batchSize)
            
            fun_args = []
            for ilat in range(len(lat_batch)):
                for ilon in range(len(lon_batch)):
                    
                    fun_args.append(
                        [
                            tmax_files,
                            pr_files,
                            var_tmax,
                            var_pr,
                            lat_batch[ilat],
                            lon_batch[ilon],
                            pr_annual
                        ]
                    )  
            
            start = time.time()
            with multiprocessing.Pool(inargs.nworkers) as p:
                kbdi_list = p.starmap(calc_kbdi, fun_args)
            
            print("Multiprocessing done.")

            kbdi = xr.combine_by_coords(kbdi_list)
            kbdi.KBDI = kbdi.KBDI.astype("float32")
            end = time.time()
            print(f"Time needed for calculation: {(end - start)/60} min")

            print(kbdi)
            start = time.time()
            compression = dict(zlib=True, complevel=1)
            kbdi.to_netcdf(f"{inargs.outputDir}kbdi_test.nc", encoding={'KBDI': compression})
            end = time.time()
            print(f"Time needed for writing to disk: {(end - start)/60} min")
            
    return


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(#description=description,
                                     #epilog=extra_info,
                                     argument_default=argparse.SUPPRESS,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--RCM", type=str, choices=['BARPA-R', 'CCAM-v2203-SN'], help="Choose RCM either 'BARPA-R' or 'CCAM-v2203-SN'.")
    parser.add_argument("--bc", type=str, choices=['raw','input','output'], help="Choose either 'raw' (py18/hq89 raw BARPA/CCAM), 'input' (ia39 bc input, 5km) or 'output' (ia39 bc output, 5km).")
    parser.add_argument("--bcSource", type=str, default='AGCD', choices=['AGCD', 'BARRA'], help="Choose either 'AGCD', 'BARRA'. Default is 'AGCD'")
    parser.add_argument("--bcMethod", type=str, default='QME', choices=['QME', 'MRNBE'], help="Choose either 'MRNBC', 'QME'. Default is 'QME'")
    parser.add_argument("--index", type=str, choices=['pet_model','pet_thornthwaite'], help="Choose either 'pet_model' or 'pet_thornthwaite' basedy index.")
    parser.add_argument("--outputDir", type=str, default='/g/data/ia39/ncra/drought_aridity/kbdi/acs_downscaled_BC_5km/', help="Output directory on Gadi. Default is'/g/data/ia39/ncra/drought_aridity/kbdi/acs_downscaled_BC_5km/'")
    parser.add_argument("--nworkers", type=int, default=7, help="Number of workers in dask distributed client.")
    parser.add_argument("--batchSizeMB", type=int, default=50, help="Size in megabyte of batch to pass to workers. Aim for ~5-100MB") # time*lat*lon*4/1e6 =~ 5MB-100MB batch_size = sqrt(20*1e6/(4*time_len)) -> for 20MB chunk size


    args = parser.parse_args()
    main(args)
    
    
    
    
    # pr_annual = pr.sel(time=slice('1971','2000')).resample(time='YE').sum('time').mean('time')
