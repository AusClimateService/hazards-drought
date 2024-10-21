#< import modules
import xarray as xr
import numpy as np
import cmdline_provenance as cmdprov
import git
import dask.distributed
from dask.distributed import Client
import tempfile
from dask.diagnostics import ProgressBar
import argparse
import os
import sys
import re
import xclim.indices

sys.path.append('/g/data/mn51/users/dh4185/hazards-drought/aridity/')
import lib_david

# Ignore warnings
import warnings
import logging
warnings.filterwarnings('ignore') 
logging.getLogger("distributed.worker.memory").setLevel(logging.ERROR)

#################< Functions ########################

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


# def get_git_hash():
#     """Returns the git hash for the working repository"""
#     git_repo = git.Repo(sys.argv[0], search_parent_directories=True)
#     git_root = git_repo.git.rev_parse("--show-toplevel")
#     git_hash = git.Repo(git_root).heads[0].commit
#     git_text = " (Git hash: %s)" %(str(git_hash)[0:7])
#     return git_text

################< Main ###############################

def main(inargs):
    """Calculate the Keetch-Byram Drought Index (KBDI)"""
    print("Starting...")
    
    dask.config.set({
        "distributed.scheduler.worker-saturation": 1, #< This should use the new behaviour which helps with memory pile up
        })
    client = Client(n_workers=inargs.nworkers, threads_per_worker=1, local_directory = tempfile.mkdtemp()) #need to make this dynamic depending on cpus/memory requested

    #< Keep attributes from input datasets
    xr.set_options(keep_attrs = True)

    RCM = inargs.RCM
    #MPI-ESM1-2-HR // NorESM2-MM for BARPA only and CNRM-ESM2-1 for CCAM only
    model_list = ['ACCESS-CM2']#,'CMCC-ESM2', 'ACCESS-ESM1-5', 'EC-Earth3', 'CESM2'] 
    # model_list = model_list + ['CNRM-ESM2-1'] if RCM == 'CCAM-v2203-SN' else model_list + ['MPI-ESM1-2-HR', 'NorESM2-MM']


    for model in model_list:
        print('========= '+RCM+'_'+model+' =========')
        bc_string = '_ACS-{}-{}-{}-{}.nc'.format(inargs.bcMethod, inargs.bcSource, '1960' if inargs.bcSource == 'AGCD' else '1979', '2022') if inargs.bc == 'output' else '.nc'
        variant_id = lib_david.data_source['CMIP6'][model]['variant-id']
        file_name = "{}/KBDI_{}_{}_{}_{}_{}_{}_{}{}".format(inargs.outputDir, inargs.index,'AGCD-05i',model,'ssp370',variant_id,'BOM' if RCM == 'BARPA-R' else 'CSIRO','v1-r1','_raw.nc' if inargs.bc == 'raw' else bc_string)

        ###############################  Compute KBDI ############################################

        if os.path.exists(file_name)==False:
            print(f"Computing annual time series. File: {file_name}...")
        
            # read input data for AI calculation
            syear = 1970
            eyear = 2014

            if inargs.index == 'pet_thornthwaite':
                pet_method = "Using Thornthwaite method to estimate potential evapotranspiration."
                input_array_p, var_pr = lib_david.load_target_files('var_p', RCM, model, bc=inargs.bc, bc_method=inargs.bcMethod, bc_source=inargs.bcSource)
                input_array_tmax, var_tmax = lib_david.load_target_files('var_tmax', RCM, model, bc=inargs.bc, bc_method=inargs.bcMethod, bc_source=inargs.bcSource)
      
                print("================ input_files_p =============")
                pr_files = filter_files_by_year(input_array_p, start_year=syear, end_year=eyear)
                print(var_pr)
                print(pr_files)
                # print("================ input_files_p_annual =============")
                # pr_files_annual = filter_files_by_year(input_array_p, start_year=1971, end_year=2000)
                # print(pr_files_annual)
                print("================ input_array_tmax =============")
                tmax_files = filter_files_by_year(input_array_tmax, start_year=syear, end_year=eyear)
                print(var_tmax)
                print(tmax_files)            
                    
                # Calculate KBDI
                print("================ Calculate KBDI =============")

                tmax = xr.open_mfdataset(tmax_files,combine='nested', concat_dim='time', parallel=True).sel(lat=slice(-29,-28),lon=slice(149,150))[var_tmax]#.chunk({"time": -1,"lat":10,"lon":10})
                pr = xr.open_mfdataset(pr_files,combine='nested', concat_dim='time', parallel=True).sel(lat=slice(-29,-28),lon=slice(149,150))[var_pr]#.chunk({"time": -1,"lat":10,"lon":10})

                if inargs.bc == 'raw':
                    pr = (pr * 86400).sel(lat=slice(-44,-10),lon=slice(112,157))
                    tmax = (tmax - 273.15).sel(lat=slice(-44,-10),lon=slice(149,157))
  
                pr_annual = pr.sel(time=slice('1971','2000')).resample(time='YE').sum('time').mean('time')
                pr_annual = pr_annual.compute()
                pr = pr.chunk({"time":-1,"lat":21,"lon":21})  # aim for chunks ~4MB
                tmax = tmax.chunk({"time":-1,"lat":21,"lon":21})
                print(pr_annual)
                print(pr)
                print(tmax)
                
                kbdi = xclim.indices.fire._ffdi.keetch_byram_drought_index(pr, tmax, pr_annual, kbdi0=None).rename("KBDI").compute()
                print(kbdi)
            
                # Save the computed KBDI data to a compressed NetCDF file
                compression = dict(zlib=True, complevel=1)  # Enable zlib compression with level 5 (balance between size and speed)
                
                # Write to NetCDF with compression
                saver = kbdi.to_netcdf(f'{inargs.outputDir}kbdi_output_compressed.nc', encoding={'KBDI': compression}, compute=False)
                future = client.persist(saver)
                dask.distributed.progress(future)
                future.compute()
                
        #     elif inargs.index == 'pet_model':
        #         pet_method = "Using potential evapotranspiration directly from model output."
        #         input_array_p = lib_david.load_target_variable('var_p', RCM, model, bc=inargs.bc, bc_method=inargs.bcMethod, bc_source=inargs.bcSource)
        #         input_array_pet = lib_david.load_target_variable('var_pet', RCM, model, bc=inargs.bc, bc_method=inargs.bcMethod, bc_source=inargs.bcSource)

        #         if inargs.bc == 'raw':
        #             input_array_p = input_array_p * 86400
        #             input_array_pet = input_array_pet * 86400
                    
        #         print("================ input_array_p =============")
        #         print(input_array_p)
        #         print("================ input_array_pet =============")
        #         print(input_array_pet)

        #         # Calculate KBDI
        #         print("================ Calculate KBDI =============")
        #         kbdi = calculate_kbdi(pet_data=input_array_pet, precip_data=input_array_p, init_kbdi=None)
            
            
        #     kbdi.attrs['description'] = f"Keetch-Byram Drought Index (KBDI). Produced for ACS."
        #     kbdi.attrs['method']  = pet_method
        #     # kbdi.attrs['history'] = cmdprov.new_log(extra_notes=[get_git_hash()])
        #     kbdi.attrs['bias correction'] = f"None. Raw {RCM} output data." if inargs.bc == 'raw' else f"{RCM} {bc_string} data";
        #     print(kbdi)        

        #     # Save the computed KBDI data to a compressed NetCDF file
        #     compression = dict(zlib=True, complevel=1)  # Enable zlib compression with level 5 (balance between size and speed)
                
        #     # Write to NetCDF with compression
        #     # kbdi.to_netcdf('kbdi_output_compressed.nc', encoding=encoding)
        #     saver = kbdi.to_netcdf('kbdi_output_compressed.nc', encoding={'KBDI': compression}, compute=False)
        #     future = client.persist(saver)
        #     dask.distributed.progress(future)
        #     future.compute()

        # else:
        #     print("{name} exists. Pass.".format(name=file_name_ann))
           


    #< Close the client
    client.close()


if __name__ == '__main__':
    extra_info =""" 
author:
    David Hoffmann, david.hoffmann@bom.gov.au
    Created: 10/10/2024
"""
    description = """
    Calculate the KBDI.     
        """    
    parser = argparse.ArgumentParser(description=description,
                                     epilog=extra_info,
                                     argument_default=argparse.SUPPRESS,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--RCM", type=str, choices=['BARPA-R', 'CCAM-v2203-SN'], help="Choose RCM either 'BARPA-R' or 'CCAM-v2203-SN'.")
    parser.add_argument("--bc", type=str, choices=['raw','input','output'], help="Choose either 'raw' (py18/hq89 raw BARPA/CCAM), 'input' (ia39 bc input, 5km) or 'output' (ia39 bc output, 5km).")
    parser.add_argument("--bcSource", type=str, default='AGCD', choices=['AGCD', 'BARRA'], help="Choose either 'AGCD', 'BARRA'. Default is 'AGCD'")
    parser.add_argument("--bcMethod", type=str, default='QME', choices=['QME', 'MRNBE'], help="Choose either 'MRNBC', 'QME'. Default is 'QME'")
    parser.add_argument("--index", type=str, choices=['pet_model','pet_thornthwaite'], help="Choose either 'pet_model' or 'pet_thornthwaite' basedy index.")
    parser.add_argument("--outputDir", type=str, default='/g/data/ia39/ncra/drought_aridity/kbdi/acs_downscaled_notBC_5km/', help="Output directory on Gadi. Default is'/g/data/ia39/ncra/drought_aridity/kbdi/'")
    parser.add_argument("--nworkers", type=int, default=5, help="Number of workers in dask distributed client.")

    args = parser.parse_args()
    main(args)
    

