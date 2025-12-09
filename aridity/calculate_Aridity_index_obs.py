#< import modules
import xarray as xr
import numpy as np
import cmdline_provenance as cmdprov
import git
import dask.distributed
from dask.distributed import Client
import tempfile
from dask.diagnostics import ProgressBar
import lib_david
import argparse
import os
import sys
import glob
import subprocess

# Ignore warnings
import warnings
import logging
warnings.filterwarnings('ignore') 
logging.getLogger("distributed.worker.memory").setLevel(logging.ERROR)

#################< Functions ########################
def calc_AI(ds_pr,ds_e0,syear,eyear):
    # Open data since it's a single file and resample annually
    # rechunk
    ds_e0 = ds_e0.chunk({'time':12,'lat':'auto','lon':'auto'})
    ds_pr = ds_pr.chunk({'time':12,'lat':'auto','lon':'auto'})

    da_e0 = ds_e0['e0'].resample(time='Y').sum('time')
    da_pr = (ds_pr['precip']).resample(time='Y').sum('time')
   
    #Compute AI
    computed_AI = (da_pr/da_e0).astype("float32").rename("AI")

    # add attributes
    computed_AI.attrs['units']         = ''
    computed_AI.attrs['long_name']     = 'Aridity index, ratio of precipitation to (potential)evepotranspiration'
    computed_AI.attrs['standard_name'] = 'Aridity Index'
    computed_AI.attrs['description']   = 'Aridity index calculated from monthly AWAP v2.0.2 precipitation data (/g/data/zv2/agcd/v2-0-1/precip/total/r005/01month/) and monthly aggregated (using CDO monsum) daily AWRA v7 potential evapotranspiration data (/g/data/fj8/BoM/AWRA/DATA/v7/simulation).'
    ds = computed_AI.to_dataset(name = 'AI')

    return computed_AI
    
def get_git_hash():
    """Returns the git hash for the working repository"""
    git_repo = git.Repo(sys.argv[0], search_parent_directories=True)
    git_root = git_repo.git.rev_parse("--show-toplevel")
    git_hash = git.Repo(git_root).heads[0].commit
    git_text = " (Git hash: %s)" %(str(git_hash)[0:7])
    return git_text

################< Main ###############################


def main(inargs):
    """Calculate the aridity index from gridded National Hydrological Projections (NHP1.0) data"""
    print("Starting...")
    
    dask.config.set({
        "distributed.scheduler.worker-saturation": 1, #< This should use the new behaviour which helps with memory pile up
        })
    client = Client(n_workers=inargs.nworkers, threads_per_worker=1, local_directory = tempfile.mkdtemp()) #need to make this dynamic depending on cpus/memory requested

    #< Keep attributes from input datasets
    xr.set_options(keep_attrs = True)

    #< Set time period
    start_y = inargs.sYear
    end_y = inargs.eYear

    file_name_ann = f"{inargs.outputDir}/AI-atmospheric-based_obs-AUS-5_{inargs.obsData}_{start_y}_{end_y}.nc"
    print(f'Computing {file_name_ann}...')

    if os.path.exists(file_name_ann):
        print(f"File for this period exists already in {inargs.outputDir}. Aborting process.")
    else:
        
        path_awap_pr = "/g/data/zv2/agcd/v2-0-1/precip/total/r005/01month/agcd_v2-0-1_precip_total_r005_monthly_*.nc"
        path_awra_e0 = "/scratch/mn51/dh4185/e0_avg_*_monsum.nc"
    
        # Loop through each year in the specified period
        for year in range(start_y, end_y + 1):
            # Build the file pattern for each year
            file_pattern = f'/scratch/mn51/dh4185/e0_avg_{year}_monsum.nc'
            
            # Check if the file exists
            if not glob.glob(file_pattern):
                print(f"Error: File for year {year} is missing.")
                print(f"Suggested action: Run the script /home/565/dh4185/mn51-dh4185/resample_daily_AWRA-L_e0_to_monthly.sh")
                
                # Automatically running the script
                if inargs.runBashScript:
                    try:
                        # Run the bash script to generate the missing data
                        subprocess.run(['/home/565/dh4185/mn51-dh4185/resample_daily_AWRA-L_e0_to_monthly.sh'], check=True)
                        print(f"Script has been executed successfully. Rerun .py script to compute AI.")
                    except subprocess.CalledProcessError as e:
                        print(f"An error occurred while running the script: {e}")
                    break  # Stop processing after running the script
                else:
                    print("Processing will stop. Please run the script manually.")
                    break  # Stop processing if user doesn't want to run the script
        else:
            # If all files are found
            print(f"All data files for the period {start_y} to {end_y} are present.")
                
        #< Open datasets
        def _preprocess(ds):
            return ds.sel(time=slice(str(start_y),str(end_y)))
            
        ds_awap_pr = xr.open_mfdataset(path_awap_pr, combine='nested',concat_dim='time',parallel=True, preprocess=_preprocess)
        ds_awra_e0 = xr.open_mfdataset(path_awra_e0, combine='nested',concat_dim='time',parallel=True, preprocess=_preprocess).rename({'latitude':'lat', 'longitude':'lon'})
    
        #< Modify grids to have same size and order
        ds_awap_pr = ds_awap_pr.sel(lat=slice(ds_awra_e0.lat[-1],ds_awra_e0.lat[0]))
        ds_awap_pr = ds_awap_pr.sel(lon=slice(ds_awra_e0.lon[0],ds_awra_e0.lon[-1]))
        ds_awra_e0 = ds_awra_e0.sortby(ds_awra_e0.lat)
        ds_awra_e0['lat'] = ds_awap_pr.lat
        ds_awra_e0['lon'] = ds_awap_pr.lon
        
        ds_AI_ann = calc_AI(ds_awap_pr,ds_awra_e0,start_y,end_y)
        ds_AI_ann.attrs['history'] = cmdprov.new_log(extra_notes=[get_git_hash()])
        
        #< Save output
        saver_ann = ds_AI_ann.to_netcdf(file_name_ann,compute=False)
        future_ann = client.persist(saver_ann)
        dask.distributed.progress(future_ann)
        future_ann.compute()


    #< Close the client
    client.close()


if __name__ == '__main__':
    extra_info =""" 
author:
    David Hoffmann, david.hoffmann@bom.gov.au
    Created: 08/01/2025
"""
    description = """
    Calculate the aridity index from monthly AWAP v2.0.2 precipitation data (/g/data/zv2/agcd/v2-0-1/precip/total/r005/01month/) and monthly aggregated (using CDO monsum) daily AWRA v7 potential evapotranspiration data (/g/data/fj8/BoM/AWRA/DATA/v7/simulation).     
        """    
    parser = argparse.ArgumentParser(description=description,
                                     epilog=extra_info,
                                     argument_default=argparse.SUPPRESS,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--sYear", type=int, default=1960, help="Choose start year for calculation.")
    parser.add_argument("--eYear", type=int, default=2014, help="Choose end year for calculation.")
    parser.add_argument("--obsData", type=str, default='AWAPv2-0-1-AWRAv7', choices=['AWAPv2-0-1-AWRAv7', 'ERA5'], help="Choose either 'AWAPv2-0-1-AWRAv7', 'ERA5'. Default is 'AGCD/AWRA'")
    # parser.add_argument("--index", type=str, choices=['atmospheric-based','plant-based'], help="Choose either 'atmospheric-based' or 'plant-based' aridity index.")
    parser.add_argument("--outputDir", type=str, default='/g/data/ia39/ncra/drought_aridity/ai/observations/', help="Output directory on Gadi. Default is'/g/data/ia39/ncra/drought_aridity/ai/'")
    parser.add_argument("--nworkers", type=int, default=8, help="Number of workers in dask distributed client. Default: 8")
    parser.add_argument("--runBashScript", type=bool, default=True, help="Choose 'True' to run bash script with CDO to produce missing annual files for PET. Default: True")

    args = parser.parse_args()
    main(args)