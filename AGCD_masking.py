#< import libraries
import xarray as xr
import os
import dask
from dask.distributed import Client
import tempfile
import sys
import numpy as np
import git
import argparse
import cmdline_provenance as cmdprov
from dask.diagnostics import ProgressBar



#< Functions                   
def create_agcd_q_mask(ClimStartYr,ClimEndYr):
    """ Creates a mask by having a weight of at least 1 (minimal influence to influence of nearby stations on climatology for grid cell) for at least 80% of days in the climatological period from ClimStartYr to ClimEndYr.
    Args:
        ClimStartYr
        ClimEndYr
    Returns:
        xarray of fractions
    """
    print("Opening files.")
    da_weights            = xr.open_mfdataset('/g/data/zv2/agcd/v1-0-2/precip/weight/r005/01day/agcd_v1_precip_weight_r005_daily_*.nc').sel(time=slice(ClimStartYr,ClimEndYr))
    da_weights_count      = (da_weights['weight'] >= 1.0).sum('time', keep_attrs=True)
    ds                    = (da_weights_count / len(da_weights['time'])).to_dataset(name='fraction')
    ds.attrs["comment"]   = f"The amount of influence that nearby station measurements have on a grid point (1 means 'minimal influence' and 3 means 'influenced by nearby stations'), set to 0 when a grid point has not been changed from its base climatology by nearby station measurements. This mask shows where nearby stations influenced the data, meaning a value of at least 1. It counts the fraction of days where the weight value is >=1 in the period {ClimStartYr}-{ClimEndYr}. When applying the mask, aim for a fraction of 0.8 or greater."
    ds.attrs['history']   = cmdprov.new_log(extra_notes=[get_git_hash()])
    ds.attrs['agcd_version']= da_weights.attrs['agcd_version']
    
    return ds
    
def get_git_hash():
    """Returns the git hash for the working repository"""
    git_repo = git.Repo(sys.argv[0], search_parent_directories=True)
    git_root = git_repo.git.rev_parse("--show-toplevel")
    git_hash = git.Repo(git_root).heads[0].commit
    git_text = " (Git hash: %s)" %(str(git_hash)[0:7])
    return git_text


def main(inargs):
    """Create AGCD quality mask for precip observations using weights from AGCD v1-0-2 on Gadi (zv2)"""

    #< Set up dask
    dask.config.set({
        "distributed.scheduler.worker-saturation": 1, #< This should use the new behaviour which helps with memory pile up
        })
    client = Client(n_workers=inargs.nworkers, threads_per_worker=1, local_directory = tempfile.mkdtemp()) #need to make this dynamic depending on cpus/memory requested

    #< Keep attributes from input datasets
    xr.set_options(keep_attrs = True)
    
    print(f'Creating mask for period {inargs.ClimStartYr} to {inargs.ClimEndYr}.')
    
    #< Masking
    mask_file  = f'{inargs.OutputDir}agcd_v1-0-2_precip_weight_r005_daily_{inargs.ClimStartYr}_{inargs.ClimEndYr}fraction_ge1.nc'

    if os.path.exists(mask_file):
        print(f"Mask already exists for period {inargs.ClimStartYr} to {inargs.ClimEndYr} in specified directory. Check {mask_file}")
    else:
        qmask = create_agcd_q_mask(inargs.ClimStartYr,inargs.ClimEndYr)
        
        #< Save output
        saver = qmask.to_netcdf(mask_file,compute=False)
        future = client.persist(saver)
        dask.distributed.progress(future)
        future.compute()

    print(f"Done! Mask file {mask_file} created.")
    
    #< Close the client
    client.close()


if __name__ == '__main__':
    extra_info  =""" author: David Hoffmann, david.hoffmann@bom.gov.au """
    description = """ Creates a mask by having a weight of at least 1 (minimal influence to influence of nearby stations on climatology for grid cell) for at least 80% of days in the climatological period from ClimStartYr to ClimEndYrfrom gridded AGCDv1-0-2 precip weight data from '/g/data/zv2/agcd/v1-0-2/precip/weight/r005/01day/'. Guideline: ACS CMIP6 downscaled data uses 1960-2022 (https://doi.org/10.5194/gmd-17-731-2024), NHP1 data use 1975-2005 for bias correction (http://www.bom.gov.au/research/publications/researchreports/BRR-061.pdf)."""    
    
    parser = argparse.ArgumentParser(description=description,
                                     epilog=extra_info,
                                     argument_default=argparse.SUPPRESS,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
                                     
    parser.add_argument("--ClimStartYr", type=str, help="Calculate climatology from this year")
    parser.add_argument("--ClimEndYr", type=str, help="Calculate climatology to this year")
    parser.add_argument("--OutputDir", type=str, default='/g/data/mn51/users/dh4185/', help="Output directory on Gadi.")
    parser.add_argument("--nworkers", type=int, default=5, help="Number of workers in dask distributed client.")

    args = parser.parse_args()
    main(args)
    



                