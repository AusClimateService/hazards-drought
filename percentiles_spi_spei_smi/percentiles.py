'''
######################################################################
## BUREAU OF METEOROLOGY
## AUSTRALIAN CLIMATE SERVICE
## DROUGHT HAZARD TEAM
##
## DATE:         June-2024
## SCRIPT:       percentiles.py
## AUTHOR:       jessica.bhardwaj@bom.gov.au
##
## PURPOSE:      Script to generate rainfall percentiles for ACS work package 3 deliverables.
##
######################################################################
'''

#< import modules
import xarray as xr
import numpy as np
import cmdline_provenance as cmdprov
import git
import dask.distributed
from dask.distributed import Client
import tempfile
from datetime import datetime
import argparse
import os
import sys
from scipy.stats import gamma, norm
import dask.array as da
import utils

# Import GWL script
sys.path.append('/g/data/mn51/users/jb6465/drought-github/submodules/gwls')
from gwl import get_GWL_timeslice

# Ignore warnings
import warnings
import logging
warnings.filterwarnings('ignore') 
logging.getLogger("distributed.worker.memory").setLevel(logging.ERROR)

def compute_percentile(base_period, percentile_threshold):
    
    percentile_array = base_period.quantile(percentile_threshold/100, dim='time')

    return percentile_array


def get_git_hash():
    """Returns the git hash for the working repository"""
    git_repo = git.Repo(sys.argv[0], search_parent_directories=True)
    git_root = git_repo.git.rev_parse("--show-toplevel")
    git_hash = git.Repo(git_root).heads[0].commit
    git_text = " (Git hash: %s)" %(str(git_hash)[0:7])
    return git_text

################< Main ###############################

def main(inargs):
    """Calculate the Standardised Precipitation Index"""
    
    dask.config.set({
        'array.chunk-size': "256 MiB",
        'array.slicing.split_large_chunks': True, #This can make AXIOM very slow
        'distributed.comm.timeouts.connect': '120s',
        'distributed.comm.timeouts.tcp': '120s',
        'distributed.comm.retry.count': 10,
        'distributed.scheduler.allowed-failures': 20,
        "distributed.scheduler.worker-saturation": 1.1, #This should use the new behaviour which helps with memory pile up
        })
    client = Client(n_workers=10, threads_per_worker=1, local_directory = tempfile.mkdtemp(), memory_limit = "63000mb", dashboard_address=":8787") #need to make this dynamic depending on cpus/memory requested
    print("Dask dashboard is available at:", client.dashboard_link)

    
    # compute spi for each model
    for RCM in ['BARPA-R', 'CCAM-v2203-SN']:
        #MPI-ESM1-2-HR // NorESM2-MM for BARPA only and CNRM-ESM2-1 for CCAM only
        model_list = ['CMCC-ESM2', 'ACCESS-ESM1-5', 'ACCESS-CM2', 'EC-Earth3', 'CESM2'] 
        model_list = model_list + ['CNRM-ESM2-1'] if RCM == 'CCAM-v2203-SN' else model_list + ['MPI-ESM1-2-HR', 'NorESM2-MM']
        
        for model in model_list:
            print('========= '+RCM+'_'+model+' =========')
            bc_string = '_ACS-{}-{}-{}-{}.nc'.format(inargs.bcMethod, inargs.bcSource, '1960' if inargs.bcSource == 'AGCD' else '1979', '2022')
            variant_id = utils.data_source['CMIP6'][model]['variant-id']
            file_name = "/scratch/mn51/jb6465/p{}_{}month_{}_{}_{}_{}_{}_{}_GWL{}{}".format(inargs.percentileThreshold, inargs.Accumulation,'AGCD-05i',model,'ssp370',variant_id,'BOM' if RCM == 'BARPA-R' else 'CSIRO','v1-r1', inargs.GWL,'.nc' if inargs.bc == 'input' else bc_string)
            
            if os.path.exists(file_name)==False:
                print("Computing {name}...".format(name=file_name))
    
                # Group by month and apply the SPI calculation
                input_array = utils.load_target_variable('var_p', RCM, model, inargs.Accumulation, bc=inargs.bc, bc_method=inargs.bcMethod, bc_source=inargs.bcSource)
                
                input_array = get_GWL_timeslice(input_array,'CMIP6',model,variant_id,'ssp370',inargs.GWL)
                input_array = input_array.astype(np.float32).chunk({'time':-1})
                ds_perc = input_array.groupby('time.month').map(lambda x: compute_percentile(x, inargs.percentileThreshold))
                
                ds_perc = ds_perc.rename('p{}_{}month'.format(inargs.percentileThreshold, inargs.Accumulation))

                ds_perc.attrs['description'] = f'Rainfall percentile threshold calculated for each GWL base period. Further details in supporting technical documentation.'
                ds_perc.attrs['created'] = (datetime.now()).strftime("%d/%m/%Y %H:%M:%S")    
                ds_perc.attrs['history'] = cmdprov.new_log(extra_notes=[get_git_hash()])
                               
                
                #< Save output
                saver = ds_perc.to_netcdf(file_name,compute=False)
                future = client.persist(saver)
                dask.distributed.progress(future)
                future.compute()
            else:
                print("{name} exists. Pass.".format(name=file_name))

    #< Close the client
    client.close()


if __name__ == '__main__':
    extra_info =""" 
author:
    Jessica Bhardwaj, jessica.bhardwaj@bom.gov.au
    Created: 19/06/2024
    Modified: --
"""
    description = """
    Calculate rainfall percentiles for ACS work package 3 deliverables.     
        """    
    parser = argparse.ArgumentParser(description=description,
                                     epilog=extra_info,
                                     argument_default=argparse.SUPPRESS,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
                                     
    parser.add_argument("--GWL", type=str, choices=['1.0', '1.2', '1.5', '2.0', '3.0', '4.0'], help="Specify GWL")
    parser.add_argument("--bc", type=str, choices=['raw','input','output'], help="Choose either 'raw' (py18/hq89 raw BARPA/CCAM), 'input' (ia39 bc input, 5km) or 'output' (ia39 bc output, 5km).")
    parser.add_argument("--bcSource", type=str, default='AGCD', choices=['AGCD', 'BARRA'], help="Choose either 'AGCD', 'BARRA'. Default is 'AGCD'")
    parser.add_argument("--bcMethod", type=str, default='QME', choices=['QME', 'MRNBE'], help="Choose either 'MRNBC', 'QME'. Default is 'QME'")
    parser.add_argument("--percentileThreshold", type=int, default=15, help="Specify rainfall percentile threshold. Default is 15 as it corresponds to SPI=-1.")
    parser.add_argument("--Accumulation", type=int, default=3, help="Choose accumulation i.e. 1-month, 3-months, 6-months, etc. Default is 3")
    parser.add_argument("--basePeriodStart", type=int, default=1965, help="Specify start year for base period. Minimum 30 years advised. 50 years helps avoid NaN errors on upper tails. Default = 1965.")
    parser.add_argument("--basePeriodEnd", type=int, default=2014, help="Specify end year for base period. Minimum 30 years advised. 50 years helps avoid NaN errors on upper tails. Default = 2014")
    parser.add_argument("--outputDir", type=str, default='/g/data/mn51/projects/work_package_4/climate_hazard_indices/drought/', help="Output directory on Gadi. Default is'/g/data/mn51/projects/work_package_4/climate_hazard_indices/drought/'")
    parser.add_argument("--nworkers", type=int, default=10, help="Number of workers in dask distributed client.")

    args = parser.parse_args()
    main(args)
    

