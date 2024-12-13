'''
######################################################################
## BUREAU OF METEOROLOGY
## AUSTRALIAN CLIMATE SERVICE
## DROUGHT HAZARD TEAM
##
## DATE:             Jun-2024
## SCRIPT:           spi.py
## AUTHOR:           jessica.bhardwaj@bom.gov.au
##
## DESCRIPTION:      Script to compute SPI for ACS work package 3 deliverables
##
#######################################################################
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

# Import utils
sys.path.append('/g/data/mn51/users/jb6465/drought-github')
import utils

# Ignore warnings
import warnings
import logging
warnings.filterwarnings('ignore') 
logging.getLogger("distributed.worker.memory").setLevel(logging.ERROR)

def compute_spi(full_period, base_period_start_year, base_period_end_year, dataset_source):

    import numpy as np
    import xarray as xr
    from scipy.stats import gamma, norm
    
    ### Gamma distribution requires non-zero values - assign random value to avoid NaN gamma fit errors
    base_period = full_period.sel(time=full_period['time'].dt.year.isin(range(base_period_start_year, base_period_end_year+1)))
    
    
    if dataset_source == 'AGCD': #deal with AGCD offshore lat/lon that have the same value assigned to them at all time slices (gives NaN errors when fitting). Masking doesn't work, just assign random non zero value.
        base_period = base_period.where(~da.all(base_period == base_period.isel(time=0), axis=0), np.random.uniform(0.1, 1, size=base_period.shape))
    non_zero_base_period = base_period.where(base_period > 0, np.random.uniform(0.1, 1, size=base_period.shape))

    ### Determine fit parameters along the time axis
    alpha, loc, beta = xr.apply_ufunc(
        gamma.fit,
        non_zero_base_period,
        input_core_dims=[['time']],
        output_core_dims=[[], [], []],
        vectorize=True,
        dask='parallelized',
        dask_gufunc_kwargs={'allow_rechunk': True},
        output_dtypes=[float, float, float],
        kwargs={'floc': 0, 'method': 'MLE'}
    )
    
    ### Apply gamma parameters from base period to gamma CDFs for full period accounting for zero observations
    p_zero = ((full_period == 0).sum(dim='time') / full_period.sizes['time']).values
    gamma_cdf = gamma.cdf(full_period, alpha, scale=beta)
    cdf = p_zero + (1 - p_zero) * gamma_cdf
    
    ### Apply inverse normal distribution
    norm_ppf = norm.ppf(cdf)

    return xr.DataArray(norm_ppf, coords=full_period.coords, dims=full_period.dims)


def get_git_hash():
    """Returns the git hash for the working repository"""
    git_repo = git.Repo(sys.argv[0], search_parent_directories=True)
    git_root = git_repo.git.rev_parse("--show-toplevel")
    git_hash = git.Repo(git_root).heads[0].commit
    git_text = " (Git hash: %s)" %(str(git_hash)[0:7])
    return git_text

################ Main ###############################

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
    client = Client(n_workers=inargs.nworkers, threads_per_worker=1, local_directory = tempfile.mkdtemp(), memory_limit = "63000mb", dashboard_address=":8787") #need to make this dynamic depending on cpus/memory requested
    print("Dask dashboard is available at:", client.dashboard_link)

    RCM = inargs.RCM
        
    # Compute SPI for input GCM/RCM/bc file
    GCM = inargs.GCM
    print('========= '+RCM+'_'+GCM+' =========')
    bc_string = '_ACS-{}-{}-{}-{}.nc'.format(inargs.bcMethod, inargs.bcSource, '1960' if inargs.bcSource == 'AGCD' else '1979', '2022') if inargs.bc == 'output' else '.nc'
    variant_id = utils.data_source['CMIP6'][GCM]['variant-id']
    file_name = "{}/SPI{}_{}_{}_{}_{}_{}_{}_{}{}".format(inargs.outputDir, inargs.spiAccumulation,'AGCD-05i',GCM,'ssp370',variant_id,'BOM' if 'BARPA' in RCM else 'CSIRO' if 'CCAM' in RCM else RCM if 'NARCliM' in RCM else 'UQ-DES','v1-r1','baseperiod'+(str(inargs.basePeriodStart)+str(inargs.basePeriodEnd)),'_raw.nc' if inargs.bc == 'raw' else bc_string) if inargs.datasetSource == 'CMIP6' else "{}/SPI{}_{}_{}.nc".format(inargs.outputDir, inargs.spiAccumulation,'AGCD','baseperiod'+(str(inargs.basePeriodStart)+str(inargs.basePeriodEnd)))
    
    if os.path.exists(file_name)==False:
        print("Computing {name}...".format(name=file_name))
        # Load in files as dask xarray
        input_array = utils.load_target_variable(inargs.datasetSource, 'var_p', RCM, GCM, inargs.spiAccumulation, bc=inargs.bc, bc_method=inargs.bcMethod, bc_source=inargs.bcSource)
        input_array = input_array.astype(np.float32).chunk({'time':-1, 'lat':'auto', 'lon':'auto'})
        # Group by month and apply the SPI calculation
        spi_array_grouped = input_array.groupby('time.month').map(lambda x: compute_spi(x, inargs.basePeriodStart, inargs.basePeriodEnd, inargs.datasetSource))
        
        ds_SPI = spi_array_grouped.assign_coords(time=input_array['time'])
        ds_SPI = ds_SPI.rename('SPI{}'.format(inargs.spiAccumulation))

        ds_SPI.attrs['description'] = "Standardised Precipitation Index computed using method of McKee et al. 1993 using a base period of {}-{}. Further details in supporting technical documentation.".format(inargs.basePeriodStart, inargs.basePeriodEnd)
        ds_SPI.attrs['created'] = (datetime.now()).strftime("%d/%m/%Y %H:%M:%S")    
        ds_SPI.attrs['history'] = cmdprov.new_log(extra_notes=[get_git_hash()])
        if inargs.datasetSource == 'CMIP6':
            ds_SPI.attrs['driving GCM'] = GCM
            ds_SPI.attrs['variant ID'] = variant_id
            ds_SPI.attrs['RCM'] = RCM, utils.data_source['CMIP6'][GCM][RCM] if RCM=='UQ-DES' else RCM

        # Specify compression options to avoid large file sizes
        encoding = {'SPI{}'.format(inargs.spiAccumulation): {'zlib': True, 'complevel': 5, 'dtype':'float32'}}
        
        # Save output
        ds_SPI.to_netcdf(file_name, encoding=encoding)
        
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
    Calculate the Standardised Precipitation Index for ACS work package 3 deliverables.     
        """    
    parser = argparse.ArgumentParser(description=description,
                                     epilog=extra_info,
                                     argument_default=argparse.SUPPRESS,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
                                     
    parser.add_argument("--datasetSource", type=str, choices=['AGCD', 'CMIP6'], help="Choose dataset source from ['AGCD', 'CMIP6']") #soon to add ERA5/ERA5-Land
    parser.add_argument("--RCM", type=str, choices=['BARPA-R', 'CCAM-v2203-SN', 'NARCliM2-0-WRF412R3', 'NARCliM2-0-WRF412R5', 'UQ-DES'], help="Choose RCM from ['BARPA-R', 'CCAM-v2203-SN', 'NARCliM2-0-WRF412R3', 'NARCliM2-0-WRF412R5', 'UQ-DES']")
    parser.add_argument("--GCM", type=str, choices=['ACCESS-CM2', 'ACCESS-ESM1-5', 'CMCC-ESM2', 'CESM2', 'CNRM-CM6-1-HR', 'CNRM-ESM2-1', 'EC-Earth3', 'EC-Earth3-Veg', 'FGOALS-g3', 'GFDL-ESM4', 'GISS-E2-1-G',  'MPI-ESM1-2-HR', 'MPI-ESM1-2-LR', 'MRI-ESM2-0', 'NorESM2-MM', 'UKESM1-0-LL'], help="Input GCM - check ia39 for correct GCM/RCM availabilies")
    parser.add_argument("--bc", type=str, choices=['raw','input','output'], help="Choose either 'raw' (py18/hq89 raw BARPA/CCAM), 'input' (ia39 bc input, 5km) or 'output' (ia39 bc output, 5km).")
    parser.add_argument("--bcSource", type=str, default='AGCD', choices=['AGCD', 'BARRA'], help="Choose either 'AGCD', 'BARRA'. Default is 'AGCD'")
    parser.add_argument("--bcMethod", type=str, default='QME', choices=['QME', 'MRNBC'], help="Choose either 'MRNBC', 'QME'. Default is 'QME'")
    parser.add_argument("--spiAccumulation", type=int, default=3, help="Choose SPI accumulation i.e. SPI-1, SPI-3, SPI-6, SPI-12. Default is 3")
    parser.add_argument("--basePeriodStart", type=int, default=1965, help="Specify start year for base period. Minimum 30 years advised. 50 years helps avoid NaN errors on upper tails. Default = 1965.")
    parser.add_argument("--basePeriodEnd", type=int, default=2014, help="Specify end year for base period. Minimum 30 years advised. 50 years helps avoid NaN errors on upper tails. Default = 2014")
    parser.add_argument("--outputDir", type=str, default='/g/data/ia39/ncra/drought_aridity/spi', help="Output directory on Gadi. Default is'/g/data/ia39/ncra/drought_aridity/spi'")
    parser.add_argument("--nworkers", type=int, default=10, help="Number of workers in dask distributed client.")

    args = parser.parse_args()
    main(args)
    

