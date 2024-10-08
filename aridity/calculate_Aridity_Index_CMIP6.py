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
import pickle
import argparse
import os
import sys

sys.path.append('/g/data/mn51/users/jb6465/drought-github/submodules/gwls')
from gwl import get_GWL_timeslice
#sys.path.append('/g/data/mn51/users/jb6465/drought-github/percentiles_spi')
#import utils

# Ignore warnings
import warnings
import logging
warnings.filterwarnings('ignore') 
logging.getLogger("distributed.worker.memory").setLevel(logging.ERROR)

#################< Functions ########################
def calc_AI(pr_file,e0_file,syear,eyear):
    # Resample to annually
    da_e0 = e0_file.resample(time='YE').sum('time')
    da_pr = (pr_file*86400).resample(time='YE').sum('time')

    # Select GWL periods
    e0 = da_e0.sel(time=slice(syear,eyear)).chunk({'time':-1,'lat':'auto','lon':'auto'})
    pr = da_pr.sel(time=slice(syear,eyear)).chunk({'time':-1,'lat':'auto','lon':'auto'})

    #Compute AI
    computed_AI = (pr/e0).astype("float32").rename("AI")
        
    # add attributes
    computed_AI.attrs['units']         = ''
    computed_AI.attrs['long_name']     = 'Aridity index, ratio of precipitation to (potential)evepotranspiration'
    computed_AI.attrs['standard_name'] = 'Aridity Index'
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

    RCM = inargs.RCM
    #MPI-ESM1-2-HR // NorESM2-MM for BARPA only and CNRM-ESM2-1 for CCAM only
    model_list = ['CMCC-ESM2', 'ACCESS-ESM1-5', 'ACCESS-CM2', 'EC-Earth3', 'CESM2'] 
    model_list = model_list + ['CNRM-ESM2-1'] if RCM == 'CCAM-v2203-SN' else model_list + ['MPI-ESM1-2-HR', 'NorESM2-MM']

    #< Variable for (P)ET component
    if inargs.index == 'atmospheric-based':
        var_e = "var_pet"
    elif inargs.index == 'plant-based':
        var_e = "var_et"




    for model in model_list:
        print('========= '+RCM+'_'+model+' =========')
        bc_string = '_ACS-{}-{}-{}-{}.nc'.format(inargs.bcMethod, inargs.bcSource, '1960' if inargs.bcSource == 'AGCD' else '1979', '2022') if inargs.bc == 'output' else '.nc'
        variant_id = lib_david.data_source['CMIP6'][model]['variant-id']
        file_name_ann = "{}/AI_{}_{}_{}_{}_{}_{}_{}{}".format(inargs.outputDir, inargs.index,'AGCD-05i',model,'ssp370',variant_id,'BOM' if RCM == 'BARPA-R' else 'CSIRO','v1-r1','_raw.nc' if inargs.bc == 'raw' else bc_string)

        ###############################  Compute annual time series ############################################

        if os.path.exists(file_name_ann)==False:
            print("Computing annual time series. File: {name}...".format(name=file_name_ann))
        
            # read input data for AI calculation
            input_array_p = lib_david.load_target_variable('var_p', RCM, model, bc=inargs.bc, bc_method=inargs.bcMethod, bc_source=inargs.bcSource)
            input_array_e = lib_david.load_target_variable(var_e, RCM, model, bc=inargs.bc, bc_method=inargs.bcMethod, bc_source=inargs.bcSource)

            print("================ input_array_p =============")
            print(input_array_p)
            print("================ input_array_e =============")
            print(input_array_e)

        
            ds_AI_ann = calc_AI(input_array_p,input_array_e,'1960','2100')
    
            ds_AI_ann.attrs['description'] = f"Aridity Index defined as the ratio of precipitation to (potential) evapotranspiration. Produced for ACS."
            ds_AI_ann.attrs['method']  = 'Using  {} aridity: pr/{}'.format(inargs.index,var_e)
            ds_AI_ann.attrs['history'] = cmdprov.new_log(extra_notes=[get_git_hash()])
            ds_AI_ann.attrs['bias correction'] = f"None. Raw {RCM} output data." if inargs.bc == 'raw' else f"{RCM} {bc_string} data";
            print("Computing {name}...".format(name=file_name_ann))
                        
            #< Save output
            saver_ann = ds_AI_ann.to_netcdf(file_name_ann,compute=False)
            future_ann = client.persist(saver_ann)
            dask.distributed.progress(future_ann)
            future_ann.compute()
        else:
            print("{name} exists. Pass.".format(name=file_name_ann))
                        
    #             for gwl in models_gwl[model][rcp]:
    #                 print(gwl)
    #                 syear = str(models_gwl[model][rcp][gwl][0])
    #                 eyear = str(models_gwl[model][rcp][gwl][1])
    #                 print(f"Start year for {model}, {rcp}, {gwl}: {syear}")
    #                 print(f"End year for {model}, {rcp}, {gwl}: {eyear}")

    #                 ###############################  Truncate to annual GWL time series IF previous file exists #############################
    #                 file_name_ann_gwl = "{}AI-{}_NHP1-AUS-5_{}_{}_{}_{}_{}_{}.nc".format(inargs.OutputDir,inargs.index,model,rcp,run,bc,'annual','GWL'+str(int(float(gwl)*10))) 

    #                 if os.path.exists(file_name_ann_gwl)==False:

    #                     ds_AI_ann_gwl = xr.open_dataset(file_name_ann).sel(time=slice(syear,eyear))

    #                     ds_AI_ann_gwl.attrs['description'] = f'Ratio of precipitation to (potential)evepotranspiration produced from National Hydrological Projections (NHP1.0) on /g/data/wj02/COMPLIANT_PUBLISHED/. Produced for ACS. '
    #                     ds_AI_ann_gwl.attrs['method']  = 'Using  {} aridity: pr/{}'.format(inargs.index,var_e)
    #                     ds_AI_ann_gwl.attrs['history'] = cmdprov.new_log(extra_notes=[get_git_hash()])
    #                     ds_AI_ann_gwl.attrs['comment'] = "Using data on /g/data/wj02/COMPLIANT_PUBLISHED/" ;
    #                     ds_AI_ann_gwl.attrs['cell_methods'] = "time: mean" ;
    #                     print("Computing {name}...".format(name=file_name_ann_gwl))

    #                     #< Save output
    #                     ds_AI_ann_gwl.to_netcdf(file_name_ann_gwl)

    #                 else:
    #                     print("{name} exists. Pass.".format(name=file_name_ann_gwl))

    #                 ###############################  Create 2D GWL IF time series file exists ############################################
    #                 file_name_2D_gwl = "{}AI-{}_NHP1-AUS-5_{}_{}_{}_{}_{}_{}.nc".format(inargs.OutputDir,inargs.index,model,rcp,run,bc,'2D','GWL'+str(int(float(gwl)*10)))

    #                 if os.path.exists(file_name_2D_gwl)==False:
                            
    #                     ds_AI_2D_gwl = xr.open_dataset(file_name_ann).sel(time=slice(syear,eyear)).mean('time')

    #                     ds_AI_2D_gwl.attrs['description'] = f'Ratio of precipitation to (potential)evepotranspiration produced from National Hydrological Projections (NHP1.0) on /g/data/wj02/COMPLIANT_PUBLISHED/. Produced for ACS. '
    #                     ds_AI_2D_gwl.attrs['method']  = 'Using  {} aridity: pr/{}'.format(inargs.index,var_e)
    #                     ds_AI_2D_gwl.attrs['history'] = cmdprov.new_log(extra_notes=[get_git_hash()])
    #                     ds_AI_2D_gwl.attrs['comment'] = "Using data on /g/data/wj02/COMPLIANT_PUBLISHED/" ;
    #                     ds_AI_2D_gwl.attrs['cell_methods'] = "time: mean" ;
    #                     print("Computing {name}...".format(name=file_name_2D_gwl))

    #                     #< Save output
    #                     ds_AI_2D_gwl.to_netcdf(file_name_2D_gwl)
    #                 else:
    #                     print("{name} exists. Pass.".format(name=file_name_2D_gwl))




    #< Close the client
    client.close()


if __name__ == '__main__':
    extra_info =""" 
author:
    David Hoffmann, david.hoffmann@bom.gov.au
    Created: 29/05/2024
    Modified: 17/06/2024
"""
    description = """
    Calculate the aridity index from gridded National Hydrological Projections (NHP1.0) data from '/g/data/wj02/COMPLIANT_PUBLISHED/'.     
        """    
    parser = argparse.ArgumentParser(description=description,
                                     epilog=extra_info,
                                     argument_default=argparse.SUPPRESS,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--RCM", type=str, choices=['BARPA-R', 'CCAM-v2203-SN'], help="Choose RCM from ['BARPA-R', 'CCAM-v2203-SN']")
    parser.add_argument("--bc", type=str, choices=['raw','input','output'], help="Choose either 'raw' (py18/hq89 raw BARPA/CCAM), 'input' (ia39 bc input, 5km) or 'output' (ia39 bc output, 5km).")
    parser.add_argument("--bcSource", type=str, default='AGCD', choices=['AGCD', 'BARRA'], help="Choose either 'AGCD', 'BARRA'. Default is 'AGCD'")
    parser.add_argument("--bcMethod", type=str, default='QME', choices=['QME', 'MRNBE'], help="Choose either 'MRNBC', 'QME'. Default is 'QME'")
    parser.add_argument("--index", type=str, choices=['atmospheric-based','plant-based'], help="Choose either 'atmospheric-based' or 'plant-based' aridity index.")
    parser.add_argument("--outputDir", type=str, default='/g/data/ia39/ncra/drought_aridity/ai/acs_downscaled_notBC_5km/', help="Output directory on Gadi. Default is'/g/data/ia39/ncra/drought_aridity/ai/'")
    parser.add_argument("--nworkers", type=int, default=15, help="Number of workers in dask distributed client.")

    args = parser.parse_args()
    main(args)
    

