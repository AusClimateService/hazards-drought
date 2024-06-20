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
from datetime import datetime
import argparse
import os
import sys

# Ignore warnings
import warnings
warnings.filterwarnings('ignore') 

#################< Functions ########################
def calc_AI(pr_file,e0_file,syear,eyear):
    # Open data since it's a single file and resample annually
    da_e0 = xr.open_mfdataset(e0_file).resample(time='YE').sum('time')
    da_pr  = xr.open_mfdataset(pr_file).resample(time='YE').sum('time')

    # Select GWL periods
    e0 = da_e0.sel(time=slice(syear,eyear))['e0'].chunk({'time':-1,'lat':'auto','lon':'auto'})
    pr = da_pr.sel(time=slice(syear,eyear))['pr'].chunk({'time':-1,'lat':'auto','lon':'auto'})

    #Compute AI
    computed_AI = (pr/e0).mean('time').astype("float32").rename("AI")
        
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
    
    dask.config.set({
        "distributed.scheduler.worker-saturation": 1, #< This should use the new behaviour which helps with memory pile up
        })
    client = Client(n_workers=inargs.nworkers, threads_per_worker=1, local_directory = tempfile.mkdtemp()) #need to make this dynamic depending on cpus/memory requested

    #< Keep attributes from input datasets
    xr.set_options(keep_attrs = True)

    # Specify the path to .pkl file
    file_path = '/g/data/mn51/users/dh4185/hazards-drought/gwl_years.pkl'
    
    # Open the file in binary read mode
    with open(file_path, 'rb') as file:
        # Load the data from the file
        models_gwl = pickle.load(file)
    
    print("From pkl fike:",models_gwl)
    #< need to add 10 years to GWLs to calculate AI over 30 years?

    if inargs.index == 'atmospheric-based':
        var_e = "e0"
    elif inargs.index == 'plant-based':
        var_e = "etot"

    path_templ_pr = "/g/data/wj02/COMPLIANT_PUBLISHED/HMINPUT/output/AUS-5/BoM/"#{}/{}/{}/{}/latest/day/pr/"
    path_templ_e0 = "/g/data/wj02/COMPLIANT_PUBLISHED/HMOUTPUT/output/AUS-5/BoM/"#{}/{}/{}/{}/latest/day/pr/"
    
    files_e0 = lib_david.get_file_paths(path_templ_e0,".nc",include=["rcp45",var_e]) + lib_david.get_file_paths(path_templ_e0,".nc",include=["rcp85",var_e])
    files_pr = lib_david.get_file_paths(path_templ_pr,".nc",include=["rcp45","pr"],exclude=["BEFORE"]) + lib_david.get_file_paths(path_templ_pr,".nc",include=["rcp85","pr"],exclude=["BEFORE"])
    
    bc_method = ["r240x120-QME","CSIRO-CCAM-r3355-r240x120-ISIMIP2b","_r240x120-ISIMIP2b","r240x120-MRNBC"]

    
    for model in inargs.GCM:
        print('========= '+model+' =========')
        
        for rcp in models_gwl[model]:
            files_e0 = lib_david.get_file_paths(path_templ_e0,".nc",include=[rcp,model,"e0"])
            files_pr = lib_david.get_file_paths(path_templ_pr,".nc",include=[rcp,model,"pr"],exclude=["BEFORE"])
            print(rcp)
    
            for bc in bc_method:
                infile_e0 = [filename_e0 for filename_e0 in files_e0 if bc in filename_e0]
                infile_pr = [filename_pr for filename_pr in files_pr if bc in filename_pr]
                print(bc)
    
                for gwl in models_gwl[model][rcp]:
                    print(gwl)
                    syear = str(models_gwl[model][rcp][gwl][0])
                    eyear = str(models_gwl[model][rcp][gwl][1])
    
                    file_name = "/scratch/mn51/dh4185/AI-{}_NHP10_{}_{}_{}_{}.nc".format(inargs.index,model,bc,rcp,gwl)
    
                    if os.path.exists(file_name)==False:
    
                        now = datetime.now() # get date and time
    
                        print("Computing {name}...".format(name=file_name))
                        ds_AI = calc_AI(infile_pr,infile_e0,syear,eyear)
    
                        ds_AI.attrs['description'] = f'Ratio of precipitation to (potential)evepotranspiration produced from National Hydrological Projections (NHP1.0) on /g/data/wj02/COMPLIANT_PUBLISHED/. Produced for ACS. '
                        ds_AI.attrs['method']  = 'Using  {} aridity: pr/{}'.format(inargs.index,var_e)
                        ds_AI.attrs['created'] = now.strftime("%d/%m/%Y %H:%M:%S")    
                        ds_AI.attrs['history'] = cmdprov.new_log(extra_notes=[get_git_hash()])
                        ds_AI.attrs['comment'] = "Using data on /g/data/wj02/COMPLIANT_PUBLISHED/" ;
                        ds_AI.attrs['cell_methods'] = "time: mean" ;
                        
                         #< Save output
                        saver = ds_AI.to_netcdf(file_name,compute=False)
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
    David Hoffmann, david.hoffmann@bom.gov.au
    Created: 29/05/2024
    Modified: --
"""
    description = """
    Calculate the aridity index from gridded National Hydrological Projections (NHP1.0) data from '/g/data/wj02/COMPLIANT_PUBLISHED/'.     
        """    
    parser = argparse.ArgumentParser(description=description,
                                     epilog=extra_info,
                                     argument_default=argparse.SUPPRESS,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
                                     
    parser.add_argument("--index", type=str, choices=['atmospheric-based','plant-based'], help="Choose either 'atmospheric-based' or 'plant-based' aridity index.")
    parser.add_argument("--GCM", nargs='+', type=str, choices=['ACCESS1-0','GFDL-ESM2M','MIROC5','CNRM-CM5'], default = ['ACCESS1-0','GFDL-ESM2M','MIROC5','CNRM-CM5'], help="Choose GCM: 'ACCESS1-0','GFDL-ESM2M','MIROC5','CNRM-CM5'. If none are selected, the default is all GCMs")
    parser.add_argument("--OutputDir", type=str, default='/g/data/ia39/ncra/drought_aridity/ai/', help="Output directory on Gadi. Default is'/g/data/mn51/projects/work_package_4/climate_hazard_indices/drought/'")
    parser.add_argument("--nworkers", type=int, default=15, help="Number of workers in dask distributed client.")

    args = parser.parse_args()
    main(args)
    

