#< frequently used functions
import os
import xarray as xr
import dask.distributed
import glob
from dask.distributed import client
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def get_file_paths(root_directory, file_extension,include=None,exclude=None):
    '''
    Collects files names recursively and uses file extension and list of custom terms to narrow search.
    In:
        root_directory: string of directory pointing to the desired files
        file_extension: Usually .nc
        include: list of strings to include in the search
        exclude: list of strings to exclude in the search

    Out:
        file_paths: list of file paths including files names.
    '''
    file_paths = []  # List to store file paths

    # Traverse through the directory structure recursively
    for root, directories, files in os.walk(root_directory):
        for filename in files:
            # Check if the file has the desired file extension
            if filename.endswith(file_extension):
                # Construct the full path to the file by joining the directory path and file name
                file_path = os.path.join(root, filename)
                if include:
                    if all(term in filename for term in include):
                        if exclude:
                            if all (term not in filename for term in exclude):
                                file_paths.append(file_path)
                        else:
                            file_paths.append(file_path)
                else:
                    file_path = os.path.join(root, filename)
                    file_paths.append(file_path)

    return file_paths


def write_netcdf(da,file_name):
    """  Write data to netcdf if it doesn't exist yet. Otherwise pass.
    Args:
        da - data to be written
        file_name - file name including output directory
    """
    if os.path.exists(file_name)==False:
        saver = da.to_netcdf(file_name,compute=False)
        future = client.persist(saver)
        dask.distributed.progress(future)
        future.compute()
        print('\n')
    else:
        print("{name} exists. Pass.".format(name=file_name))


def get_agcd_q_mask(ClimStartYr,ClimEndYr,mask_file,ref_ds):
    """ Creates a mask by having a weight of at least 1 (minimal influence to influence of nearby stations on climatology for grid cell) for at least 80% of days in the climatological period from ClimStartYr to ClimEndYr.
    Args:
        ClimStartYr
        ClimEndYr
    Returns:
        numpy masked array
    """
    if os.path.exists(mask_file):
        ds = xr.open_dataset(mask_file)
    else:
        da_weights            = xr.open_mfdataset('/g/data/zv2/agcd/v1-0-1/precip/weight/r005/01day/agcd_v1-0-1_precip_weight_r005_daily_*.nc').sel(time=slice(ClimStartYr,ClimEndYr))
        da_weights            = da_weights.sel(lon=slice(int(ref_ds.longitude[0]),int(ref_ds.longitude[-1]))).astype("float32")
        da_weights            = da_weights.sortby(da_weights["lat"], ascending=False).sel(lat=slice(int(ref_ds.latitude[0]),int(ref_ds.latitude[-1]))).astype("float32").chunk(dict(time="auto"))
        da_weights_count      = (da_weights['weight'] >= 1.0).sum('time', keep_attrs=True)
        ds                    = (da_weights_count / len(da_weights['time'])).to_dataset(name='fraction')
        ds.attrs["comment"]   = f"Fraction of days where observations influence the data. It counts the fraction of days where the weight value is >=1 in the period {ClimStartYr}-{ClimEndYr}"
        ds.attrs['long_name'] = mask_file
        # ds.attrs['history']   = cmdprov.new_log(extra_notes=[get_git_hash()])
        ds.attrs['agcd_version']= da_weights.attrs['agcd_version']
        write_netcdf(ds,mask_file)
        
    q_mask = np.ma.masked_greater(ds.fraction,0.8)
    
    return q_mask

def get_git_hash(): # not working in Jupyter notebook but add to exe python script
    """Returns the git hash for the working repository"""
    git_repo = git.Repo(sys.argv[0], search_parent_directories=True)
    git_root = git_repo.git.rev_parse("--show-toplevel")
    git_hash = git.Repo(git_root).heads[0].commit
    git_text = " (Git hash: %s)" %(str(git_hash)[0:7])
    return git_text


def generate_diverging_colors(colormap_name, num_intervals):
    """
    Generate hex color codes for a diverging color bar with specified intervals.

    Parameters:
    - colormap_name: Name of the colormap to use.
    - num_intervals: Number of intervals for the color bar.

    Returns:
    - List of hex color codes.
    """
    # Generate color map
    cmap = plt.cm.get_cmap(colormap_name, num_intervals)

    # Define colors for intervals
    colors = [cmap(i) for i in range(num_intervals)]

    # Function to convert RGB tuples to hex codes
    def rgb_to_hex(rgb):
        return '#{:02x}{:02x}{:02x}'.format(int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))

    # Convert colors to hex codes
    hex_codes = [rgb_to_hex(color) for color in colors]

    return hex_codes


def plot_AI(data_dict,method,plot_title,plot_file=False,):

    # Check if the method is one of the allowed values
    allowed_methods = ['absolute', 'categories', 'change_prct', 'change_simple']

    if method not in allowed_methods:
        raise ValueError("Invalid method. Please choose from 'absolute', 'change_simple' or 'change_prct'.")

        
    if method == 'absolute':
        pass
        
    elif method == 'categories':
        # Define color categories and corresponding colors
        categories = ['Hyper-arid', 'Arid', 'Semi-arid', 'Sub-humid', 'Humid']
        colors = ['#e2211c','#f58e3d', '#f9cc5d', '#fffdb2', '#666666']  
    
        # Define category thresholds
        thresholds = [0, 0.05, 0.2, 0.5, 0.65, 1]
        _ticks = [0.025, 0.125, 0.35, 0.575, 0.825]

    elif method == 'change_prct':
        # Define color categories and corresponding colors
        categories = ['<-20%', '<-15-20%','<-10-15%','<-5-10%','<-2-5%','No change\n(-2-+2%)','<2-5%','>5-10%', '>10-15%','>15-20%','>20%']
        colors = generate_diverging_colors('BrBG',11)#['#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#003c30'] # ['#b9141d','#ce4629','#de7149','#ef9c6d','#f8cb92','#f8f6c0','#d3df9b','#aacd7b','#81ba5d','#5aa73d','#37922d']
    
        # Define category thresholds
        thresholds = [-100,-20,-15,-10,-5,-2,2,5,10,15,20,100]
        _ticks = [-60,-17.5,-12.5,-7.5,-3.5,0,3.5,7.5,12.5,17.5,60]

    elif method == 'change_simple':
        # Define color categories and corresponding colors
        categories = ['drying', 'wetting']
        colors = ['#D2B48C','#87CEEB']  # Using hex color code for light brown
    
        # Define category thresholds
        thresholds = [-5, 0, 5]
        _ticks = [-2.5,2.5]
        
    # Create a custom colormap
    cmap = plt.cm.colors.ListedColormap(colors, name='custom_cmap', N=len(colors))
    
    # Create subplots for each model
    num_models = len(data_dict)
    num_rows = (num_models + 1) // 2
    if len(data_dict) == 1:
        num_cols = 1
        plt_width = 7
    else: 
        num_cols = 3
        plt_width = 4
        
    # Create subplots with two columns if more than two files to plot
    fig, axs = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(10, plt_width * num_rows), subplot_kw={'projection': ccrs.PlateCarree()}, sharex=True, sharey=True)
    
    # Flatten the 2D array of subplots for easy indexing
    axs = np.ravel([axs]) 
    
    # Iterate over each model in the dictionary and plot the data
    for i, (model_name, da) in enumerate(data_dict.items()):
        # Create a normalization instance for mapping the data values to the colormap
        norm = plt.cm.colors.BoundaryNorm(thresholds, cmap.N, clip=True)
    
        # Plot the data with the custom colormap and color categories
        if 'AGCD'in model_name:
            plot = da.sel(lat=slice(-5,-45),lon=slice(110,155)).plot.pcolormesh(ax=axs[i], cmap=cmap, norm=norm, add_colorbar=False)
        else:
            plot = da.sel(lat=slice(-45,-5),lon=slice(110,155)).plot.pcolormesh(ax=axs[i], cmap=cmap, norm=norm, add_colorbar=False)
      
        # Set subplot title and labels
        axs[i].set_title(model_name) if len(data_dict) > 1 else axs[i].set_title(plot_title)
        axs[i].set_xlabel('Longitude')
        axs[i].set_ylabel('Latitude')

        # Add coastlines using cartopy
        extent = [110, 155, -45, -5]
        axs[i].add_feature(cfeature.COASTLINE)
        axs[i].set_extent(extent)
        # Add gridlines and ticks
        # axs[i].gridlines(draw_labels=True)
        axs[i].set_xticks(np.arange(extent[0], extent[1] + 1, 5), crs=ccrs.PlateCarree())
        axs[i].set_yticks(np.arange(extent[2], extent[3] + 1, 5), crs=ccrs.PlateCarree())

    
    # Set the main title of the plot if multi_plot
    if len(data_dict) > 1: plt.suptitle(plot_title, fontsize=16)
    
    # Adjust the layout to make room for the colorbar at the bottom
    plt.subplots_adjust(bottom=0.2, right=0.8, top=0.9, hspace=0.3)
    
    # Add a colorbar with category labels to the last subplot
    if len(data_dict) > 1:
        cbar_ax = fig.add_axes([0.18, 0.02, 0.7, 0.01])  # Adjust position and size as needed
        cbar = plt.colorbar(plot, cax=cbar_ax, pad=0.1, ticks=_ticks, boundaries=thresholds, orientation='horizontal')
    else:
        # Add a colorbar with category labels
        cbar = plt.colorbar(plot, ax=axs, ticks=_ticks, boundaries=thresholds, orientation='vertical')

    cbar.set_ticklabels(categories)
    if method == 'absolute': cbar.set_label('Categories')
    
    # Hide any remaining empty subplots
    for i in range(len(data_dict), len(axs)):
        axs[i].axis('off')
    
    # Ensure the main title is not cut off
    if len(data_dict) > 1: plt.tight_layout(rect=[0, 0.03, 1, 1])

    # Save the plot as a PNG file in the data directory
    if plot_file: plt.savefig(plot_file)

    plt.show()


def mean_bias_corr (data, ref, group_by='time.month'):
    #< Format and crop lat and lon to be identical
    ref = ref.rename({"latitude": "lat"})
    ref = ref.rename({"longitude": "lon"})
    ref = ref.sel(lon=slice(int(data.lon[0]),int(data.lon[-1]))).astype("float32")
    ref = ref.sortby(ref["lat"], ascending=False).sel(lat=slice(int(data.lat[0]),int(data.lat[-1]))).astype("float32")
    
    if group_by:
        bias_corr = data - (data.groupby(group_by).mean('time') - ref.groupby(group_by).mean('time'))
    else:
        bias_corr = data - (data.mean('time') - ref.mean('time'))
    return bias_corr
    
def get_files(source,source_char,model=None):
    if source == 'AGCD':
        #< AWRA-L and AGCD directory (AWRA-L only until 2017 on fj8. Using AWO product from iu04)
        files = {'pet_files':sorted(get_file_paths("/g/data/iu04/australian-water-outlook/historical/v1/AWRALv7/",".nc",["e0"])), # What's the difference to 'public'?
                 'p_files':sorted(get_file_paths("/g/data/zv2/agcd/v2-0-1/precip/total/r005/01month/",".nc",["agcd_v2-0-1_precip_total_r005_monthly"]))}

    elif source == 'ERA5-land':
        #< ERA5-Land
        files = {'pet_files':sorted(get_file_paths("/g/data/zz93/era5-land/reanalysis/pev/","nc",["pev_era5-"])), 
                 'p_files':sorted(get_file_paths("/g/data/zz93/era5-land/reanalysis/tp/","nc",["tp_era5-"]))}

    elif source == 'ERA5':
        #< ERA5
        era5_path = "/g/data/ia39/australian-climate-service/release/CORDEX-CMIP6/output/AUS-15/BOM/ECMWF-ERA5/evaluation/r1i1p1f1/BOM-BARPA-R/v1/mon/{}/"
        
        files = {'pet_files':sorted(get_file_paths(era5_path.format('evspsblpot'),'.nc')),
                 'p_files':sorted(get_file_paths(era5_path.format('pr'),'.nc'))}
                
    elif source == 'CMIP6':
            
        print(model)

        # cmip6_dir = "/g/data/py18/BARPA/output/CMIP6/DD/AUS-15/BOM/{}/{}/{}/BARPA-R/v1-r1/mon/pr/v20231001/"   
        cmip6_dir = "/g/data/ia39/australian-climate-service/test-data/CORDEX-CMIP6/bias-adjustment-output/AGCD-05i/BOM/{}/{}/{}/BARPA-R/v1-r1/mon/pr/v20231001/"
        cmip6_hist = cmip6_dir.format(source_char[source][model]['name'],'historical',source_char[source][model]['variant-id'])
        cmip6_ssp370 = cmip6_dir.format(source_char[source][model]['name'],'ssp370',source_char[source][model]['variant-id'])
        # cmip6_ssp126 = cmip6_dir.format(source_char[source][model]['name'],'ssp126',source_char[source][model]['variant-id']) # how to add?
            
        p_files = sorted(get_file_paths(cmip6_hist,'.nc'))+sorted(get_file_paths(cmip6_ssp370,'.nc'))
 
        cmip6_dir = "/g/data/py18/BARPA/output/CMIP6/DD/AUS-15/BOM/{}/{}/{}/BARPA-R/v1-r1/mon/evspsblpot/v20231001/"   
        cmip6_hist = cmip6_dir.format(source_char[source][model]['name'],'historical',source_char[source][model]['variant-id'])
        cmip6_ssp370 = cmip6_dir.format(source_char[source][model]['name'],'ssp370',source_char[source][model]['variant-id'])
        # cmip6_ssp126 = cmip6_dir.format(source_char[source][model]['name'],'ssp126',source_char[source][model]['variant-id'])
            
        pet_files = sorted(get_file_paths(cmip6_hist,'.nc'))+sorted(get_file_paths(cmip6_ssp370,'.nc'))

        files = {'pet_files': pet_files,
                 'p_files': p_files}
    return files

def calc_AI(files,source,var_e0,var_p,model=None):
    
    #< Open datasets
    print("=== opened datasets ===")
    ds_e0 = xr.open_mfdataset(files['e0'],combine='nested', concat_dim='time',parallel=True)
    ds_p  = xr.open_mfdataset(files['p'],combine='nested', concat_dim='time',parallel=True)
    try:
        ds_e0 = ds_e0.chunk(dict(time="auto"))#,latitude=128,longitude=128)) #< Open daily PET (e0) data 
        ds_p = ds_p.chunk(dict(time="auto")) #< Open monthly rainfall (p) data 
    except NotImplementedError:
        print("Continue without auto-rechunking")
    
    if source == 'AGCD': #< Create/read AGCD mask
        y_start = int(ds_e0.time.dt.strftime("%Y")[0])
        y_end = int(ds_e0.time.dt.strftime("%Y")[-1])
        
        mask_file = out_dir+f"agcd_v1-0-1_precip_weight_r005_monthly_{y_start}_{y_end}fraction_ge1.nc"
        q_mask = get_agcd_q_mask(str(y_start),str(y_end),mask_file,ds_e0)

    #< Resample daily/monthly data to annual to achieve equal temporal length and resolution
    print("=== resample data ===")
    if source == 'CMIP':
        da_e0 = ds_e0.resample(time='YE').sum('time')[var_e0]
        da_p  = ds_p.resample(time='YE').sum('time')[var_p]
    else:
        da_e0 = ds_e0.resample(time='YE').sum('time')[var_e0]
        da_p  = ds_p.resample(time='YE').sum('time')[var_p]

    
    if source == 'AGCD' or source == 'ERA5-land':
        print('Rename latitude and longitude.')
        da_e0 = da_e0.rename({"latitude": "lat"})
        da_e0 = da_e0.rename({"longitude": "lon"})

    da_e0 = da_e0.chunk(dict(time="auto",lat=-1,lon=-1))
    da_p  = da_p.chunk(dict(time="auto",lat=-1,lon=-1))


    if source == 'AGCD':
        # #< Format time array to be identical
        da_e0["time"] = da_e0["time"].dt.strftime("%Y")
        da_p["time"] = da_p["time"].dt.strftime("%Y")
        
        #< Format and crop lat and lon to be identical
        da_e0["lat"] = da_e0["lat"].astype("float32")
        da_e0["lon"] = da_e0["lon"].astype("float32")
        da_p = da_p.sel(lon=slice(int(da_e0.lon[0]),int(da_e0.lon[-1]))).astype("float32")
        da_p = da_p.sortby(da_p["lat"], ascending=False).sel(lat=slice(int(da_e0.lat[0]),int(da_e0.lat[-1]))).astype("float32")

    return (da_p/da_e0).mean('time').astype("float32").rename("AI")


'''
######################################################################
## BUREAU OF METEOROLOGY
## AUSTRALIAN CLIMATE SERVICE
## DROUGHT HAZARD TEAM
##
## DATE:         June-2024
## SCRIPT:       utils.py
## AUTHOR:       jessica.bhardwaj@bom.gov.au
##
## PURPOSE:      Core functions and dictionaries for spi and percentiles script.
##
######################################################################
'''

# Dictionary that contains the CMIP6 model metadata. 
data_source = {
    'ERA5':{'var_p':'tp', 'var_sm':'swvl{}', 'var_pet':'evspsblpot','var_lat':'lat','var_lon':'lon'}, #swvl1+swvl2+swvl3 is volume of water in 1m soil column
    'AGCD':{'var_p':'precip', 'var_lat':'latitude','var_lon':'longitude'},
    'AWRA':{'var_sm':'s{}', 'var_pet':'e0','var_lat':'latitude','var_lon':'longitude'}, #s0 or ss
    'CMIP6':{'var_p':'pr', 'var_sm':'mrsos', 'var_et':'evspsbl', 'var_pet':'evspsblpot', 'var_lat':'lat','var_lon':'lon',
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

# Dictionary that contains GADI project file path bases for relevant datasets. 
file_paths = {'AGCD': "/g/data/zv2/agcd/v2-0-1/precip/total/r005/01month",
              'AWRA': "/g/data/iu04/australian-water-outlook/historical/v1/AWRALv7",
              'ERA5': "/g/data/zz93/era5-land/reanalysis/{}/{}", #variable #year
              'BARPA-R': "/g/data/py18/BARPA/output/CMIP6/DD/AUS-15/BOM/{}/{}/{}/BARPA-R/v1-r1/mon/{}/v20231001", #model #hist/ssp #variantid #variable
              'CCAM-v2203-SN': "/g/data/hq89/CCAM/output/CMIP6/DD/AUS-10i/CSIRO/{}/{}/{}/CCAM-v2203-SN/v1-r1/mon/{}/v20231206", #model #hist/ssp #variantid #variable
              'bias-correction': "/g/data/ia39/australian-climate-service/test-data/CORDEX-CMIP6/bias-adjustment-{}/AGCD-05i/{}/{}/{}/{}/{}/{}/day/{}" #output/input #BOM/CSIRO #model #hist/ssp #variantid #RCM #BCdetails #variable
             }


def load_target_variable(target_variable, RCM, model, bc=False, bc_method=None, bc_source= None):
    """
    Function to create dictionaries with tree structure of period/RCMs/models with relevant target variable grids as dictionary values.

    Parameters:
    - target_variable: string from ['var_p', 'var_sm', 'var_pet'] 
        (variable to accumulate - NOTE function only tested for var_pr so far.)
    - RCM: string from ['BARPA-R', 'CCAM-v2203-SN'] 
        (BARPA/CCAM RCM choice.)
    - model: string of relevant GCM 
        (CMIP6 Global Climate Model)
    - bc: string from ['input', 'output', 'raw'] 
        (string switch to access raw (hq89/py18) preprocessed bc (ia39) or bc data (ia39)
    - bc_method: string from ['QME', MRNBC'] 
        (bias correction method - NOTE MRNBC is not live yet for all ensembles)
    - bc_source: string from ['AGCD', 'BARRA'] 
        (bias correction source)

    Returns:
    - output_xr: xarray of target variable for specified model and RCM
    """
    import xarray as xr
    import glob

    climstart = 1960
    climend = 2100

    if bc == 'raw':
        file_path_base = file_paths[RCM]
        files=[]
        cmip6_hist = file_path_base.format(model,'historical',data_source['CMIP6'][model]['variant-id'], data_source['CMIP6'][target_variable])
        cmip6_ssp370 = file_path_base.format(model,'ssp370',data_source['CMIP6'][model]['variant-id'], data_source['CMIP6'][target_variable])
        for i in range(climstart,climend+1):
            files.extend(sorted(glob.glob("{}/*{}12.nc".format(cmip6_hist, i))))
            files.extend(sorted(glob.glob("{}/*{}12.nc".format(cmip6_ssp370, i))))
        for file in files:
            print(f"{file}")
        output_ds = xr.open_mfdataset(files, combine='nested', concat_dim='time', parallel=True).chunk({"time": -1})[data_source['CMIP6'][target_variable]]
        
    else:
        target_variable_key = data_source['CMIP6'][target_variable] if bc == 'input' else data_source['CMIP6'][target_variable]+'Adjust'
        file_path_base = file_paths['bias-correction']
        files=[]
        cmip6_hist = file_path_base.format(bc,\
                                            'BOM' if 'BARPA' in RCM else 'CSIRO',\
                                            model,\
                                            'historical',\
                                            data_source['CMIP6'][model]['variant-id'],\
                                            RCM,\
                                            'v1-r1' if bc == 'input' else 'v1-r1-ACS-{}-{}-{}-2022'.format(bc_method, bc_source, '1960' if bc_source == 'AGCD' else '1979'),
                                            target_variable_key)
        cmip6_ssp370 = cmip6_hist.replace('historical', 'ssp370')
        for i in range(climstart,climend+1):
            files.extend(sorted(glob.glob("{}/*{}".format(cmip6_hist, str(i)+'12.nc' if bc == 'raw' else str(i)+'1231.nc'))))
            files.extend(sorted(glob.glob("{}/*{}".format(cmip6_ssp370, str(i)+'12.nc' if bc == 'raw' else str(i)+'1231.nc'))))

        for file in files:
            print(f"{file}")
        #bc files are of daily frequency 
        output_ds = xr.open_mfdataset(files)[target_variable_key]
                  

    return output_ds
