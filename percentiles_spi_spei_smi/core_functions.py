# get GWL start and end years from https://github.com/mathause/cmip_warming_levels/blob/main/warming_levels/cmip6_all_ens/csv/cmip6_warming_levels_all_ens_1850_1900_grid.csv
def get_GWL_syear_eyear(CMIP,GCM,pathway,variant_id,GWL,start_end):
    import pandas as pd
    df = pd.read_csv('https://raw.githubusercontent.com/mathause/cmip_warming_levels/main/warming_levels/%s_all_ens/csv/%s_warming_levels_all_ens_1850_1900_grid.csv' % (CMIP, CMIP), skiprows=4)
    output_year = df[(df['model'] == GCM) & (df[' ensemble'] == ' ' + variant_id) & (df[' exp'] == ' ' + pathway) & (df[' warming_level'] == GWL)][' %s_year' % (start_end)].values[0]
    return output_year

    
def mask_ocean(array, lat, lon):
    from mpl_toolkits.basemap import Basemap, maskoceans
    import numpy as np
    lon_m = np.repeat((lon,), len(lat), axis=0)
    lat_m = np.repeat((lat,), len(lon),axis=0).transpose()
    array_m = maskoceans(lon_m, lat_m, array, inlands = False, resolution='f', grid=1.25) 
    return array_m, lat, lon
    
def load_target_variable(target_variable, accumulation, period_list, bc=False, bc_method=None, bc_source= None):
    """
    Function to create dictionaries with tree structure of period/RCMs/models with relevant target variable grids as dictionary values.

    Parameters:
    - target_variable: string from ['var_p', 'var_sm', 'var_pet'] (variable to accumulate - NOTE function only tested for var_pr so far.)
    - accumulation: int from [1, 3, 6, 9, 12] (no. of months to accumulate)
    - period_list: list like ['full', 'recent', 'GWL1.2-ssp370', 'GWL1.5-ssp370', 'GWL2.0-ssp370', 'GWL3.0-ssp370'] (periods to load in, consult dictionaries.climatology for details)
    - bc: True/False (switch to access bias-corrected data, true-> function will pull from ia39 otherwise py18 and hq89)
    - bc_method: string from ['QME', MRNBC'] (bias correction method - NOTE MRNBC is not live yet)
    - bc_source: string from ['AGCD', 'BARRA'] (bias correction source)

    Returns:
    - target_variable_dict: dictionary of target variable with structure: target_variable_dict[period][RCM][model]
    """
    import dictionaries
    import xarray as xr
    import glob


    print('---> BC SWITCH ON: USING DATA FROM ia39' if bc==True else '---> BC SWITCH OFF: USING DATA FROM py18 AND hq89')
    target_variable_dict = {}
    
    for period in period_list:
        print()
        print('---> LOADING %s PERIOD' % (period.upper()))
        target_variable_dict[period]={}
    
        #LOAD AGCD
        if period in ['recent', 'full']:
            print('- AGCD')
            climstart = dictionaries.climatology[period]['start'] if period == 'recent' else 1961
            climend = dictionaries.climatology[period]['end'] if period == 'recent' else 2022
            file_path_base = dictionaries.file_paths['AGCD']
            files = [file_path_base + "/agcd_v2-0-1_precip_total_r005_monthly_{}.nc".format(i) for i in range(climstart-1,climend+1)] # include one year before to allow accumulations up to 12 months - target period is still only on the defined years
            target_period  = xr.open_mfdataset(files ,combine='nested', concat_dim='time',parallel=True).chunk(dict(time="auto"))[dictionaries.data_source['AGCD']['var_p']]
            
            # accumulate and write to dictionary
            target_variable_dict[period]={}
            target_variable_dict[period]['AGCD'] = target_period.rolling(time=accumulation, center=False).sum().sel(time=slice(target_period.time[12], None))
        
        
        #LOAD BARPA-R // CCAM
        for RCM in ['BARPA-R', 'CCAM-v2203-SN']:
            #MPI-ESM1-2-HR // NorESM2-MM for BARPA only and CNRM-ESM2-1 for CCAM only
            model_list = ['CMCC-ESM2', 'ACCESS-ESM1-5', 'ACCESS-CM2', 'EC-Earth3', 'CESM2'] 
            model_list = model_list + ['CNRM-ESM2-1'] if RCM == 'CCAM-v2203-SN' else model_list + ['MPI-ESM1-2-HR', 'NorESM2-MM']
            target_variable_dict[period][RCM]={}
            for model in model_list:
                print('- %s-%s (%s)' % (RCM, model, dictionaries.data_source['CMIP6'][model]['variant-id']))
                chunk_sizes = {"time": -1} 
                # include one year before to allow accumulations up to 12 months - target period is still only on the defined years
                climstart = dictionaries.climatology[period][model]['start'] if 'GWL' in period else dictionaries.climatology[period]['start']-1
                climend = dictionaries.climatology[period][model]['end'] if 'GWL' in period else dictionaries.climatology[period]['end']
                if bc:
                    file_path_base = dictionaries.file_paths['bias-correction']
                    files=[]
                    cmip6_hist = file_path_base.format('BOM' if 'BARPA' in RCM else 'CSIRO',model,'historical',dictionaries.data_source['CMIP6'][model]['variant-id'], RCM, bc_method, bc_source, '1960' if bc_source == 'AGCD' else '1979', dictionaries.data_source['CMIP6'][target_variable])
                    cmip6_ssp370 = file_path_base.format('BOM' if 'BARPA' in RCM else 'CSIRO',model,'ssp370',dictionaries.data_source['CMIP6'][model]['variant-id'], RCM, bc_method, bc_source, '1960' if bc_source == 'AGCD' else '1979', dictionaries.data_source['CMIP6'][target_variable])
                    for i in range(climstart,climend+1):
                        files.extend(sorted(glob.glob("{}/*{}1231.nc".format(cmip6_hist, i))))
                        files.extend(sorted(glob.glob("{}/*{}1231.nc".format(cmip6_ssp370, i))))
                    #bc files are of daily frequency - load in and resample to calendar month
                    target_period = xr.open_mfdataset(files)[dictionaries.data_source['CMIP6'][target_variable]+'Adjust'].resample(time='ME').sum()
                    # accumulate and write to dictionary
                    target_variable_dict[period][RCM][model] = target_period.rolling(time=accumulation, center=False).sum().sel(time=slice(target_period.time[12], None))
                else:
                    file_path_base = dictionaries.file_paths[RCM]
                    files=[]
                    cmip6_hist = file_path_base.format(model,'historical',dictionaries.data_source['CMIP6'][model]['variant-id'], dictionaries.data_source['CMIP6'][target_variable])
                    cmip6_ssp370 = file_path_base.format(model,'ssp370',dictionaries.data_source['CMIP6'][model]['variant-id'], dictionaries.data_source['CMIP6'][target_variable])
                    for i in range(climstart,climend+1):
                        files.extend(sorted(glob.glob("{}/*{}12.nc".format(cmip6_hist, i))))
                        files.extend(sorted(glob.glob("{}/*{}12.nc".format(cmip6_ssp370, i))))
            
                    target_period = xr.open_mfdataset(files, combine='nested', concat_dim='time', parallel=True).chunk(chunk_sizes)[dictionaries.data_source['CMIP6'][target_variable]]
                    # convert precip flux (kg m-2 s-1) into mm/month 
                    target_period = target_period * 86400 * 30
                    # accumulate and write to dictionary
                    target_variable_dict[period][RCM][model] = target_period.rolling(time=accumulation, center=False).sum().sel(time=slice(target_period.time[12], None))

    return target_variable_dict


def calculate_spi(target_period, base_period, accumulation, fitting_method):
    """
    Function to calculate SPI for a nominated accumulation using monthly rainfall data.

    Parameters:
    - target_period_monthly_rain: xarray DataArray (monthly rainfall array for period of interest)
    - base_period_monthly_rain: xarray DataArray (monthly rainfall array for base period/climatology)
    - accumulation: integer from [1, 3, 6, 9, 12] (no. of months to accumulate)
    - fitting_method: string from ['MLE', 'MM'] (gamma distribution fit method - Maximum Likelihood Estimation(MLE) or Method of Moments (MM))  
    Returns:
    - spi_data: xarray DataArray (array containing SPI values)
    """
    import numpy as np
    import xarray as xr
    from scipy.stats import gamma
    import scipy
    
    # Step 1: Accumulate
    accumulated_base_period = base_period.rolling(time=accumulation, center=False).sum().sel(time=slice(base_period.time[12], None))
    accumulated_target_period = target_period.rolling(time=accumulation, center=False).sum().sel(time=slice(target_period.time[12], None))
    
    # Step 2: Compute the gamma distribution parameters (alpha, beta) using the climatology period
    alpha, loc, beta = xr.apply_ufunc(
    gamma.fit, 
    (accumulated_base_period.where(accumulated_base_period > 0)).dropna(dim='time'),
    input_core_dims=[['time']], 
    output_core_dims=[[], [], []], 
    vectorize=True,  
    dask='parallelized',
    dask_gufunc_kwargs={'allow_rechunk': True},
    output_dtypes=[float, float, float],  
    kwargs={'floc': 0, 'method':fitting_method}  
    )
    
    # Step 3: Calculate the gamma cumulative distribution function (CDF)
    gamma_cdf = gamma.cdf(accumulated_target_period, alpha, scale=beta)
    
    # Step 4: Calculate SPI values based on the inverse CDF of the normal distribution
    # (accounting for observed zero precipitation values)
    zeros = (accumulated_base_period == 0).sum(dim='time')
    probabilities_of_zero = zeros / len(base_period.time)
    probabilities = probabilities_of_zero.values + ((1 - probabilities_of_zero.values) * gamma_cdf)
    probabilities = np.where(accumulated_base_period == 0, probabilities_of_zero, probabilities)
    spi_values = scipy.stats.norm.ppf(probabilities)
    
    # Step 5: Convert SPI values to xarray DataArray
    spi_data = xr.DataArray(spi_values, coords=accumulated_target_period.coords, dims=accumulated_target_period.dims)

    return spi_data


    
