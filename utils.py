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
    'AWRA':{'var_sm':'s{}', 'var_pet':'e0','var_lat':'latitude','var_lon':'longitude'}, 
    'CMIP6':{'var_p':'pr', 'var_tasmax':'tasmax', 'var_tasmin':'tasmin', 'var_lat':'lat','var_lon':'lon',
             'ACCESS-CM2': {'variant-id':'r4i1p1f1', 'UQ-DES':'CCAMoc-v2112'},
             'ACCESS-ESM1-5': {'variant-id':'r6i1p1f1', 'UQ-DES':'CCAM-v2105'},
             'CMCC-ESM2': {'variant-id':'r1i1p1f1', 'UQ-DES':'CCAM-v2105'}, 
             'CESM2':{'variant-id':'r11i1p1f1'},
             'CNRM-CM6-1-HR': {'variant-id':'r1i1p1f2', 'UQ-DES':'CCAM-v2112'},
             'CNRM-ESM2-1': {'variant-id':'r1i1p1f2'},
             'EC-Earth3': {'variant-id':'r1i1p1f1', 'UQ-DES':'CCAM-v2105'},
             'EC-Earth3-Veg': {'variant-id':'r1i1p1f1'},
             'FGOALS-g3': {'variant-id':'r4i1p1f1', 'UQ-DES':'CCAM-v2105'},
             'GFDL-ESM4': {'variant-id':'r1i1p1f1', 'UQ-DES':'CCAM-v2105'},
             'GISS-E2-1-G': {'variant-id':'r2i1p1f2', 'UQ-DES':'CCAM-v2105'},
             'MPI-ESM1-2-HR': {'variant-id':'r1i1p1f1',}, 
             'MPI-ESM1-2-LR': {'variant-id':'r9i1p1f1', 'UQ-DES':'CCAM-v2105'}, 
             'MRI-ESM2-0': {'variant-id':'r1i1p1f1', 'UQ-DES':'CCAM-v2105'}, 
             'NorESM2-MM': {'variant-id':'r1i1p1f1', 'UQ-DES':'CCAM-v2112'},
             'UKESM1-0-LL': {'variant-id':'r1i1p1f2'}
        }
    }

# Dictionary that contains GADI project file path bases for relevant datasets. 
file_paths = {'AGCD': "/g/data/zv2/agcd/v2-0-1/precip/total/r005/01month",
              'AWRA': "/g/data/iu04/australian-water-outlook/historical/v1/AWRALv7",
              'ERA5': "/g/data/zz93/era5-land/reanalysis/{}/{}", #variable #year
              'BARPA-R': "/g/data/py18/BARPA/output/CMIP6/DD/AUS-15/BOM/{}/{}/{}/BARPA-R/v1-r1/mon/{}/v20231001", #model #hist/ssp #variantid #variable
              'CCAM-v2203-SN': "/g/data/hq89/CCAM/output/CMIP6/DD/AUS-10i/CSIRO/{}/{}/{}/CCAM-v2203-SN/v1-r1/mon/{}/v20231206", #model #hist/ssp #variantid #variable
              'NARCliM': "/g/data/zz63/NARCliM2-0/output/CMIP6/DD/AUS-18/NSW-Government/{}/{}/{}/{}/v1-r1/mon/{}/v20240312", #model #hist/ssp #variantid #NARCliM configuration #variable
              'UQ-DES': "/g/data/ig45/QldFCP-2/output/CMIP6/DD/AUS-10i/UQ-DEC/{}/{}/{}/{}/v1-r1/mon/{}/v20240709", #model #hist/ssp #variantid #UQ-DES downscaling RCM #variable
              'bias-correction': "/g/data/ia39/australian-climate-service/release/CORDEX/output-Adjust/CMIP6/bias-adjusted-{}/AUST-05i/{}/{}/{}/{}/{}/{}/day/{}/v20241216" #output/input #RCM #model #hist/ssp #variantid #RCM #BCdetails #variable
             }


def load_target_variable(dataset_source, target_variable, ssp, RCM, model, accumulation, bc=False, bc_method=None, bc_source= None):
    """
    Function to create dictionaries with tree structure of period/RCMs/models with relevant target variable grids as dictionary values.

    Parameters:
    - target_variable: string from ['var_p', 'var_sm', 'var_pet'] 
        (variable to accumulate - NOTE function only tested for var_pr so far.)
    - RCM: string from ['BARPA-R', 'CCAM-v2203-SN'] 
        (BARPA/CCAM RCM choice.)
    - model: string of relevant GCM 
        (CMIP6 Global Climate Model)
    - accumulation: int from [1, 3, 6, 9, 12] 
        (no. of months to accumulate)
    - bc: string from ['input', 'output', 'raw'] 
        (string switch to access raw (hq89/py18/zz63/ig45) or bc input/output data (ia39)
    - bc_method: string from ['QME', MRNBC'] 
        (bias correction method - NOTE MRNBC is not live yet for all ensembles)
    - bc_source: string from ['AGCD', 'BARRA'] 
        (bias correction source)

    Returns:
    - output_xr: xarray of target variable for specified model and RCM
    """
    import xarray as xr
    import glob
    import logging
    (logging.getLogger('flox')).setLevel(logging.WARNING)
    
    climstart = 1960
    climend = 2100
    data_source['CMIP6']['ACCESS-CM2']['variant-id']='r2i1p1f1' if RCM=='UQ-DES' else data_source['CMIP6']['ACCESS-CM2']['variant-id']
    if dataset_source =='AGCD':
        climstart = 1901
        climend = 2022
        file_path_base = file_paths['AGCD']
        files = [file_path_base + "/agcd_v2-0-1_precip_total_r005_monthly_{}.nc".format(i) for i in range(climstart-1,climend+1)] 
        target_period  = xr.open_mfdataset(files ,combine='nested', concat_dim='time',parallel=True).chunk(dict(time="auto"))[data_source['AGCD'][target_variable]]
        output_xr = target_period.rolling(time=accumulation, center=False).sum().sel(time=slice(target_period.time[12], None))
    elif dataset_source =='CMIP6' and bc == 'raw':
        file_path_base = file_paths[RCM]
        files=[]
        cmip6_hist = file_path_base.format(model,'historical',data_source['CMIP6'][model]['variant-id'], data_source['CMIP6'][target_variable])
        cmip6_ssp = file_path_base.format(model, ssp, data_source['CMIP6'][model]['variant-id'], data_source['CMIP6'][target_variable])
        for i in range(climstart,climend+1):
            files.extend(sorted(glob.glob("{}/*{}12.nc".format(cmip6_hist, i))))
            files.extend(sorted(glob.glob("{}/*{}12.nc".format(cmip6_ssp, i))))
        target_period = xr.open_mfdataset(files, combine='nested', concat_dim='time', parallel=True).chunk({"time": -1})[data_source['CMIP6'][target_variable]]
        # convert precip flux (kg m-2 s-1) into mm/month 
        target_period = target_period * 86400 * 30
        # accumulate and write to dictionary
        output_xr = target_period.rolling(time=accumulation, center=False).sum().sel(time=slice(target_period.time[12], None))
    elif dataset_source =='CMIP6' and bc in ['input', 'output']:
        target_variable_key = data_source['CMIP6'][target_variable] if bc == 'input' else data_source['CMIP6'][target_variable]+'Adjust'
        file_path_base = file_paths['bias-correction']
        files=[]
        cmip6_hist = file_path_base.format(bc,\
                                            'BOM' if 'BARPA' in RCM else 'CSIRO' if 'CCAM' in RCM else 'NSW-Government' if 'NARCliM' in RCM else 'UQ-DES',\
                                            model,\
                                            'historical',\
                                            data_source['CMIP6'][model]['variant-id'],\
                                            data_source['CMIP6'][model][RCM] if RCM=='UQ-DES' else RCM,\
                                            'v1-r1' if bc == 'input' else 'v1-r1-ACS-{}-{}-{}-2022'.format(bc_method, bc_source, '1960' if bc_source == 'AGCDv1' else '1980'),
                                            target_variable_key)
        cmip6_ssp = cmip6_hist.replace('historical', ssp)
        for i in range(climstart,climend+1):
            files.extend(sorted(glob.glob("{}/*{}".format(cmip6_hist, str(i)+'12.nc' if bc == 'raw' else str(i)+'1231.nc'))))
            files.extend(sorted(glob.glob("{}/*{}".format(cmip6_ssp, str(i)+'12.nc' if bc == 'raw' else str(i)+'1231.nc'))))
    
        #bc files are of daily frequency - load in and resample to calendar month
        target_period = xr.open_mfdataset(files)[target_variable_key].resample(time='ME').sum()
        # accumulate and write to dictionary
        output_xr = target_period.rolling(time=accumulation, center=False).sum().sel(time=slice(target_period.time[12], None))
    return output_xr