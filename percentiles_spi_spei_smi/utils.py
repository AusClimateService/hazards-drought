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
    'CMIP6':{'var_p':'pr', 'var_sm':'evspsblpot', 'var_pet':'evspsblpot', 'var_lat':'lat','var_lon':'lon',
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


def load_target_variable(target_variable, RCM, model, accumulation, bc=False, bc_method=None, bc_source= None):
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
        target_period = xr.open_mfdataset(files, combine='nested', concat_dim='time', parallel=True).chunk({"time": -1})[data_source['CMIP6'][target_variable]]
        # convert precip flux (kg m-2 s-1) into mm/month 
        target_period = target_period * 86400 * 30
        # accumulate and write to dictionary
        output_xr = target_period.rolling(time=accumulation, center=False).sum().sel(time=slice(target_period.time[12], None))
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
            files.extend(sorted(glob.glob("{}/*{}".format(cmip6_hist, str(i)+'1231.nc' if bc == 'output' or RCM=='CCAM-v2203-SN' else str(i)+'12.nc'))))
            files.extend(sorted(glob.glob("{}/*{}".format(cmip6_ssp370, str(i)+'1231.nc' if bc == 'output' or RCM=='CCAM-v2203-SN' else str(i)+'12.nc'))))

        #bc files are of daily frequency - load in and resample to calendar month
        target_period = xr.open_mfdataset(files)[target_variable_key].resample(time='ME').sum()
        # accumulate and write to dictionary
        output_xr = target_period.rolling(time=accumulation, center=False).sum().sel(time=slice(target_period.time[12], None))
                    

    return output_xr