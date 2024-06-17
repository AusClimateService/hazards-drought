'''
Dictionary that contains the CMIP6 model metadata. 
'''
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

'''
Dictionary that contains GADI project file path bases for relevant datasets. 
'''
file_paths = {'AGCD': "/g/data/zv2/agcd/v2-0-1/precip/total/r005/01month",
              'AWRA': "/g/data/iu04/australian-water-outlook/historical/v1/AWRALv7",
              'ERA5': "/g/data/zz93/era5-land/reanalysis/{}/{}", #variable #year
              'BARPA-R': "/g/data/py18/BARPA/output/CMIP6/DD/AUS-15/BOM/{}/{}/{}/BARPA-R/v1-r1/mon/{}/v20231001", #model #hist/ssp #variantid #variable
              'CCAM-v2203-SN': "/g/data/hq89/CCAM/output/CMIP6/DD/AUS-10i/CSIRO/{}/{}/{}/CCAM-v2203-SN/v1-r1/mon/{}/v20231206", #model #hist/ssp #variantid #variable
              'bias-correction': "/g/data/ia39/australian-climate-service/test-data/CORDEX-CMIP6/bias-adjustment-output/AGCD-05i/{}/{}/{}/{}/{}/v1-r1-ACS-{}-{}-{}-2022/day/{}Adjust" #BOM/CSIRO #model #hist/ssp #variantid #RCM #BC_method #BC_source #1960/AGCD_1979/BARRA #variable
             }
'''
Dictionary that contains the start and end year for the different time periods. Populated using get_GWL_syear_eyear function.
'''
climatology = {'full':{'start':1960,'end':2100},
               'recent':{'start':1991,'end':2020},
               # 'current':{'start':2011,'end':2030},
               'GWL1.2-ssp370':{},
               'GWL1.5-ssp370':{},
               'GWL2.0-ssp370':{},
               'GWL3.0-ssp370':{},
               'GWL3.0-ssp370':{}
              }


from core_functions import get_GWL_syear_eyear
    
for model in list(data_source['CMIP6'].keys())[5:]:
    climatology['GWL1.2-ssp370'][model]={'start':get_GWL_syear_eyear('cmip6',model,'ssp370',data_source['CMIP6'][model]['variant-id'],1.20, 'start'),\
                                         'end':get_GWL_syear_eyear('cmip6',model,'ssp370',data_source['CMIP6'][model]['variant-id'],1.20, 'end')}
    climatology['GWL1.5-ssp370'][model]={'start':get_GWL_syear_eyear('cmip6',model,'ssp370',data_source['CMIP6'][model]['variant-id'],1.50, 'start'),\
                                         'end':get_GWL_syear_eyear('cmip6',model,'ssp370',data_source['CMIP6'][model]['variant-id'],1.50, 'end')}
    climatology['GWL2.0-ssp370'][model]={'start':get_GWL_syear_eyear('cmip6',model,'ssp370',data_source['CMIP6'][model]['variant-id'],2.00, 'start'),\
                                         'end':get_GWL_syear_eyear('cmip6',model,'ssp370',data_source['CMIP6'][model]['variant-id'],2.00, 'end')}
    climatology['GWL3.0-ssp370'][model]={'start':get_GWL_syear_eyear('cmip6',model,'ssp370',data_source['CMIP6'][model]['variant-id'],3.00, 'start'),\
                                         'end':get_GWL_syear_eyear('cmip6',model,'ssp370',data_source['CMIP6'][model]['variant-id'],3.00, 'end')}

