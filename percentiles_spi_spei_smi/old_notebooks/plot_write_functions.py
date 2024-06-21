def plot_write_percentile_threshold(source, period, input_grid, percentile_threshold, accumulation, **kwargs):

    import dictionaries
    import xarray as xr
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import numpy as np
    from core_functions import mask_ocean
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    BC = kwargs.get('bias_correction')
    write_netcdf = kwargs.get('write_netcdf')
    
    # load data
    if source == 'AGCD':
        percentile_thresh = xr.DataArray(np.percentile(input_grid, percentile_threshold, axis=0), coords=input_grid[0].coords, dims=input_grid[0].dims)
    else:
        input_dict = {k:v for k,v in input_grid[period][source].items() if k != 'AGCD'}
        for model in list(input_dict.keys()):
            input_dict[model] = xr.DataArray(np.percentile(input_dict[model], percentile_threshold, axis=0), coords=input_dict[model][0].coords, dims=input_dict[model][0].dims)
            if write_netcdf:
                file_str = '%s_%s_%s_%sp_%s_%s.nc' % (source,model,dictionaries.data_source['CMIP6'][model]['variant-id'],percentile_threshold,period,'RAW' if BC==False else BC)
                input_dict[model].rename('R%sp' %(percentile_threshold)).to_netcdf(file_str,mode='w')
        models_concat = xr.concat(list(input_dict.values()), dim='model_mean')
        percentile_thresh = models_concat.mean(dim='model_mean')
        percentile_thresh = percentile_thresh.interp(lat=input_grid['recent']['AGCD'].lat, lon=input_grid['recent']['AGCD'].lon)
    
    # cmap - to match OCHS colour scheme (http://www.bom.gov.au/climate/maps/rainfall/?variable=rainfall&map=totals&period=3month)
    colors = ['#FFFFFF', '#FFBF59', '#FFAC01', '#FFFF00', '#B2FF00', '#4CFF00', '#00E599', '#01A5FF', '#3F3FFF', '#B200FF', '#FF01FF', '#FF4C9B']
    ranges = [1, 5 if accumulation ==1 else 2, 10, 25, 50, 100, 200, 300, 400, 600, 900 if accumulation >=9 else 800, 1200] 
    range_idx = accumulation if accumulation <= 9 else 9  # adapt colorbar ranges depending on accumulation
    ranges = [0] + ranges[range_idx // 3:] + [1600, 2000, 2400][:range_idx // 3]
    norm = mcolors.BoundaryNorm(ranges, len(colors))
    cmap = mcolors.ListedColormap(colors, 'custom_cmap', len(colors))
    
    # plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    p15_plot = ax.contourf(percentile_thresh.lon, percentile_thresh.lat, percentile_thresh, norm=norm, cmap=cmap, levels=ranges, extend='max', transform=ccrs.PlateCarree())
    cbar = plt.colorbar(p15_plot, ax=ax, ticks=ranges)
    cbar.set_label('%s-month %s %sth percentile threshold (mm)' % (str(accumulation), source, str(percentile_threshold)))
    cbar.ax.set_yticklabels(ranges, fontsize=8)
    ax.contour(percentile_thresh.lon, percentile_thresh.lat, percentile_thresh, colors='black', linewidths=0.2, levels=ranges, transform=ccrs.PlateCarree())
    ocean_mask = cfeature.NaturalEarthFeature('physical', 'ocean', '10m', edgecolor='face')
    ax.add_feature(ocean_mask, facecolor='white', zorder=10)
    ax.coastlines(linewidth=0.2, zorder=100)
    ax.spines['geo'].set_edgecolor('white')
    
    #labels and titles
    ax.set_xlabel('Longitude', labelpad=17)
    ax.set_ylabel('Latitude', labelpad=28)
    title = '%s\n15th percentile threshold for %s period' if BC == False else 'Bias-Corrected %s\n15th percentile threshold for %s period'
    plt.title(title % (source, period), pad=5)
    if source == 'AGCD':
        plt.text(113, -43, 'Source: AGCDv2\nBase period: %s-%s' % (dictionaries.climatology[period]['start'], dictionaries.climatology[period]['end']), fontsize=10, color='black', zorder=10)
    elif 'BARPA' in source:
        plt.text(113, -43, 'Source models: CMCC-ESM2\nACCESS-ESM1-5, ACCESS-CM2, EC-Earth3,\nCESM2, MPI-ESM1-2-HR, NorESM2-MM\n \nBase period: %s' % ('Recent' if period == 'recent' else 'Global Warming Level '+period[3:6]+'°C'), fontsize=10, color='black', zorder=10)
    elif 'CCAM' in source:
        plt.text(113, -43, 'Source models: CMCC-ESM2\nACCESS-ESM1-5, ACCESS-CM2, EC-Earth3,\nCESM2, CNRM-ESM2-1\n \nBase period: %s' % ('Recent' if period == 'recent' else 'Global Warming Level '+period[3:6]+'°C'), fontsize=10, color='black', zorder=10)
    
    
    plt.show()
    plt.clf()
    return

def plot_comparison_box_whisker(comparison_contemporary, input_array, RCM, NCRA_dict):
    from matplotlib.patches import Patch
    import seaborn as sns
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    import warnings

    warnings.filterwarnings("ignore", category=RuntimeWarning, message="All-NaN slice encountered")
    
    sns.set()
    box_plot_data = []
    for NCRA_region in ['NSW_&_ACT', 'NT', 'QLD_North', 'QLD_South', 'SA', 'TAS', 'VIC', 'WA_North', 'WA_South']:
        for comparison_GWL in ['GWL1.5-ssp370', 'GWL2.0-ssp370', 'GWL3.0-ssp370']:
            gwl_list = []
            pct_change_list = []
            for model in input_array['GWL1.2-ssp370'][RCM].keys():
                GWL_median =np.nanmedian(NCRA_dict[NCRA_region][comparison_GWL][model].astype(np.float64), axis=0)
                recent_median = np.nanmedian(NCRA_dict[NCRA_region][comparison_contemporary][model].astype(np.float64), axis=0)
                pct_change = np.array(100*(GWL_median-recent_median)/recent_median).flatten()
                pct_change_list.extend(pct_change)
                gwl_list.extend([comparison_GWL] * len(pct_change))
            box_plot_data.extend(zip(gwl_list, [NCRA_region.replace('_', ' ')]*len(pct_change_list), pct_change_list))
    
    box_plot_df = pd.DataFrame(box_plot_data, columns=['GWL', 'Region', 'PctChange'])
    
    plt.figure(figsize=(10, 5))
    sns.boxplot(data=box_plot_df, x='Region', y='PctChange', hue='GWL', whis=(10, 90), showfliers=False, palette='Set2', width=0.8, dodge=True)
    comparison_contemporary = comparison_contemporary[:-7] if comparison_contemporary!= 'recent' else comparison_contemporary
    plt.title('SSP3.70 CMIP6/%s percent change in median rainfall between GWL and %s period' % (RCM, comparison_contemporary), size=12)
    plt.xlabel('Region', size=12)
    plt.subplots_adjust(bottom=0.5)
    plt.ylabel('Percent Change (%)', size=12)
    legend_labels = ['1.5°C', '2.0°C', '3.0°C']
    legend_handles = [Patch(color=color) for color in sns.color_palette('Set2', len(legend_labels))]
    plt.legend(title='Global Warming Level', labels=legend_labels, handles=legend_handles, ncol=3)
    
    plt.tight_layout()
    plt.show()
    return


def plot_comparison_heatmap(comparison_contemporary, input_array, RCM, NCRA_dict):
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    import seaborn as sns
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    import warnings

    warnings.filterwarnings("ignore", category=RuntimeWarning, message="All-NaN slice encountered")
    
    sns.set()
    plot_df_list = []
    for comparison_GWL in ['GWL1.5-ssp370', 'GWL2.0-ssp370', 'GWL3.0-ssp370']:
            plot_df = pd.DataFrame(columns=list(input_array[comparison_contemporary][RCM].keys()), index=list(NCRA_dict.keys())).astype(float)
            for model in input_array['GWL1.2-ssp370'][RCM].keys():
                for NCRA_region in ['NSW_&_ACT', 'NT', 'QLD_North', 'QLD_South', 'SA', 'TAS', 'VIC', 'WA_North', 'WA_South']:
                    GWL_median =np.nanmedian(NCRA_dict[NCRA_region][comparison_GWL][model].astype(np.float64))
                    recent_median = np.nanmedian(NCRA_dict[NCRA_region][comparison_contemporary][model].astype(np.float64))
                    plot_df[model][NCRA_region]=100*(GWL_median-recent_median)/recent_median
    
            plot_df_list.append(plot_df.drop('AGCD', axis=1) if 'AGCD' in list(plot_df.columns) else plot_df)
            
    colors = sns.color_palette("BrBG", 12)
    cmap = ListedColormap(colors)
    sns.set(font_scale=0.8)
    plt.figure(figsize=(13, 4))
    sns.heatmap(pd.concat(plot_df_list, axis=1), annot=True, cmap=cmap, center=0, vmin=-30, vmax=30, linewidth = 2, cbar_kws={'extend':'both', 'label':'Percent Change (%)',}, yticklabels=[key.replace('_', ' ') for key in NCRA_dict.keys()])
    comparison_contemporary = comparison_contemporary[:-7] if comparison_contemporary!= 'recent' else comparison_contemporary
    plt.title("SSP3.70 %s/CMIP6 percent change in median rainfall between GWL period and %s period" % (RCM, comparison_contemporary), pad=18, size=12)
    
    plt.annotate('GWL1.5', xy=((1/3)/2, -0.5), size=12, xycoords='axes fraction', ha='center', va='center')
    plt.annotate('GWL2.0', xy=((1/3)+(1/3)/2, -0.5), size=12, xycoords='axes fraction', ha='center', va='center')
    plt.annotate('GWL3.0', xy=((2/3)+(1/3)/2, -0.5), size=12, xycoords='axes fraction', ha='center', va='center')
    plt.annotate('|', xy=(0, -0.5), xycoords='axes fraction', ha='center', va='center', size=14)
    plt.annotate('|', xy=(1/3, -0.5), xycoords='axes fraction', ha='center', va='center', size=14)
    plt.annotate('|', xy=(2/3, -0.5), xycoords='axes fraction', ha='center', va='center', size=14)
    plt.annotate('|', xy=(1, -0.5), xycoords='axes fraction', ha='center', va='center', size=14)
    plt.axvline(x=0.04, ymin=-10, ymax=1, color='black', lw=1)
    no_of_models = len(input_array['GWL1.2-ssp370'][RCM].keys())
    plt.axvline(x=no_of_models , color='black', lw=1)
    plt.axvline(x=2*no_of_models, color='black', lw=1)
    plt.axvline(x=3*no_of_models, color='black', lw=1)

    plt.ylabel('NCRA Regions', size=12)
    plt.show()
    plt.clf()
    return
