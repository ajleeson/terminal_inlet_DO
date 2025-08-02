"""
Plot time series of Puget Sound hypoxic volume
"""

import matplotlib.dates as mdates
import numpy as np
import xarray as xr
import matplotlib.patches as patches
import pandas as pd
import cmocean
import matplotlib.pylab as plt

import helper_functions


def hypoxic_volume(grid_ds,hyp_vol_dict,PSbox_ds):

    years =  ['2014','2015','2016','2017','2018','2019']

    ##############################################################
    ##       Get masked land / water map of Puget Sound         ##
    ##############################################################

    # get the grid data
    z = -grid_ds.h.values
    mask_rho = np.transpose(grid_ds.mask_rho.values)
    lon = grid_ds.lon_rho.values
    lat = grid_ds.lat_rho.values
    plon, plat = helper_functions.get_plon_plat(lon,lat)
    # make a version of z with nans where masked
    # this gives us a binary map of land and water cells
    zm = z.copy()
    zm[np.transpose(mask_rho) == 0] = np.nan
    zm[np.transpose(mask_rho) != 0] = -1

    ##############################################################
    ##             Plot hypoxic volume time series              ##
    ##############################################################

    # Puget Sound bounds
    xmin = -123.29
    xmax = -122.1
    ymin = 46.95
    ymax = 48.93
    
    # colors for plotting
    colors = ['#62B6CB','#A8C256','#96031A','#957FEF','#F9627D','#476ad1','darkorange']

    # initialize figure
    fig, (ax0, ax1) = plt.subplots(1,2,figsize = (12,5.5),gridspec_kw={'width_ratios': [1, 4]})

    # format figure
    ax0.set_xlim([xmin,xmax])
    ax0.set_ylim([ymin,ymax])
    ax0.set_ylabel('Latitude', fontsize=12)
    ax0.set_xlabel('Longitude', fontsize=12)
    ax0.tick_params(axis='both', labelsize=12)
    ax0.pcolormesh(plon, plat, zm, vmin=-8, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))
    helper_functions.dar(ax0)
    # Create a Rectangle patch to omit Straits
    # get lat and lon
    lons = PSbox_ds.coords['lon_rho'].values
    lats = PSbox_ds.coords['lat_rho'].values
    lon = lons[0,:]
    lat = lats[:,0]
    # Straits
    lonmax = -122.76
    lonmin = xmin
    latmax = ymax
    latmin = 48.14
    # convert lat/lon to eta/xi
    diff = np.absolute(lon-lonmin)
    ximin = diff.argmin()
    diff = np.absolute(lon-lonmax)
    ximax = diff.argmin()
    diff = np.absolute(lat-latmin)
    etamin = diff.argmin()
    diff = np.absolute(lat-latmax)
    etamax = diff.argmin()
    rect = patches.Rectangle((lon[ximin], lat[etamin]), lon[ximax]-lon[ximin], lat[etamax]-lat[etamin],
                            edgecolor='none', facecolor='white', alpha=0.9)
    # Add the patch to the Axes
    ax0.add_patch(rect)
    ax0.text(-123.2,48.25,'Straits\nomitted', rotation=90, fontsize=12)
    ax0.set_title('(a)', fontsize = 14, loc='left', fontweight='bold')

    # Puget Sound volume with straits omitted
    PS_vol = 195.2716230839466 # [km^3]

    # create time vector
    startdate = '2020.01.01'
    enddate = '2020.12.31'
    dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
    dates_local = [helper_functions.get_dt_local(x) for x in dates]

    # plot timeseries
    for i,year in enumerate(years):
        if year == '2017':
            ax1.plot(dates_local,hyp_vol_dict[year],color='white',
                    linewidth=4,zorder=4)
            ax1.plot(dates_local,hyp_vol_dict[year],color='black',
                    linewidth=3.5,label=year,zorder=4)
        else:
            ax1.plot(dates_local,hyp_vol_dict[year],color='white',
                    linewidth=2.5)
            ax1.plot(dates_local,hyp_vol_dict[year],color=colors[i],
                    linewidth=2,label=year)

    # get median hypoxic volume
    med_vol = np.nanmedian(list(hyp_vol_dict.values()), axis=0)
    ax1.plot(dates_local,med_vol,color='k',
            linestyle='--',linewidth=2,label='median')

    # format figure
    # ax1.grid(visible=True, axis='x', color='w')
    ax1.grid(visible=True, axis='both', color='silver', linestyle='--')
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax1.tick_params(axis='both', labelsize=12)
    ax1.set_ylabel(r'Hypoxic volume [km$^3$]', fontsize=12)
    plt.legend(loc='upper left', fontsize=12)
    plt.title('(b)', fontsize = 14, loc='left', fontweight='bold')
    ax1.set_xlim([dates_local[0],dates_local[-1]])
    ax1.set_ylim([0,13])

    # convert hypoxic volume to percent hypoxic volume
    percent = lambda hyp_vol: hyp_vol/PS_vol*100
    # get left axis limits
    ymin, ymax = ax1.get_ylim()
    # match ticks
    ax2 = ax1.twinx()
    ax2.set_ylim((percent(ymin),percent(ymax)))
    ax2.plot([],[])
    for border in ['top','right','bottom','left']:
        ax2.spines[border].set_visible(False)
    ax2.set_ylabel(r'Percent of regional volume [%]', fontsize=12)

    plt.show()

    return