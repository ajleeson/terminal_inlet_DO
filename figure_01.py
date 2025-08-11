"""
Creat map of LiveOcean's Salish Sea and Puget Sound bathymetry
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
import cmocean

import helper_functions

def model_bathy(grid_ds):

    background = 'white'

    # Get LiveOcean grid info --------------------------------------------------

    # get the grid data
    z = -grid_ds.h.values
    mask_rho = np.transpose(grid_ds.mask_rho.values)
    lon = grid_ds.lon_rho.values
    lat = grid_ds.lat_rho.values
    X = lon[0,:] # grid cell X values
    Y = lat[:,0] # grid cell Y values
    plon, plat = helper_functions.get_plon_plat(lon,lat)
    # make a version of z with nans where masked
    zm = z.copy()
    zm[np.transpose(mask_rho) == 0] = np.nan

    # Create bathymetry plot --------------------------------------------------------------

    # Initialize figure
    fig = plt.figure(figsize = (11,7))
    plt.subplots_adjust(wspace=0, hspace=0)
    newcmap = cmocean.tools.crop_by_percent(cmocean.cm.thermal_r, 20, which='min')
    newcmap.set_bad(background,1.) # background color

    # Salish Sea ----------------------------------------------------------
    ax0 = fig.add_subplot(1,2,1)
    cs = ax0.pcolormesh(plon, plat, zm*-1, vmin=0, vmax=250, cmap=newcmap)
    helper_functions.dar(ax0)
    # Set axis limits
    ax0.set_xlim(-124.98549,-122)  # Salish Sea
    ax0.set_ylim(46.8165519,50.39679) # Salish Sea

    plt.xticks(rotation=30,horizontalalignment='right',fontsize=12)
    plt.yticks(fontsize=12)
    ax0.set_ylabel('Latitude', fontsize=12)
    ax0.set_xlabel('Longitude', fontsize=12)
    ax0.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax0.yaxis.set_major_locator(MaxNLocator(integer=True))
    # add title
    ax0.set_title('(a) Salish Sea',fontsize=14, loc='left', fontweight='bold')#,color='#EEEEEE')
    # add 10 km bar
    lat0 = 47
    lon0 = -124.63175
    lat1 = lat0
    lon1 = -124.36975
    distances_m = helper_functions.ll2xy(lon1,lat1,lon0,lat0)
    x_dist_km = round(distances_m[0]/1000)
    ax0.plot([lon0,lon1],[lat0,lat1],color='white',linewidth=5)
    ax0.text((lon0+lon1)/2,lat0+0.05,'{} km'.format(x_dist_km),color='w',fontsize=12,
            horizontalalignment='center')
    # draw box around Puget Sound
    bordercolor = 'k'
    ax0.add_patch(Rectangle((-123.3, 46.93), 1.2, 1.52,
                edgecolor = bordercolor, facecolor='none', lw=1))

    # add major cities
    # Seattle
    ax0.scatter([-122.3328],[47.6061],s=[250],color='pink',
                marker='*',edgecolors='darkred')
    ax0.text(-122.3328 + 0.1,47.6061,'Seattle',color='darkred', rotation=90,
            horizontalalignment='left',verticalalignment='center', size=12)
    # Vancouver
    ax0.scatter([-123.1207],[49.2827],s=[250],color='pink',
                marker='*',edgecolors='darkred')
    ax0.text(-123.1207 + 0.1,49.2827,'Vancouver',color='darkred', rotation=0,
            horizontalalignment='left',verticalalignment='center', size=12)

    # add major water bodies
    ax0.text(-124.937095,47.782238,'Pacific Ocean',color='k', rotation=-75,
            horizontalalignment='left',verticalalignment='center', size=12,fontweight='bold')

    # Puget Sound ----------------------------------------------------------
    ax1 = fig.add_subplot(1,2,2)
    cs = ax1.pcolormesh(plon, plat, zm*-1, vmin=0, vmax=250, cmap=newcmap)
    cbar = plt.colorbar(cs,ax=ax1, location='right', pad=0.05)
    cbar.ax.tick_params(labelsize=11)
    cbar.ax.set_ylabel('Depth [m]', fontsize=11)
    cbar.outline.set_visible(False)
    # format figure
    helper_functions.dar(ax1)
    # Set axis limits
    ax1.set_xlim([-123.3,-122.1])
    ax1.set_ylim([46.93,48.45])
    ax1.set_yticklabels([])
    ax1.set_xticklabels([])
    # add title
    ax1.set_title('(b) Puget Sound',fontsize=14,loc='left', fontweight='bold')# color='#EEEEEE')
    # add 10 km bar
    lat0 = 47
    lon0 = -123.2
    lat1 = lat0
    lon1 = -123.06825
    distances_m = helper_functions.ll2xy(lon1,lat1,lon0,lat0)
    x_dist_km = round(distances_m[0]/1000)
    ax1.plot([lon0,lon1],[lat0,lat1],color='k',linewidth=5)
    ax1.text((lon0+lon1)/2,lat0+0.015,'{} km'.format(x_dist_km),color='k',fontsize=12,
            horizontalalignment='center')

    plt.subplots_adjust(wspace = -0.3)
    plt.show()

    return
