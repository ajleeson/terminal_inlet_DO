"""
Plots maps of six-year mean:
(a) number of days that a grid cell experiences bottom hypoxia per year
(b) bottom DO concentration during the hypoxic period (Aug 1 -  Sep 30)
"""

# import things
import numpy as np
import matplotlib.pylab as plt

import helper_functions


def pugetsound_hyp_map(grid_ds,PSbox_ds,hyp_days_dict,hyp_seas_DO_dict):

    # Puget Sound region
    xmin = -123.29
    xmax = -122.1
    ymin = 46.95
    ymax = 48.93


    # Get LiveOcean grid info 
    z = -grid_ds.h.values
    mask_rho = np.transpose(grid_ds.mask_rho.values)
    lon = grid_ds.lon_rho.values
    lat = grid_ds.lat_rho.values
    plon, plat = helper_functions.get_plon_plat(lon,lat)
    # make a version of z with nans where masked
    zm = z.copy()
    zm[np.transpose(mask_rho) == 0] = np.nan
    zm[np.transpose(mask_rho) != 0] = -1.1


    ##############################################################
    ##               PLOT DAYS WITH BOTTOM HYPOXIA              ##
    ##############################################################

    # Create map of Puget Sound (fully grey)
    fig = plt.figure(figsize=(11,9))
    ax = fig.add_subplot(1,2,1)
    plt.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-1.5, vmax=0, cmap=plt.get_cmap('Greys'))

    # get average number of days that each grid cell experiences bottom hypoxia every year
    DO_days = hyp_days_dict['avg']

    # get lat and lon for plotting
    lons = PSbox_ds.coords['lon_rho'].values
    lats = PSbox_ds.coords['lat_rho'].values
    px, py = helper_functions.get_plon_plat(lons,lats)

    # plot average number of days that each grid cell experiences bottom hypoxia every year
    cs = ax.pcolormesh(px,py,DO_days, vmin=0, vmax=np.nanmax(DO_days), cmap='rainbow')
    cbar = fig.colorbar(cs)
    cbar.ax.tick_params(labelsize=12)
    cbar.outline.set_visible(False)

    # add 10 km bar
    lat0_10k = 47#46.94
    lon0_10k = -123.05 + 0.7
    lat1_10k = lat0_10k
    lon1_10k = -122.91825 + 0.7
    distances_m = helper_functions.ll2xy(lon1_10k,lat1_10k,lon0_10k,lat0_10k)
    x_dist_km = round(distances_m[0]/1000)
    ax.plot([lon0_10k,lon1_10k],[lat0_10k,lat1_10k],color='k',linewidth=2)
    ax.text(lon0_10k-0.04,lat0_10k+0.01,'{} km'.format(x_dist_km),color='k',fontsize=12)

    # format figure
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    helper_functions.dar(ax)
    plt.xticks(rotation=30)
    plt.yticks(rotation=30)
                                    
    # Add colormap title
    ax.set_title('(a)', fontsize=12, loc='left', fontweight='bold')

    ##############################################################
    ##             PLOT MEAN BOTTOM DO CONCENTRATION            ##
    ##############################################################

    # add subplot
    ax = fig.add_subplot(1,2,2)

    # plot six-year mean bottom DO concentration during hypoxic season
    cmap = plt.cm.get_cmap('rainbow_r', 10)
    vmin = 0
    vmax = 10
    cs = ax.pcolormesh(px,py,hyp_seas_DO_dict['avg'], vmin=vmin, vmax=vmax, cmap=cmap)
    cbar = fig.colorbar(cs, location='right')
    cbar.ax.tick_params(labelsize=12)
    cbar.outline.set_visible(False)

    # add 10 km bar
    ax.plot([lon0_10k,lon1_10k],[lat0_10k,lat1_10k],color='k',linewidth=2)
    ax.text(lon0_10k-0.04,lat0_10k+0.01,'{} km'.format(x_dist_km),color='k',fontsize=12)

    # format figure
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    helper_functions.dar(ax)
    plt.xticks(rotation=30)
    plt.yticks(rotation=30)
                                    
    # Add colormap title
    ax.set_title('(b)', fontsize=12, loc='left', fontweight='bold')

    # Generate plot
    plt.tight_layout
    plt.show()

    return