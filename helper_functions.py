"""
This module contains utility functions for filtering.

Parker MacCready
"""

import numpy as np
import pytz

def lowpass(data, f='hanning', n=40, nanpad=True):
    """
    A replacement for almost all previous filter code.
    f = 'hanning' (default) or 'godin'
    
    Input: ND numpy array, any number of dimensions, with time on axis 0.
    
    Output: Array of the same size, filtered with Hanning window of length n,
        or the Godin filter (hourly data only) padded with nan's.
    """
    if n == 1:
        return data
    else:
        if f == 'hanning':
            filt = hanning_shape(n=n)
        elif f == 'godin':
            filt = godin_shape()
        else:
            print('ERROR in filt_general(): unsupported filter ' + f)
            filt = np.nan
        npad = np.floor(len(filt)/2).astype(int)
        sh = data.shape
        df = data.flatten('F')
        dfs = np.convolve(df, filt, mode = 'same')
        smooth = dfs.reshape(sh, order='F')
        # note that the indexing below defaults to being on axis 0,
        # and correctly broadcasts without having to mention the other axes
        if nanpad:
            smooth[:npad] = np.nan
            smooth[-npad:] = np.nan
        else:
            smooth[:npad] = data[:npad]
            smooth[-npad:] = data[-npad:]
        return smooth
    
def godin_shape():
    """
    Based on matlab code of 4/8/2013  Parker MacCready
    Returns a 71 element numpy array that is the weights
    for the Godin 24-24-25 tildal averaging filter. This is the shape given in
    Emery and Thomson (1997) Eqn. (5.10.37)
    ** use ONLY with hourly data! **
    """
    k = np.arange(12)
    filt = np.NaN * np.ones(71)
    filt[35:47] = (0.5/(24*24*25))*(1200-(12-k)*(13-k)-(12+k)*(13+k))
    k = np.arange(12,36)
    filt[47:71] = (0.5/(24*24*25))*(36-k)*(37-k)
    filt[:35] = filt[:35:-1]
    return filt

def hanning_shape(n=40):
    """
    Returns a Hanning window of the specified length.
    """
    ff = np.cos(np.linspace(-np.pi,np.pi,n+2))[1:-1]
    filt = (1 + ff)/2
    filt = filt / filt.sum()
    return filt

def get_dt_local(dt, tzl='US/Pacific'):
    # take a model datetime (assumed to be UTC) and return local datetime
    tz_utc = pytz.timezone('UTC')
    tz_local = pytz.timezone(tzl)
    dt_utc = dt.replace(tzinfo=tz_utc)
    dt_local = dt_utc.astimezone(tz_local)
    return dt_local

def get_plon_plat(lon, lat):
    """
    This takes the 2-D lon and lat grids (ndarrays) and returns extended
    "psi" grids that are suitable for plotting using pcolormesh
    for any field on the original grid.
    NOTE: It checks to make sure the original grid is plaid.
    
    You would pass it G['lon_rho'] and G['lat_rho'], for
    a field on the rho grid.
    """
    # input checking
    Lon = lon[0,:]
    Lat = lat[:,0]
    if (Lon - lon[-1,:]).sum() != 0:
        print('Error from get_plon_plat: lon grid not plaid')
        sys.exit()
    if (Lat - lat[:,-1]).sum() != 0:
        print('Error from get_plon_plat: lat grid not plaid')
        sys.exit()
    plon = np.ones(len(Lon) + 1)
    plat = np.ones(len(Lat) + 1)
    dx2 = np.diff(Lon)/2
    dy2 = np.diff(Lat)/2
    Plon = np.concatenate(((Lon[0]-dx2[0]).reshape((1,)), Lon[:-1]+dx2, (Lon[-1]+dx2[-1]).reshape((1,))))
    Plat = np.concatenate(((Lat[0]-dy2[0]).reshape((1,)), Lat[:-1]+dy2, (Lat[-1]+dy2[-1]).reshape((1,))))
    plon, plat = np.meshgrid(Plon, Plat)
    return plon, plat

def dar(ax):
    """
    Fixes the plot aspect ratio to be locally Cartesian.
    """
    yl = ax.get_ylim()
    yav = (yl[0] + yl[1])/2
    ax.set_aspect(1/np.cos(np.pi*yav/180))