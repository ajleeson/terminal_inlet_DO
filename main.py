"""
Main script to process data and generate figures for
studying drivers of low oxygen in Puget Sound terminal inlets.

Aurora Leeson
August 2025
"""

import pandas as pd
import xarray as xr
import matplotlib.pylab as plt
import pickle

# import helper functions
import helper_functions
import get_monthly_means
import budget_error
import figure_01
import figure_07
import figure_08
import figure_09
import figure_10
import figure_11
import figure_12
import multiple_regression

# reload to make editing easier
from importlib import reload
reload(helper_functions)
reload(get_monthly_means)
reload(budget_error)
reload(figure_01)
reload(figure_07)
reload(figure_08)
reload(figure_09)
reload(figure_10)
reload(figure_11)
reload(figure_12)
reload(multiple_regression)

plt.close('all')

##########################################################
##                    Read in data                      ##
##########################################################

print('Reading data...')


# LiveOcean grid (cas7 version)
grid_ds = xr.open_dataset('../DATA_terminal_inlet_DO/LO_cas7_grid.nc')

# Puget Sound sub-domain within LiveOcean
PSbox_ds = xr.open_dataset('../DATA_terminal_inlet_DO/PugetSound_gridsizes.nc')

# Puget Sound hypoxic volume time series
with open('../DATA_terminal_inlet_DO/PS_hypoxic_volume_dict.pickle', 'rb') as handle:
    hyp_vol_dict = pickle.load(handle)

# Number of days that each grid cell experiences bottom hypoxia per year
with open('../DATA_terminal_inlet_DO/days_with_bottom_hypoxia_dict.pickle', 'rb') as handle:
    hyp_days_dict = pickle.load(handle)

# Mean bottom DO concentration of each grid cell during hypoxic season
with open('../DATA_terminal_inlet_DO/mean_hypoxic_season_bottom_DO_dict.pickle', 'rb') as handle:
    hyp_seas_DO_dict = pickle.load(handle)

# NOTE: data in deeplay_dict and shallowlay_dict
# are tidally-averaged daily time series
# in units of kmol O2 per second
# Values have been passed through a 71-hour lowpass Godin filter
# (Thomson & Emery, 2014)

# terminal inlet deep layer values
with open('../DATA_terminal_inlet_DO/deeplay_dict.pickle', 'rb') as handle:
    deeplay_dict = pickle.load(handle)

# terminal inlet shallow layer values
with open('../DATA_terminal_inlet_DO/shallowlay_dict.pickle', 'rb') as handle:
    shallowlay_dict = pickle.load(handle)

# terminal inlet dimensions
with open('../DATA_terminal_inlet_DO/dimensions_dict.pickle', 'rb') as handle:
    dimensions_dict = pickle.load(handle)

# terminal inlet DO concentrations [mg/L]
with open('../DATA_terminal_inlet_DO/DOconcen_dict.pickle', 'rb') as handle:
    DOconcen_dict = pickle.load(handle)

# get inlet names
inlets = list(deeplay_dict.keys())

# list of hypoxic inlets
hyp_inlets = ['penn','case','holmes','portsusan','lynchcove','dabob']

##########################################################
##                 Key values                           ##
##########################################################

# convert from kmol O2 per m3 per second to mg/L per day
kmolm3sec_to_mgLday = 1000 * 32 * 60 * 60 * 24

# yearday of drawdown period (June 15 through August 15)
minday = 164
maxday = 225

##########################################################
##   Get dates for analysis (2017.01.02 to 2017.12.30)  ##
##########################################################

year = '2017'

# set up dates
startdate = year + '.01.01'
enddate = year + '.12.31'
enddate_hrly = str(int(year)+1)+'.01.01 00:00:00'

# create time_vector
dates_hrly = pd.date_range(start= startdate, end=enddate_hrly, freq= 'h')
dates_local_hrly = [helper_functions.get_dt_local(x) for x in dates_hrly]
# crop time vector (because we only have jan 2 - dec 30)
dates_daily = pd.date_range(start= startdate, end=enddate, freq= 'd')[2::]
dates_local_daily = [helper_functions.get_dt_local(x) for x in dates_daily]

##########################################################
##                 Get monthly means                    ## 
##########################################################

# MONTHLYmean_XXXX are arrays of monthly mean values
# for all inlets, compressed into a single array

# df_MONTHLY_mean_XXX are dataframes, where each column
# is an individual inlet. All columns contain monthly
# mean values corresponding to the inlet (ie., 12 rows)

[MONTHLYmean_DOdeep,
MONTHLYmean_DOin,
MONTHLYmean_Tflush,
MONTHLYmean_perchyp,
df_MONTHLYmean_DOdeep,
df_MONTHLYmean_DOin,
df_MONTHLYmean_Tflush,
df_MONTHLYmean_perchyp] = get_monthly_means.get_monthly_means(deeplay_dict,DOconcen_dict,
                                                                dimensions_dict,inlets)

##########################################################
##               Deep Budget Error Analysis             ##
##########################################################

# calculate and print error of budget
# expressed as a % of QinDOin and biological consumption
budget_error.budget_error(inlets,shallowlay_dict,deeplay_dict,
                          dimensions_dict,kmolm3sec_to_mgLday)

##########################################################
##                   Bathymetry map                     ## 
##########################################################

figure_01.model_bathy(grid_ds)

##########################################################
##             Hypoxic volume time series               ## 
##########################################################

figure_07.hypoxic_volume(grid_ds,hyp_vol_dict,PSbox_ds)

##########################################################
##              Map of Puget Sound hypoxia              ## 
##########################################################

figure_08.pugetsound_hyp_map(grid_ds,PSbox_ds,hyp_days_dict,
                             hyp_seas_DO_dict)

##########################################################
##   Mean DOdeep vs % hyp vol and  DOdeep time series   ## 
##########################################################

figure_09.dodeep_hypvol_timeseries(MONTHLYmean_DOdeep,
                                    MONTHLYmean_perchyp,
                                    DOconcen_dict,
                                    dates_local_daily,
                                    dates_local_hrly,
                                    inlets,minday,maxday)

##########################################################
##                  Budget Bar Charts                   ##
##########################################################

figure_10.budget_barchart(inlets,shallowlay_dict,deeplay_dict,
                    dates_local_hrly,dates_local_daily,hyp_inlets,
                    minday,maxday,kmolm3sec_to_mgLday)

##########################################################
##        Net decrease (Jun 15 to Aug 15) boxplots      ## 
##########################################################

figure_11.net_decrease_boxplots(dimensions_dict,deeplay_dict,
                                minday,maxday)

#########################################################
##Plot monthly mean DOdeep, DOin, Tflush, and % hyp vol##
#########################################################
        
figure_12.plot_monthly_means(MONTHLYmean_DOdeep,
                            MONTHLYmean_DOin,
                            MONTHLYmean_Tflush,
                            MONTHLYmean_perchyp,
                            df_MONTHLYmean_DOdeep,
                            df_MONTHLYmean_DOin,
                            df_MONTHLYmean_Tflush)

##########################################################
##                 Multiple regression                  ## 
##########################################################

multiple_regression.multiple_regression(MONTHLYmean_DOdeep,
                                        MONTHLYmean_DOin,
                                        MONTHLYmean_Tflush)