"""
Main script to process data and generate figures for
studying drivers of low oxygen in Puget Sound terminal inlets.

Aurora Leeson
July 2025
"""

import pandas as pd
import xarray as xr
import matplotlib.pylab as plt
import pickle

# import helper functions
import helper_functions
import get_monthly_means
import hypoxic_volume
import budget_error
import budget_barchart
import plot_monthly_means
import dodeep_hypvol_timeseries
import net_decrease_boxplots
import multiple_regression

# reload to make editing easier
from importlib import reload
reload(helper_functions)
reload(get_monthly_means)
reload(hypoxic_volume)
reload(budget_error)
reload(budget_barchart)
reload(plot_monthly_means)
reload(dodeep_hypvol_timeseries)
reload(net_decrease_boxplots)
reload(multiple_regression)

# # DELETE THESE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# import QinDOin_correl_consumption
# import TEST1LAYER_net_decrease_boxplots
# import d_dt_DO_timeseries
# import TEST_STATISTICS

plt.close('all')

##########################################################
##                    Read in data                      ##
##########################################################

print('Reading data...')

# NOTE: data in these deeplay_dict and shallowlay_dict
# are tidally-averaged daily time series
# in units of kmol O2 per second
# Values have been passed through a 71-hour lowpass Godin filter
# (Thomson & Emery, 2014)

# terminal inlet deep layer values
with open('../data/deeplay_dict.pickle', 'rb') as handle:
    deeplay_dict = pickle.load(handle)

# terminal inlet shallow layer values
with open('../data/shallowlay_dict.pickle', 'rb') as handle:
    shallowlay_dict = pickle.load(handle)

# terminal inlet dimensions
with open('../data/dimensions_dict.pickle', 'rb') as handle:
    dimensions_dict = pickle.load(handle)

# terminal inlet DO concentrations [mg/L]
with open('../data/DOconcen_dict.pickle', 'rb') as handle:
    DOconcen_dict = pickle.load(handle)

# get inlet names
inlets = list(deeplay_dict.keys())

# list of hypoxic inlets
hyp_inlets = ['penn','case','holmes','portsusan','lynchcove','dabob']

# LiveOcean grid (cas7 version)
grid_ds = xr.open_dataset('../data/LO_cas7_grid.nc')

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
##             Hypoxic volume time series               ## 
##########################################################

hypoxic_volume.hypoxic_volume(grid_ds)

##########################################################
##               Deep Budget Error Analysis             ##
##########################################################

# calculate and print error of budget
# expressed as a % of QinDOin and biological consumption
budget_error.budget_error(inlets,shallowlay_dict,deeplay_dict,
                          dimensions_dict,kmolm3sec_to_mgLday)

##########################################################
##                  Budget Bar Charts                   ##
##########################################################

budget_barchart.budget_barchart(inlets,shallowlay_dict,deeplay_dict,
                    dates_local_hrly,dates_local_daily,hyp_inlets,
                    minday,maxday,kmolm3sec_to_mgLday)

#########################################################
##Plot monthly mean DOdeep, DOin, Tflush, and % hyp vol##
#########################################################
        
plot_monthly_means.plot_monthly_means(MONTHLYmean_DOdeep,
                                        MONTHLYmean_DOin,
                                        MONTHLYmean_Tflush,
                                        MONTHLYmean_perchyp,
                                        df_MONTHLYmean_DOdeep,
                                        df_MONTHLYmean_DOin,
                                        df_MONTHLYmean_Tflush)

##########################################################
##   Mean DOdeep vs % hyp vol and  DOdeep time series   ## 
##########################################################

dodeep_hypvol_timeseries.dodeep_hypvol_timeseries(MONTHLYmean_DOdeep,
                                                MONTHLYmean_perchyp,
                                                DOconcen_dict,
                                                dates_local_daily,
                                                dates_local_hrly,
                                                inlets,minday,maxday)

##########################################################
##        Net decrease (Jun 15 to Aug 15) boxplots      ## 
##########################################################

net_decrease_boxplots.net_decrease_boxplots(dimensions_dict,deeplay_dict,
                                            minday,maxday)

##########################################################
##                 Multiple regression                  ## 
##########################################################

multiple_regression.multiple_regression(MONTHLYmean_DOdeep,
                                        MONTHLYmean_DOin,
                                        MONTHLYmean_Tflush)







# ##########################################################
# ##   TESTING WHY IS D/DT(DO) THE SAME IN ALL INLETS?    ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ##########################################################
# # DELETE THESE

# TEST1LAYER_net_decrease_boxplots.net_decrease_boxplots(dimensions_dict,deeplay_dict,
#                                                         shallowlay_dict,minday,maxday)


# d_dt_DO_timeseries.d_dt_DO_timeseries(DOconcen_dict,
#                         dates_local_daily,
#                         dates_local_hrly,
#                         inlets,minday,maxday,
#                         dimensions_dict,deeplay_dict,
#                         shallowlay_dict,)

# TEST_STATISTICS.test_statistics(inlets,shallowlay_dict,deeplay_dict,
#                     dates_local_hrly,dates_local_daily,hyp_inlets,
#                     minday,maxday,kmolm3sec_to_mgLday)

# ########################################################
# # Correlation of Cons & QinDOin (mid-jun to mid-aug) ##
# ########################################################

# QinDOin_correl_consumption.correl(inlets,deeplay_dict,minday,maxday,kmolm3sec_to_mgLday)