import warnings
warnings.filterwarnings('ignore')
import numpy as np
import xarray as xr
import ecco_v4_py as ecco
from datetime import date
import matplotlib.pyplot as plt
import cmocean

rho_const = 1029 # density, kg/m^3

map_dx = 0.2 # lon resolution
map_dy = 0.2 # lat resolution

dir_ecco = '/Volumes/ECCO/ecco.jpl.nasa.gov/drive/files/Version4/Release4/'
dir_monthly = dir_ecco + '/nctiles_monthly/'
dir_monthly_snapshots = dir_ecco + '/nctiles_monthly_snapshots'
dir_daily = dir_ecco + '/nctiles_daily/'
dir_grid = dir_ecco + '/nctiles_grid/'


fname_grid = 'ECCOv4r3_grid.nc'

## read grid
ds_grid = ecco.load_ecco_grid_nc(dir_grid,fname_grid, dask_chunk=True)

year_start = 1992
year_end = 2017

# load one extra year worth of snapshots
ds_ecco_monthly_snaps = \
    ecco.recursive_load_ecco_var_from_years_nc(dir_monthly_snapshots,\
                                               vars_to_load=['ETAN'],\
                                               years_to_load='all',\
                                               dask_chunk=True)
date_last_record = \
    ecco.extract_yyyy_mm_dd_hh_mm_ss_from_datetime64(ds_ecco_monthly_snaps.time[-1].values)

# load monthly mean data
ds_ecco_monthly_mean = \
    ecco.recursive_load_ecco_var_from_years_nc(dir_monthly,\
                                               vars_to_load=['ETAN',\
                                                             'oceFWflx',\
                                                             'UVELMASS',\
                                                             'VVELMASS',\
                                                             'WVELMASS'],\
                                               years_to_load='all',\
                                               dask_chunk=True)

grid_ecco_xgcm = ecco.get_llc_grid(ds_grid)

# calculte total tendency
num_months = len(ds_ecco_monthly_mean.time)
G_total_tendency_month = ds_ecco_monthly_mean.ETAN.isel(time=range(1,num_months)).values - ds_ecco_monthly_mean.ETAN.isel(time=range(0,num_months-1)).values
tmp = ds_ecco_monthly_mean.ETAN.isel(time=range(0,num_months-1)).copy(deep=True)
tmp.name = 'G_total_tendency_month'
tmp.values = G_total_tendency_month
G_total_tendency_month = tmp

t = []
for year in range(year_start,year_end+1):
    for month in range(1,13):
        t = np.append(t,date(year,month,1).toordinal())

seconds_per_month = 3600 * 24 * (t[1:]-t[0:len(t)-1])
seconds_per_month = xr.DataArray(seconds_per_month,\
                                 coords={'time': G_total_tendency_month.time.values},\
                                 dims='time')

G_total_tendency = G_total_tendency_month / seconds_per_month
G_total_tendency_mean = G_total_tendency.mean('time')
month_length_weights = seconds_per_month / seconds_per_month.sum()
G_total_tendency_mean = (G_total_tendency * month_length_weights).sum('time')

plt.figure(figsize=(20,8))
ecco.plot_proj_to_latlon_grid(ds_grid.XC, ds_grid.YC,\
                              G_total_tendency_mean,\
                              show_colorbar=True,\
                              cmin=-1e-9, cmax=1e-9, \
                              cmap=cmocean.cm.balance, user_lon_0=-67,\
                              dx=map_dx,dy=map_dy)
plt.title('Average $\partial \eta / \partial t$ [m/s]', fontsize=20);

ETAN_delta_method_2 = ds_ecco_monthly_mean.ETAN.isel(time=-1).values -  ds_ecco_monthly_mean.ETAN.isel(time=0).values
plt.figure(figsize=(20,8))
ecco.plot_proj_to_latlon_grid(ds_grid.XC, ds_grid.YC,\
                              ETAN_delta_method_2,\
                              show_colorbar=True,\
                              cmin=-1e-9, cmax=1e-9, \
                              cmap=cmocean.cm.balance, user_lon_0=-67,\
                              dx=map_dx,dy=map_dy)
plt.title('Actual $\Delta \eta$ [m]', fontsize=20);