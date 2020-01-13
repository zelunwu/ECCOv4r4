import numpy as np
import xarray as xr
import sys
import matplotlib.pyplot as plt
import json

import warnings
warnings.filterwarnings('ignore')
import ecco_v4_py as ecco

dir_base = '/Volumes/ECCO/ecco.jpl.nasa.gov/drive/files/Version4'
dir_ECCO = dir_base +'/Release4/'
dir_grid = dir_ECCO +'/nctiles_grid/'

## load the grid
# grid = xr.open_dataset(dir_grid + '/ECCOv4r3_grid.nc')
grid = ecco.load_ecco_grid_nc(dir_grid, 'ECCOv4r3_grid.nc')

dir_day_mean = dir_ECCO + '/nctiles_daily/'

ecco_vars = \
    ecco.recursive_load_ecco_var_from_years_nc(dir_day_mean, \
                                               vars_to_load=['SSH'], \
                                               years_to_load=2010,\
                                               dask_chunk=True)

ecco_ds = xr.merge([grid,ecco_vars])

# ssh_arr = ecco_ds.SSH.values
#
# fig = plt.figure(figsize=(8,6.5))
# plt.imshow(ssh_arr[0,2,:,:],origin='lower')
# fig = plt.figure(figsize=(8,6.5))
# plt.imshow(ecco_ds.SSH[0,2,:,:],origin='lower')
#
# ssh_masked = ecco_ds.SSH.where(grid.hFacC > 0, np.nan)
# ssh_masked[0,2,:,:,0].plot(cmap='jet')
#
# xc = ecco_ds.XC
# time = ecco_ds.time

for n in range(0,13):
    n_tile = n
    fig = plt.figure(figsize=(14,6))
    ecco.plot_proj_to_latlon_grid(ecco_ds.XC.isel(tile=n_tile),ecco_ds.YC.isel(tile=n_tile), \
                                  ecco_ds.SSH.isel(time = 0,tile=n_tile),\
                                  dx = 0.5, dy=0.5,user_lon_0=0,\
                                  show_colorbar=True,\
                                  show_grid_labels=True,\
                                  show_grid_lines=True)
    plt.title('tile = ' +str(n))

ssh_da = ecco_ds.SSH
lat_bounds = np.logical_and(ssh_da.YC  > -20, ssh_da.YC < 60)
lon_bounds = np.logical_and(ssh_da.XC  > -50, ssh_da.XC < 10)
lat_lon_bounds = np.logical_and(lat_bounds, lon_bounds)
ssh_da_subset_space = ssh_da.where(lat_lon_bounds, np.nan)
