import netCDF4 as nc
import ecco_v4_py as ecco
import numpy as np
import xarray as xr
import sys
import json
from matplotlib import pyplot as plt

## Set top-level file directory for the ECCO NetCDF files
## =================================================================
## define a high-level directory for ECCO fields
base_dir = '/Volumes/ECCO/ecco.jpl.nasa.gov/drive/files/Version4'
ECCO_dir = base_dir + '/Release4'

## LOAD NETCDF FILE
## ================
# grid
grid_dir = base_dir + ECCO_dir +'/nctiles_grid'
## load the grid
grid_dataset = xr.open_dataset(grid_dir + '/ECCOv4r3_grid.nc')


# directory containing the file
data_dir= ECCO_dir + '/nctiles_monthly/THETA'

# filename
fname = '/1992/THETA_1992_01.nc'

# load the file
theta_dataset = xr.open_dataset(data_dir + fname).load()
# figures
fig = plt.figure(figsize=(8, 6.5))
theta_dataset.THETA.isel(k=0,tile=2,time=0).plot(vmin=-2, vmax=25, cmap='jet')



# UVEL file
# Directory of the UVEL files
data_dir= ECCO_dir + '/nctiles_monthly/UVEL/'

fname = '/1992/UVEL_1992_01.nc'
uvel_dataset = xr.open_dataset(data_dir + fname).load()

fig = plt.figure(figsize=(8, 6.5))
ud_masked = uvel_dataset.UVEL.where(grid_dataset.hFacW > 0, np.nan)
ud_masked.isel(k=0,tile=1, time=0).plot(cmap='jet', vmin=-.2,vmax=.2)
# uvel_dataset.UVEL.isel(k=0,tile=1, time=0).plot(cmap='jet', vmin=-.2,vmax=.2)

# Directory of the VVEL files
data_dir= ECCO_dir + '/nctiles_monthly//VVEL'

fname = '/1992/VVEL_1992_01.nc'
vvel_dataset = xr.open_dataset(data_dir + fname).load()
vvel_dataset.attrs = []

fig=plt.figure(figsize=(8, 6.5))
vd_masked = vvel_dataset.VVEL.where(grid_dataset.hFacS > 0, np.nan)
vd_masked.isel(k=0,tile=1, time=0).plot(cmap='jet', vmin=-.2,vmax=.2)