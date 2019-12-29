import numpy as np
import xarray as xr
import sys
import matplotlib.pyplot as plt
import ecco_v4_py as ecco

## define a high-level directory for ECCO fields
dir_base = '/Volumes/ECCO/ecco.jpl.nasa.gov/drive/files/Version4'
dir_ECCO = dir_base +'/Release4'
dir_grid = dir_ECCO +'/nctiles_grid/'

## load the grid
grid = xr.open_dataset(dir_grid + '/ECCOv4r3_grid.nc')
# grid = ecco.load_ecco_grid_nc(dir_grid, 'ECCOv4r3_grid.nc')

# Method 1: using xarray.open_dataset
dir_ETAN = dir_ECCO + '/nctiles_monthly/ETAN/'
dataset_ETAN = xr.open_dataset(dir_ETAN + '/1992/ETAN_1992_01.nc')

## plot
ETAN = dataset_ETAN.ETAN.isel(time=0)
## mask to nan where hFacC(k=0) = 0
ETAN = ETAN.where(grid.hFacC.isel(k=0))

fig = plt.figure(figsize=(16,7))
ecco.plot_proj_to_latlon_grid(grid.XC, grid.YC, ETAN, show_colorbar=True, cmin = -1.5, cmax = 1.5)
plt.title('SSH (m)')

# Method 2, using load_ecco_var_from_years_nc
dataset_ETAN = ecco.load_ecco_var_from_years_nc(dir_ETAN,'ETAN',years_to_load=1992).load()
dataset_ETAN.attrs = []
dataset_ETAN

# Method 3, using dask
## Example 1a: Load two years monthly-mean ECCO fields into memory without Dask
dir_mon_mean = dir_ECCO + '/nctiles_monthly/'
import time
t_0 = time.time()

large_subset_no_dask = \
    ecco.recursive_load_ecco_var_from_years_nc(dir_mon_mean, \
                                               vars_to_load= ['THETA','ETAN'], \
                                               years_to_load=[2010, 2011], \
                                               less_output=True)
delta_t_2_yrs_no_dask = time.time() - t_0
print(delta_t_2_yrs_no_dask)

## Example 1b: Load two years of monthly-mean ECCO fields into memory using Dask
t_0 = time.time()
large_subset_with_dask = \
    ecco.recursive_load_ecco_var_from_years_nc(dir_mon_mean, \
                                               vars_to_load= ['THETA','ETAN'], \
                                               years_to_load=[2010, 2011])
delta_t_2_yrs_with_dask = time.time() - t_0
print(delta_t_2_yrs_with_dask)
## Example 2: Load two years of monthly-mean ECCO fields into memory with Dask
t_0 = time.time()
large_subset_with_dask1 = \
    ecco.recursive_load_ecco_var_from_years_nc(dir_mon_mean, \
                                               vars_to_load= ['THETA','ETAN'], \
                                               years_to_load=[2010, 2011],\
                                               dask_chunk=True)
delta_t_2_yrs_with_dask1 = time.time() - t_0
print(delta_t_2_yrs_with_dask1)