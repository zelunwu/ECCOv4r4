import warnings
warnings.filterwarnings('ignore')

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import ecco_v4_py as ecco
import  socket
hostName = socket.gethostname()

import ecco_v4_py.vector_calc as vc
import ecco_v4_py.scalar_calc as sc

## define a high-level directory for ECCO fields
dir_base = '/Volumes/ECCO/ecco.jpl.nasa.gov/drive/files/Version4'
dir_ECCO = dir_base +'/Release4'
if hostName[-1] == '3':
    dir_ECCO = '/Volumes/15394457571/ecco_eg'
dir_grid = dir_ECCO + '/nctiles_grid'
dir_daily = dir_ECCO + '/nctiles_daily/'
dir_monthly = dir_ECCO + '/nctiles_monthly/'

## load the grid
# grid = xr.open_dataset(dir_grid + '/ECCOv4r3_grid.nc')
grid = ecco.load_ecco_grid_nc(dir_grid, 'ECCOv4r3_grid.nc')
ecco_ds = ecco.recursive_load_ecco_var_from_years_nc(dir_monthly,\
                                                     vars_to_load=['ADVx_TH', 'ADVy_TH', \
                                                        'DFxE_TH', 'DFyE_TH'],\
                                                     years_to_load='all',\
                                                     dask_chunk=True)
ecco_ds = xr.merge([ecco_ds,grid])

grid = ecco.get_llc_grid(ecco_ds)
rapid_mask_W, rapid_mask_S = ecco.vector_calc.get_latitude_masks(lat_val=26,yc=ecco_ds.YC,grid=grid)
rapid_mask_C = ecco.scalar_calc.get_latitude_mask(lat_val=26,yc=ecco_ds.YC, grid=grid)


lat = 26
ones = xr.ones_like(ecco_ds.YC)
dome_mask_C = ones.where(ecco_ds.YC>=lat,0)

plt.figure(figsize=(12,6))
ecco.plot_proj_to_latlon_grid(ecco_ds.XC,ecco_ds.YC,dome_mask_C,\
                              projection_type='robin',cmap='jet',\
                              user_lon_0=0,show_colorbar=True,show_grid_lines=True,show_grid_labels=True)