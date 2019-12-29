import numpy as np
import xarray as xr
import sys
import matplotlib.pyplot as plt
import ecco_v4_py as ecco
import warnings
warnings.filterwarnings('ignore')

## define a high-level directory for ECCO fields
# dir_base = '/Volumes/ECCO/ecco.jpl.nasa.gov/drive/files/Version4'
# dir_ECCO = dir_base +'/Release4'
dir_ECCO = '/Volumes/15394457571/ecco_eg'
dir_grid = dir_ECCO +'/nctiles_grid/'
dir_daymean = dir_ECCO + '/nctiles_daily/'

## load the grid
# grid = xr.open_dataset(dir_grid + '/ECCOv4r3_grid.nc')
grid = ecco.load_ecco_grid_nc(dir_grid, 'ECCOv4r3_grid.nc')

ecco_ds = ecco.recursive_load_ecco_var_from_years_nc(dir_daymean,\
                                                       vars_to_load=['SSH'],\
                                                       years_to_load=1992,\
                                                       dask_chunk=True)
ecco_ds = xr.merge([ecco_ds,grid])
# plot a single tile with imshow
plt.figure(figsize=(6,5),dpi=90)
plt.imshow(ecco_ds.SSH.isel(tile=2,time=0).where(grid.hFacC.isel(tile=2,k=0) !=0, np.nan))
plt.colorbar()
plt.xlabel('x -->')
plt.ylabel('y -->')
plt.title('SSH (m')
# plot a single tile with pcolor and contourf
tile_num=2
lons = grid.XC.sel(tile=tile_num)
lats = grid.YC.sel(tile=tile_num)
tile_to_plot = ecco_ds.SSH.isel(tile=tile_num, time=0)
tile_to_plot= tile_to_plot.where(ecco_ds.hFacC.isel(tile=tile_num,k=0) !=0, np.nan)

plt.figure(figsize=(10,10))
plt.subplot(121)
plt.pcolor(lons, lats, tile_to_plot,vmin=-1,vmax=1,cmap='jet')
plt.colorbar()

plt.subplot(122)
plt.contourf(lons, lats, tile_to_plot, np.linspace(-1,1, 20,endpoint=True),\
             vmin=-1, vmax=1, cmap='jet',shading='interp')
plt.colorbar()