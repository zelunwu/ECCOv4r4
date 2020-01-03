import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import ecco_v4_py as ecco
import  socket
hostName = socket.gethostname()

## define a high-level directory for ECCO fields
dir_base = '/Volumes/ECCO/ecco.jpl.nasa.gov/drive/files/Version4'
dir_ECCO = dir_base +'/Release4'
if hostName[-1] == '3':
    dir_ECCO = '/Volumes/15394457571/ecco_eg'
dir_grid = dir_ECCO + '/nctiles_grid'
dir_daily = dir_ECCO + '/nctiles_daily/'

## load the grid
# grid = xr.open_dataset(dir_grid + '/ECCOv4r3_grid.nc')
grid = ecco.load_ecco_grid_nc(dir_grid, 'ECCOv4r3_grid.nc')
ecco_ds = ecco.recursive_load_ecco_var_from_years_nc(dir_daily,\
                                                     vars_to_load=['SSH'],\
                                                     years_to_load=1992,
                                                     dask_chunk=True)
ecco_ds = xr.merge([grid,ecco_ds])

# Area_global
ocean_mask = np.ceil(ecco_ds.hFacC)
ocean_mask = ocean_mask.where(ocean_mask==1, np.nan)

# plt.figure(figsize=(12,9),dpi=90)
ax = ecco.plot_tiles(ocean_mask.sel(k=0),layout='latlon', rotate_to_latlon=True,fig_size=6)
plt.suptitle('test')

area_sum = np.nansum(ocean_mask.isel(k=0)*ecco_ds.rA)
ssh_global_mean = (ecco_ds.SSH*ecco_ds.rA*ocean_mask.isel(k=0)).sum(dim=['i','j','tile'])/area_sum
ssh_global_mean = ssh_global_mean - ssh_global_mean.mean(dim='time')
ssh_global_mean.['units']='m'
ssh_global_mean.plot()
