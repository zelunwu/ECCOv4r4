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
# grid2 = ecco.load_ecco_grid_nc(grid_dir, 'ECCOv4r3_grid.nc')
## load the grid
grid = xr.open_dataset(dir_grid + '/ECCOv4r3_grid.nc')

# plot hFac
ecco.plot_tiles(grid.hFacC.sel(k=0), show_colorbar=True, cmap='gray',figsize=(10, 9));

ecco.plot_tiles(grid.rA, show_colorbar=True);

## Method 2: Loading the model grid parameters using load_ecco_grid_nc
grid_subset = ecco.load_ecco_grid_nc(dir_grid, 'ECCOv4r3_grid.nc', tiles_to_load = 'all', k=0)
ecco.plot_tiles(grid_subset.hFacC.sel(k=0), show_colorbar=True, cmap='gray');