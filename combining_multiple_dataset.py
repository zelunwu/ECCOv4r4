# Show how to combine multiple ECCO v4 state estimate Datasets after loading.
import xarray as xr
# import sys
import matplotlib.pyplot as plt
import ecco_v4_py as ecco

## define a high-level directory for ECCO fields
dir_base = '/Volumes/ECCO/ecco.jpl.nasa.gov/drive/files/Version4'
dir_ECCO = dir_base +'/Release4'
dir_grid = dir_ECCO +'/nctiles_grid/'

## load the grid
# grid = xr.open_dataset(dir_grid + '/ECCOv4r3_grid.nc')
grid = ecco.load_ecco_grid_nc(dir_grid, 'ECCOv4r3_grid.nc')

# load SSH
dir_ETAN = dir_ECCO + '/nctiles_monthly/ETAN/'
dataset_ETAN = ecco.load_ecco_var_from_years_nc(dir_ETAN,'ETAN',years_to_load=1992,dask_chunk=True)
# load U point ADVx_TH
dir_ADVx_TH = dir_ECCO + '/nctiles_monthly/ADVx_TH/'
dataset_ADVx_TH = ecco.load_ecco_var_from_years_nc(dir_ADVx_TH,'ADVx_TH',years_to_load=1992,dask_chunk=True)
# load V point ADVy_TH
dir_ADVy_TH = dir_ECCO + '/nctiles_monthly/ADVy_TH/'
dataset_ADVy_TH = ecco.load_ecco_var_from_years_nc(dir_ADVy_TH,'ADVy_TH',years_to_load=1992,dask_chunk=True)

# Merging multiple Datasets from state estimate variables
dataset_ecco = xr.merge([dataset_ETAN,dataset_ADVx_TH,dataset_ADVy_TH])
dataset_ecco_with_grid = xr.merge([dataset_ecco,grid])