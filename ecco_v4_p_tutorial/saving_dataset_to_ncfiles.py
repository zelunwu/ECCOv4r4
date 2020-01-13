import xarray as xr
import json
import ecco_v4_py as ecco

dir_base = '/Volumes/ECCO/ecco.jpl.nasa.gov/drive/files/Version4'
dir_ECCO = dir_base +'/Release4'
dir_grid = dir_ECCO +'/nctiles_grid/'

## load the grid
# grid = xr.open_dataset(dir_grid + '/ECCOv4r3_grid.nc')
grid = ecco.load_ecco_grid_nc(dir_grid, 'ECCOv4r3_grid.nc')

# load a single tile, monthly average THETA for March 2010 for model tile 2
dir_theta = dir_ECCO + '/nctiles_monthly/THETA/2010/'
fname = 'THETA_2010_03.nc'
dataset_theta = xr.open_dataset(dir_theta+fname).load()

fname_save = './test_output.nc'
print('saving to ',fname_save)
dataset_theta.to_netcdf(path=fname_save)
print('finished saving')

# Saving a new custom Dataset to NetCDF
dataset_ecco.to_netcdf(path=fname_save)

# Saving the results of calculations
dir_ETAN = dir_ECCO + '/nctiles_monthly/ETAN/'
dataset_ETAN = xr.open_dataset(dir_ETAN + '/1992/ETAN_1992_01.nc')

SSH_THETA_201003 = \
    ecco.recursive_load_ecco_var_from_years_nc(dir_ECCO+'/nctiles_monthly/', \
                                              ['SSH', 'THETA'], \
                                              tiles_to_load = [0,1,2],
                                              years_to_load = 2010)

SSH_sq = SSH_THETA_201003.SSH * SSH_THETA_201003.SSH
SSH_sq.name = 'SSH^2'
SSH_sq.attrs['long_name'] = 'Square of Surface Height Anomaly'
SSH_sq.attrs['units'] = 'm^2'
SSH_sq.to_netcdf(path='./test.nc')

