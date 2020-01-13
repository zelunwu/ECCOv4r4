import ecco_v4_py as ecco 
import xarray as xr 

dir_ecco = '/Volumes/ECCO/ecco.jpl.nasa.gov/drive/files/Version4/Release4'
dir_daily = dir_ecco + '/nctiles_daily/'
dir_grid = dir_ecco + '/nctiles_grid/'

ds_grid = ecco.load_ecco_grid_nc(dir_grid,'ECCOv4r3_grid.nc')

vars_to_load = []