import warnings
warnings.filterwarnings('ignore')

import matplotlib.pyplot as plt
import xarray as xr
import ecco_v4_py as ecco
dir_ecco = '/Volumes/ECCO/ecco.jpl.nasa.gov/drive/files/Version4/Release4/'
dir_monthly = dir_ecco + '/nctiles_monthly/'
dir_daily = dir_ecco + '/nctiles_daily/'
dir_grid = dir_ecco + '/nctiles_grid/'

fname_grid = 'ECCOv4r3_grid.nc'

## read grid
ds_grid = ecco.load_ecco_grid_nc(dir_grid,fname_grid, dask_chunk=True)
vars_ecco = ecco.recursive_load_ecco_var_from_years_nc(dir_monthly,\
                                                      vars_to_load=['ADVx_TH', 'ADVy_TH',\
                                                                    'DFxE_TH', 'DFyE_TH',\
                                                                    'ADVx_SLT','ADVy_SLT',\
                                                                    'DFxE_SLT','DFyE_SLT',\
                                                                    'UVELMASS','VVELMASS'],\
                                                      years_to_load='all',\
                                                      dask_chunk=True)
ds_ecco = xr.merge([vars_ecco,ds_grid])
del(ds_grid,vars_ecco)

lat = 26
ds_ecco['maskC'] = ds_ecco.hFacC
ds_ecco['maskW'] = ds_ecco.hFacW
ds_ecco['maskS'] = ds_ecco.hFacS

grid = ecco.get_llc_grid(ds_ecco)
rapid_maskW, rapid_maskS = ecco.vector_calc.get_latitude_masks(lat_val=lat, yc=ds_ecco.YC, grid=grid)
rapid_maskC = ecco.scalar_calc.get_latitude_mask(lat_val=lat, yc=ds_ecco.YC, grid=grid)

ones = xr.ones_like(ds_ecco.YC)
dome_maskC = ones.where(ds_ecco.YC>=lat, 0)

# plt.figure(figsize=(12,5))
# ecco.plot_proj_to_latlon_grid(dome_maskC.XC,dome_maskC.YC,dome_maskC,\
#                               projection_type='robin',\
#                               cmap='hot',
#                               user_lon_0=0,\
#                               show_colorbar=True, show_grid_labels=True, show_grid_lines=True)

maskC = ds_ecco.maskC
maskW = ds_ecco.maskW
maskS = ds_ecco.maskS

plt.figure(figsize=(12,5))
ecco.plot_tiles(dome_maskC+maskC.isel(k=0), cmap='viridis')
plt.show()


lat_maskW = grid.diff(dome_maskC,'X',boundary='fill')
lat_maskS = grid.diff(dome_maskC,'Y',boundary='fill')

atl_maskW = ecco.get_basin_mask(basin_name='atlExt', mask=maskW.isel(k=0))
atl_maskS = ecco.get_basin_mask(basin_name='atlExt', mask=maskS.isel(k=0))

plt.figure(figsize=(12,5))
ecco.plot_proj_to_latlon_grid(ds_ecco.XC, ds_ecco.YC, atl_maskW, projection_type='robin', cmap='viridis', user_lon_0=0)
plt.show()

mvt = ecco.calc_meridional_vol_trsp(ds_ecco,lat_vals=lat, basin_name='atlExt')
mht = ecco.calc_meridional_heat_trsp(ds_ecco,lat_vals=lat, basin_name='atlExt')
mst = ecco.calc_meridional_salt_trsp(ds_ecco,lat_vals=lat, basin_name='atlExt')

plt.figure(figsize=(12,5))
plt.plot(mvt.time,mvt.vol_trsp)
plt.ylabel('Meridional Volume Transport (Sv)')
plt.show()
plt.savefig('mvt_26N_Altlantic.png')

plt.figure(figsize=(12,5))
plt.plot(mht.time,mht.heat_trsp)
plt.ylabel('Meridional Heat Transport (PW)')
plt.show()
plt.savefig('mht_26N_Altlantic.png')

plt.figure(figsize=(12,5))
plt.plot(mst.time,mst.salt_trsp)
plt.ylabel('Meridional Salt Transport (psu*Sv)')
plt.show()
plt.savefig('mst_26N_Altlantic.png')