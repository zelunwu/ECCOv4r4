import ecco_v4_py as ecco
# import matplotlib.pyplot as plt
import numpy as np
# import numpy.ma as ma
import xarray as xr

# def isleap(year):
#     if (year % 4) == 0:
#        if (year % 100) == 0:
#            if (year % 400) == 0:
#                leap = True
#            else:
#                leap = False
#        else:
#            leap = True
#     else:
#        leap = False
#     return leap

dir_ecco = '/Volumes/ECCO/ecco.jpl.nasa.gov/drive/files/Version4/Release4/'
dir_monthly = dir_ecco + 'nctiles_monthly/'
dir_grid = dir_ecco + '/nctiles_grid/'
dir_daily = dir_ecco + 'nctiles_daily/'

dir_output = '/Volumes/ECCO/ecco.jpl.nasa.gov/drive/files/Version4/Release4/interp_daily/SST/chinasea/'

ds_grid = ecco.load_ecco_grid_nc(dir_grid,'ECCOv4r3_grid.nc')
mask_C = ds_grid.hFacC.isel(k=0)
mask_C = mask_C.where(mask_C>0, np.nan)

# resample
dx = 0.25; dy = 0.25
min_lat = 0; max_lat = 50; min_lon = 100; max_lon = 140;
n_lon = len(np.arange(min_lon+dx,max_lon-dx,dx))+1; n_lat = len(np.arange(min_lat+dy,max_lat-dy,dy))+1;
#%% load data

for year in range(1992,2018):
    ds_data = ecco.load_ecco_var_from_years_nc(dir_daily + '/THETA/'+ str(year), var_to_load='THETA', years_to_load=[year],
                                               k_subset=[0])
    n_time = len(ds_data.time)
    sst = np.full([n_time, n_lat, n_lon], np.nan)

    for in_time in range(0,n_time):
        print(year,in_time)
        lon,lat,sst[in_time,:,:] = ecco.resample_to_latlon(ds_data.XC,ds_data.YC,ds_data.THETA.isel(time=in_time, k=0)*mask_C,\
                                              min_lat+dy, max_lat-dy, dy,\
                                              min_lon+dx, max_lon-dx, dx,\
                                              mapping_method='nearest_neighbor',fill_value=np.nan)
    ds_sst = xr.Dataset({})
    ds_sst.coords['lon'] = lon[0,:]
    ds_sst.coords['lat'] = lat[:,0]
    ds_sst.coords['time'] = ds_data.time
    ds_sst['sst'] = (('time','lat','lon'),sst)
    # name
    ds_sst.lon['long_name'] = 'longitude (degree)'
    ds_sst.lat['long_name'] = 'latitude (degree)'
    ds_sst.sst['long_name'] = 'sea surface temperature (degree C)'
    # output
    file_output = dir_output + 'chinasea_sst_' + str(year) + '.nc'
    print('Saving '+file_output + '...')
    ds_sst.to_netcdf(path=file_output)
    print('Finish ' + file_output + '.')

