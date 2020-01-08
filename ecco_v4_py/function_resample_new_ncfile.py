def resample_to_latlon_and_save_to_nc(dir_data,dir_grid,dir_output,\
                                      var_to_load,years_to_load, k_sub = [0], \
                                      min_lat=-90, max_lat=90,min_lon=0, max_lon=360, dx = 0.25, dy = 0.25,\
                                      mapping_method='nearest_neighbor',fill_value=np.nan):
    import ecco_v4_py as ecco
    import numpy as np
    import xarray as xr

    ds_grid = ecco.load_ecco_grid_nc(dir_grid)
    mask_C = ds_grid.hFacC.isel(k=k_sub)
    mask_C = mask_C.where(mask_C>0, np.nan)

    # resample
    n_lon = len(np.arange(min_lon+dx,max_lon-dx,dx))+1; n_lat = len(np.arange(min_lat+dy,max_lat-dy,dy))+1;
    #%% load data

    for year in list(years_to_load):
        ds_data = ecco.load_ecco_var_from_years_nc(dir_data + '/' + str(year), var_to_load=var_to_load, years_to_load=[year],
                                                   k_subset=[0])
        n_time = len(ds_data.time)
        sst = np.full([n_time, n_lat, n_lon], np.nan)

        for in_time in range(0,n_time):
            print(year,in_time+1)
            lon,lat,sst[in_time,:,:] = ecco.resample_to_latlon(ds_data.XC,ds_data.YC,ds_data.THETA.isel(time=in_time, k=0)*mask_C,\
                                                  min_lat+dy, max_lat-dy, dy,\
                                                  min_lon+dx, max_lon-dx, dx,\
                                                  mapping_method=mapping_method,fill_value=fill_value)
        ds_data_save = xr.Dataset({})
        ds_data_save.coords['lon'] = lon[0,:]
        ds_data_save.coords['lat'] = lat[:,0]
        ds_data_save.coords['time'] = ds_data.time
        ds_data_save[var_to_load] = (('time','lat','lon'),sst)
        # name
        ds_data_save.lon['long_name'] = 'longitude (degree)'
        ds_data_save.lat['long_name'] = 'latitude (degree)'
        ds_data_save[var_to_load]['long_name'] = var_to_load
        # output
        file_output = dir_output + var_to_load + str(year) + '.nc'
        print('Saving '+file_output + '...')
        ds_data_save.to_netcdf(path=file_output)
        print('Finish ' + file_output + '.')