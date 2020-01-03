dir_ecco = '/Volumes/ECCO/ecco.jpl.nasa.gov/drive/files/Version4/Release4/';
dir_grid = '/Volumes/ECCO/ecco.jpl.nasa.gov/drive/files/Version4/Release3/nctiles_grid/';
% dir_grid = '/Volumes/ECCO/ecco.jpl.nasa.gov/drive/files/Version4/Release4/nctiles_grid/ECCOv4r3_grid.nc';
dir_monthly = [dir_ecco,'nctiles_monthly/'];
dir_daily = [dir_ecco,'nctiles_daily/'];
dir_init = [dir_ecco,'input_init/'];

global mygrid, mygrid = [];
grid_load(dir_grid,5,'nctiles',0,1);

%% Monthly volume budget
year_start = 1992; year_end = 2017;
n_month = (year_end - year_start + 1) * 12;
time = [year_start*ones(n_month,1) [2:(n_month+1)]' [1*ones(n_month-1,1);0.5]];
time = 24 * (datenum(time) - datenum(1992,1,1,12,0,0));
dt = diff([0,time']);
time_Units = 'hours_since_1992-1-1-12:00:00';
secPerHour = 3600;

dz_matF = mk3D(mygrid.DRF,mygrid.hFacC);
dz_mat = dz_matF.*mygrid.hFacC;
RAC_mat = mk3D(mygrid.RAC,mygrid.hFacC); % hFacC: mask of tiles, RAC: should be r coordinates
volume_grid = mygrid.mskC.*RAC_mat.*dz_mat; % volume
n_level = numel(mygrid.RC);

% load 2d thermal flux
f_id = fopen([dir_init,'geothermalFlux.bin'],'r','b');
geo_flux_2d = fread(f_id,'float32'); fclose(f_id);
geo_flux_2d = convert2gcmfaces(reshape(geo_flux_2d,[90,1170]));

% create 3d version
mask_C = mygrid.mskC;
mask_C(isnan(mask_C)) = 0;
mask_C_p1 = mask_C;
mask_C_p1(:,:,n_level+1) = 0;
mask_C_p1(:,:,1) = [];
mask_B = mask_C - mask_C_p1;
geo_flux_3d = mk3D(geo_flux_2d,mask_C) .* mask_B .* mask_C;
%% Loop through time step
for in_time = 1:n_month
    disp(num2str(in_time));
    % 2d
    
    for var = {'SSH','ETAN','oceFWflx'}
        str_cmd = [var,' = read_nctiles(strcat(dir_monthly,"',var,'"),"',var,'",in_time);'];
        eval(str_cmd);