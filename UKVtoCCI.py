import numpy as np
from pyproj import Transformer
import xarray as xr
    
# input and output directories for meteorological data 
idir = '/exports/csce/datastore/geos/groups/boreal/SnowPatches/data/UKV_nc/2016_2021_FSM/'
odir = 'UKV_Cairngorms/'

# ancillary file defining snow cci grid for Cairngorms domain
CCI = xr.open_dataset(odir+'ancil.nc')
lons = CCI.lon[:]
lats = CCI.lat[:]

# transformation between UKV and snow cci grids
lons,lats = np.meshgrid(lons,lats)
cci_to_UKV = Transformer.from_crs('EPSG:4326','EPSG:27700')
x_cci,y_cci = cci_to_UKV.transform(lats,lons)
x_cci = xr.DataArray(x_cci,dims=('lat','lon'))
y_cci = xr.DataArray(y_cci,dims=('lat','lon'))

# interpolate UKV data onto snow cci grid
method = 'nearest'
for y in range(6):
    year = str(2016+y)
    for m in range(12):
        month = str(m+1)
        if m<9: month = '0'+month
        try:
            print(year,month)
            UKV = xr.open_dataset(idir+'LWdown'+year+month+'.nc')
            v = UKV.LWdown[:,:,38:] 
            v = v.interp({'x':x_cci,'y':y_cci},method=method)
            v.to_netcdf(odir+'LWdown'+year+month+'.nc')
            UKV = xr.open_dataset(idir+'PSurf'+year+month+'.nc')
            v = UKV.PSurf[:,:,38:] 
            v = v.interp({'x':x_cci,'y':y_cci},method=method)
            v.rename('Psurf')
            v.to_netcdf(odir+'Psurf'+year+month+'.nc')
            UKV = xr.open_dataset(idir+'Qair'+year+month+'.nc')
            v = UKV.Qair[:,:,38:] 
            v = v.interp({'x':x_cci,'y':y_cci},method=method)
            v.to_netcdf(odir+'Qair'+year+month+'.nc')
            UKV = xr.open_dataset(idir+'Rainf'+year+month+'.nc')
            v = UKV.Rainf[:,:,38:] 
            v = v.interp({'x':x_cci,'y':y_cci},method=method)
            v.to_netcdf(odir+'Rainf'+year+month+'.nc')
            UKV = xr.open_dataset(idir+'Snowf'+year+month+'.nc')
            v = UKV.Snowf[:,:,38:] 
            v = v.interp({'x':x_cci,'y':y_cci},method=method)
            v.to_netcdf(odir+'Snowf'+year+month+'.nc')
            UKV = xr.open_dataset(idir+'SWdown'+year+month+'.nc')
            v = UKV.SWdown[:,:,38:] 
            v = v.interp({'x':x_cci,'y':y_cci},method=method)
            v.to_netcdf(odir+'SWdown'+year+month+'.nc')
            UKV = xr.open_dataset(idir+'Tair'+year+month+'.nc')
            v = UKV.Tair[:,:,38:] 
            v = v.interp({'x':x_cci,'y':y_cci},method=method)
            v.to_netcdf(odir+'Tair'+year+month+'.nc')
            UKV = xr.open_dataset(idir+'Wind'+year+month+'.nc')
            v = UKV.Wind[:,:,38:] 
            v = v.interp({'x':x_cci,'y':y_cci},method=method)
            v.to_netcdf(odir+'Wind'+year+month+'.nc')
        except:
            pass

