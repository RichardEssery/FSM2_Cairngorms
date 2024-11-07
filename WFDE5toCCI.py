import numpy as np
from pyproj import Transformer
import xarray as xr
    
# input and output directories for meteorological data 
idir = '/exports/csce/datastore/geos/groups/boreal/WFDE5/'
odir = 'WFDE5_Cairngorms/'

# ancillary file defining snow cci grid for Cairngorms domain
CCI = xr.open_dataset(odir+'ancil.nc')
lons = CCI.lon[:]
lats = CCI.lat[:]

WFDE5 = xr.open_dataset(idir+'2017/LWdown_WFDE5_CRU_201701_v2.0.nc')
v = WFDE5.LWdown[:,292:296,351:355] 
v = v.interp({'lon':lons,'lat':lats},method='linear')
v.to_netcdf('LWdown.nc')
print(v)

# interpolate WFDE5 data onto snow cci grid
method = 'nearest'
for y in range(6):
    year = str(2016+y)
    for m in range(12):
        month = str(m+1)
        if m<9: month = '0'+month
        try:
            print(year,month)
            WFDE5 = xr.open_dataset(idir+year+'/LWdown_WFDE5_CRU_'+year+month+'_v2.0.nc')
            v = WFDE5.LWdown[:,292:296,351:355] 
            v = v.interp({'lon':lons,'lat':lats},method=method)
            v.to_netcdf(odir+'LWdown'+year+month+'.nc')
            WFDE5 = xr.open_dataset(idir+year+'/PSurf_WFDE5_CRU_'+year+month+'_v2.0.nc')
            v = WFDE5.PSurf[:,292:296,351:355] 
            v = v.interp({'lon':lons,'lat':lats},method=method)
            v.rename('Psurf')
            v.to_netcdf(odir+'Psurf'+year+month+'.nc')
            WFDE5 = xr.open_dataset(idir+year+'/Qair_WFDE5_CRU_'+year+month+'_v2.0.nc')
            v = WFDE5.Qair[:,292:296,351:355] 
            v = v.interp({'lon':lons,'lat':lats},method=method)
            v.to_netcdf(odir+'Qair'+year+month+'.nc')
            WFDE5 = xr.open_dataset(idir+year+'/Rainf_WFDE5_CRU+GPCC_'+year+month+'_v2.0.nc')
            v = WFDE5.Rainf[:,292:296,351:355] 
            v = v.interp({'lon':lons,'lat':lats},method=method)
            v.to_netcdf(odir+'Rainf'+year+month+'.nc')
            WFDE5 = xr.open_dataset(idir+year+'/Snowf_WFDE5_CRU+GPCC_'+year+month+'_v2.0.nc')
            v = WFDE5.Snowf[:,292:296,351:355] 
            v = v.interp({'lon':lons,'lat':lats},method=method)
            v.to_netcdf(odir+'Snowf'+year+month+'.nc')
            WFDE5 = xr.open_dataset(idir+year+'/SWdown_WFDE5_CRU_'+year+month+'_v2.0.nc')
            v = WFDE5.SWdown[:,292:296,351:355] 
            v = v.interp({'lon':lons,'lat':lats},method=method)
            v.to_netcdf(odir+'SWdown'+year+month+'.nc')
            WFDE5 = xr.open_dataset(idir+year+'/Tair_WFDE5_CRU_'+year+month+'_v2.0.nc')
            v = WFDE5.Tair[:,292:296,351:355] 
            v = v.interp({'lon':lons,'lat':lats},method=method)
            v.to_netcdf(odir+'Tair'+year+month+'.nc')
            WFDE5 = xr.open_dataset(idir+year+'/Wind_WFDE5_CRU_'+year+month+'_v2.0.nc')
            v = WFDE5.Wind[:,292:296,351:355] 
            v = v.interp({'lon':lons,'lat':lats},method=method)
            v.to_netcdf(odir+'Wind'+year+month+'.nc')
        except:
            pass

