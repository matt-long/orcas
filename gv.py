#! /usr/bin/env python
import os
from subprocess import call
import numpy as np
import xarray as xr
from config_calc import dataroot

fltdata = os.path.join(dataroot,'orcas/gv_merged_data')
datestr = '20170526'
freq = 10

flight = ['ORCASrf01', 'ORCASrf02', 'ORCASrf03', 'ORCASrf04','ORCASrf05',
          'ORCASrf06', 'ORCASrf07', 'ORCASrf08', 'ORCASrf09', 'ORCASrf10',
          'ORCASrf11','ORCASrf13','ORCASrf14', 'ORCASrf15','ORCASrf16',
          'ORCASrf17', 'ORCASrf18','ORCASrf19']

flight_file = [os.path.join(fltdata,'.'.join([name,'merge%d'%freq,datestr,'nc']))
               for name in flight]

flight_groups = {'ORCASrf01-rf19' : [flight.index('ORCASrf01'),flight.index('ORCASrf19')],
                 'ORCASrf01-rf06' : [flight.index('ORCASrf01'),flight.index('ORCASrf06')],
                 'ORCASrf07-rf11' : [flight.index('ORCASrf07'),flight.index('ORCASrf11')],
                 'ORCASrf13-rf19' : [flight.index('ORCASrf13'),flight.index('ORCASrf19')]}

non_research_flights = ['tf01', 'tf02', 'ff01', 'ff02','ff03' ]



#-------------------------------------------------------------------------------
#-- function
#-------------------------------------------------------------------------------

def region_quality_mask(x,y,z,
                        lat_rgn = [-90.,-44],
                        lon_rgn = [-180.,180.],
                        named_points = {
                            'SCCI' : [-53.01062,-70.85168,42],
                            'SCAR' : [-18.3483,-70.3386, 167],
                            'SCTE' : [-41.438611, -73.0939, 294],
                            'SCVD' : [-39.649722,-73.086111,59]}):

    from earth_geometry import points_in_range

    airport_lon = np.array([v[1] for v in named_points.values()])
    airport_lat = np.array([v[0] for v in named_points.values()])

    #-- land_mask = within 10km of airport and below 4 km
    land_mask = ~( (points_in_range(airport_lon,airport_lat,x,y,10.)) & (z < 4.) )

    #-- region_mask
    region_mask = ( (lat_rgn[0] <= y) & (y <= lat_rgn[1]) & (lon_rgn[0] <= x) & (x <= lon_rgn[1]) )

    return ( land_mask & region_mask )

#-------------------------------------------------------------------------------
#-- function
#-------------------------------------------------------------------------------

def open_flightdata(case,mask=True):

    diri = os.path.join(dataroot,'orcas','cesm_flight_data')
    model_files = [os.path.join(diri,'.'.join([case,os.path.basename(f)]))
                  for f in flight_file]

    obs = xr.open_mfdataset(flight_file)
    mdl = xr.open_mfdataset(model_files)

    obs['GGALT'] = obs.GGALT * 1e-3
    obs.GGALT.attrs['units'] = 'km'
    obs = obs.drop(['UTC','DOY'])

    mdl['GGALT'] = obs.GGALT.copy()
    mdl['GGLAT'] = obs.GGLAT.copy()
    mdl['GGLON'] = obs.GGLON.copy()
    mdl = mdl.drop(['UTC','DOY'])

    if mask:
        mdl = mdl.where(region_quality_mask(mdl.GGLON.values,mdl.GGLAT.values,mdl.GGALT.values))
        obs = obs.where(region_quality_mask(obs.GGLON.values,obs.GGLAT.values,obs.GGALT.values))

    return obs,mdl

if __name__ == '__main__':

    #---------------------------------------------
    #--- process flight data
    #---------------------------------------------

    #--- pick out individual flights
    fltdata_merge = os.path.join(fltdata,'ORCASall.merge%d.'%freq+datestr+'.nc')

    ds = xr.open_dataset(fltdata_merge,decode_times=False,decode_coords=False)
    uniflt,fltI = np.unique(ds.flt.values,return_index=True)

    for flt in uniflt:
        if flt <= 0 or flt > 19: continue
        name = 'ORCASrf%02d'%flt
        file_out = flight_file[flight.index(name)]
        print(file_out)
        fnx = np.nonzero(flt == ds.flt.values)[0]
        dsi = ds.isel(time=slice(fnx[0],fnx[-1]+1))
        dsi.to_netcdf(file_out,unlimited_dims='time')
