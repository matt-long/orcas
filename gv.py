#! /usr/bin/env python
import os
from subprocess import call
import numpy as np
import xarray as xr

fltdata = '/glade/p/work/mclong/orcas/gv_merged_data'
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
