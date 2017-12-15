#! /usr/bin/env python
import os
import sys
import socket
import numpy as np
import xarray as xr

from datetime import datetime

path_tools = ['../']
for p in path_tools:
    sys.path.insert(0,os.path.abspath(os.path.expanduser(p)))

import grid_tools

nowstr = datetime.now().strftime('%Y%m%d')

'''
HTTPSITE=http://bluemoon.ucsd.edu/publications/ralph/airseaflux/Data_files
file=(
    o2flux_ann_global.dat
    o2flux_sea_global.dat
    o2thermal_ann_global.dat
    o2thermal_sea_global.dat
    o2flux_bioann_global.dat
    o2flux_biosea_global.dat)

name=(
    "O2 flux (annual)"
    "O2 flux (monthly anomaly)"
    "O2 flux (annual, thermal component)"
    "O2 flux (monthly anomaly, thermal component)"
    "O2 flux (annual, biological component)"
    "O2 flux (monthly anomaly, biological component)")
'''

droot = '/glade/p/work/mclong/garcia-keeling'
dpm  = np.array([31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.])
eom = np.cumsum(dpm)
bom = np.concatenate((np.array([1]),eom[0:11]+1))

shift_time = 10.
scaleby=0.82


#-------------------------------------------------------------------------------
#-- function
#-------------------------------------------------------------------------------

def dat2nc(dat_file_name):

    nt = 12
    ny = 160
    nx = 320

    dx = 360./nx
    dy = 180./ny

    dss = grid_tools.generate_latlon_grid(nx=320,ny=160,lon0=-180.)

    time = xr.DataArray(np.vstack((bom,eom)).mean(axis=0) + shift_time,
                        dims=('time'),
                        attrs = {'units':'day of year'})

    dss = dss.assign_coords(time = time)

    date = np.round(2000*10000 + np.arange(1,13,1) * 100. + dpm/2.)+shift_time
    dss['date'] = xr.DataArray(date,dims=('time'))


    data = np.loadtxt(dat_file_name).ravel().reshape((nt,ny,nx))

    data = data / dpm[:,None,None] * 365.
    data = data * scaleby
    data[data==0.] = np.nan
    dss['O2_FLUX'] = xr.DataArray(data,dims=('time','lat','lon'),
                                  attrs={'long_name':'O2 flux',
                                         'units':'mol/m^2/yr'})

    print(dss)
    dss.to_netcdf(dat_file_name.replace('.dat','.'+'.'.join(['adjusted',nowstr,'nc'])))
    dss.close()

#-------------------------------------------------------------------------------
#-- main
#-------------------------------------------------------------------------------

if __name__ == '__main__':


    for dat_file_in in ['o2flux_ann_global.dat',
                        'o2flux_sea_global.dat',
                        'o2thermal_ann_global.dat',
                        'o2thermal_sea_global.dat',
                        'o2flux_bioann_global.dat',
                        'o2flux_biosea_global.dat']:

        print('generating %s'%dat_file_in)
        dat2nc(os.path.join(droot,dat_file_in))
        print
