#! /usr/bin/env python
import xarray as xr
from workflow.argpass import pickleparse

#-------------------------------------------------------------------------------
#-- set defaults
#-------------------------------------------------------------------------------

control = pickleparse(default={'file_in':None,
                               'file_out':None,
                               'function':None,
                               'kwargs':{},
                               'isel': {}},
                      description='Transform dataset',
                      required_parameters=['file_in','file_out',
                                           'function'])
#-------------------------------------------------------------------------------
#-- parse input args
#-------------------------------------------------------------------------------

file_in = control['file_in']
file_out = control['file_out']
func_str = control['function']
isel = control['isel']
kwargs = control['kwargs']

#-- file_in is should be a list
if isinstance(file_in,str):
    file_in = [file_in]

#-------------------------------------------------------------------------------
#-- assign the func_str to an actual function
#-------------------------------------------------------------------------------

if func_str == 'compute_potential_temperature':
    import cam
    func = lambda ds: cam.potential_temperature(cam.pres_hybrid(ds)['Pm'],
                                                ds['T'],return_type='Dataset')
elif func_str == 'compute_layer_pressure':
    import cam
    func = lambda ds: cam.pres_hybrid(ds,layer_center=True,layer_interface=False)

elif func_str == 'compute_interface_pressure':
    import cam
    func = lambda ds: cam.pres_hybrid(ds,layer_center=False,layer_interface=True)

elif func_str == '80W':
    func = lambda ds: ds.sel(lon=-80.+360.,lat=slice(-90,-30))

elif func_str == '170E':
    func = lambda ds: ds.sel(lon=170.,lat=slice(-90,-30))

elif func_str == 'scargo_profiles':
    def func(ds):
        import cam
        import numpy as np
        waypoints = {'Christchurch':[172.+32./60.,-43+29./60.,],
                     '170E,60S' :[170.,-60.],
                     '170E,65S' :[170.,-65.],
                     '170E,70S' :[170.,-70.],
                     'McMurdo' : [166+28./60.,-77+51/60.],
                     'Beardmore' : [164+23./60.,-84.],
                     'South Pole': [139+16/60.,-90.]}

        profile  = [n for n in waypoints.keys()]
        lon = np.array([])
        lat = np.array([])
        for n in profile:
            lonlat = waypoints[n]
            lon = np.concatenate((lon[:],[lonlat[0]]))
            lat = np.concatenate((lat[:],[lonlat[1]]))
        lon = xr.DataArray(lon,dims=('profile'),coords={'profile':profile})
        lat = xr.DataArray(lat,dims=('profile'),coords={'profile':profile})

        dss = []
        for l in range(len(ds.time)):
            var = [ds[v].isel(time=l) for v in ds.variables if 'lat' in ds[v].dims and 'lon' in ds[v].dims]
            dss.append(cam.interp_columns_esmf(lon,lat,ds.lon,ds.lat,*var))
        dss = xr.concat(dss,dim='time')
        dss['time'] = ds.time
        return dss


elif func_str == 'so_ocean_mean':
    import cam
    func = cam.so_ocean_mean

else:
    raise ValueError('unknown function: {0}'.format(func_str))


#-------------------------------------------------------------------------------
#-- open datasets
#-------------------------------------------------------------------------------

xrods = {'decode_coords':False,'decode_times':False,
         'drop_variables':'time_written'}
ds = {}
for f in file_in:
    print('reading {0}\n'.format(f))
    ds = xr.merge((ds,xr.open_dataset(f,**xrods)))
if isel:
    ds = ds.isel(**isel)

print('-'*80)
print('Input dataset')
ds.info()
print()

#-------------------------------------------------------------------------------
#-- apply transformation
#-------------------------------------------------------------------------------

dso = func(ds,**kwargs)
if ds.time.bounds not in dso.variables:
    dso[ds.time.bounds] = ds[ds.time.bounds]

print('-'*80)
print('Output dataset')
dso.info()
print()

#-------------------------------------------------------------------------------
#-- write result
#-------------------------------------------------------------------------------

print('writing to file: {0}'.format(file_out))
dso.to_netcdf(file_out,unlimited_dims='time')
