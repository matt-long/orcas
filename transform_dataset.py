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
elif func_str == 'compute_pressure':
    import cam
    func = cam.pres_hybrid

elif func_str == '80W':
    func = lambda ds: ds.sel(lon=-80.+360.,lat=slice(-90,-30))

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
print

#-------------------------------------------------------------------------------
#-- apply transformation
#-------------------------------------------------------------------------------

dso = func(ds,**kwargs)
print('-'*80)
print('Output dataset')
dso.info()
print

#-------------------------------------------------------------------------------
#-- write result
#-------------------------------------------------------------------------------

print('writing to file: {0}'.format(file_out))
dso.to_netcdf(file_out,unlimited_dims='time')
