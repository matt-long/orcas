#! /usr/bin/env python
import xarray as xr
from workflow.argpass import pickleparse

#-- set defaults
control = pickleparse(default={'file_in':None,
                               'file_out':None,
                               'function':None,
                               'kwargs':{},
                               'isel': {}},
                      description='Transform dataset',
                      required_parameters=['file_in','file_out',
                                           'function'])

file_in = control['file_in']
file_out = control['file_out']
func_str = control['function']
isel = control['isel']
kwargs = control['kwargs']

if isinstance(file_in,str):
    file_in = [file_in]

if func_str == 'compute_potential_temperature':
    import cam
    func = lambda ds: cam.potential_temperature(cam.pres_hybrid(ds)['Pm'],
                                                ds['T'],return_type='Dataset')
elif func_str == 'compute_pressure':
    import cam
    func = cam.pres_hybrid
    
elif func_str == 'regional_mean':
    landfrac = xr.open_dataset('/glade/p/work/mclong/grids/f09_f09.nc')['LANDFRAC'].isel(time=0)
    rmask = landfrac.where(landfrac<0.9).fillna(0.).where(landfrac>=0.9).fillna(1.).where(landfrac.lat<-44.).fillna(0.)

else:
    raise ValueError('unknown function')


xrods = {'decode_coords':False,'decode_times':False,
         'drop_variables':'time_written'}
ds = {}
for f in file_in:
    print('reading {}\n'.format(f))
    ds = xr.merge((ds,xr.open_dataset(f,**xrods)))
if isel:
    ds = ds.isel(**isel)

print('-'*80)
print('Input dataset')
ds.info()

dso = func(ds,**kwargs)
print('-'*80)
print('Output dataset')
dso.info()

print('writing to file')
dso.to_netcdf(file_out,unlimited_dims='time')
