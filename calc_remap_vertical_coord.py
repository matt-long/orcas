#! /usr/bin/env python
import xarray as xr
import numpy as np

#-------------------------------------------------------------------------------
#--- FUNCTION
#-------------------------------------------------------------------------------

def remap_vertical_coord(coord_field,*variables,**kwargs):
    '''Interpolate to new vertical coordinate.

    Parameters
    ----------

    coord_field : xarray DataArray
        4d coord_field field
    *variables : DataArrays
    **kwargs : options
       method = {'linear','log'}
       new_levels = np.array
    '''

    import metpy.calc as mcalc

    geopotential_height_levels = xr.DataArray(
        np.array([ 42482.4296875,40789.765625,39157.46875,
                   37586.35546875,36074.7578125,34620.71875,
                   33221.05078125,31872.421875,30571.47265625,
                   29308.29492188,28075.45507812,26870.96875,
                   25694.5390625,24545.7109375,23423.09570312,
                   22326.65625,21252.5234375,20197.5078125,
                   19162.2578125,18146.8984375,17149.96484375,
                   16160.68554688,15167.20898438,14164.18066406,
                   13147.97851562,12116.08984375,11066.44824219,
                   9995.7109375,8912.61816406,8007.54882812,
                   7296.19238281,6637.14892578,6022.51513672,
                   5446.35351562,4903.87304688,4391.23974609,
                   3905.2668457,3443.27319336,3074.66967773,
                   2790.08691406,2513.83837891,2245.44848633,
                   1984.4967041,1755.61950684,1581.37988281,
                   1434.59033203,1290.03381348,1147.62548828,
                   1007.2913208,868.96289062,732.5602417,
                   598.00524902,465.22851562,334.16360474,
                   204.7640686,77.2164535]),
        dims=('zlev'),
        attrs={'long_name':'Geopotential height (above sea level)',
               'units':'m'})

    pressure_levels = xr.DataArray(
        np.array([50.,100.,150.,200.,250.,290.,330.,370.,
                  400.,420.,480.,500.,520.,540.,560.,580.,
                  600.,620.,630.,640.,650.,660.,670.,680.,
                  690.,700.,710.,720.,730.,740.,750.,760.,
                  780.,790.,800.,810.,820.,830.,840.,850.,
                  860.,870.,880.,890.,900.,910.,920.,930.,
                  940.,950.,960.,970.,980.,990.,1000.,1005.,
                  1010.,1015.,1020.]),
          dims=('plev'),
          attrs={'long_name':'Pressure',
                 'units':'hPa'})


    #-- parse arguments
    if coord_field.name == 'Z3':
        new_levels = kwargs.pop('new_levels',geopotential_height_levels)
        levdim = u'zlev'
    else:
        new_levels = kwargs.pop('new_levels',pressure_levels)
        levdim = u'plev'

    if len(new_levels) != len(coord_field.lev):
        raise ValueError('new_levels is required to have some number of levels as lev')

    method = kwargs.pop('method','log')
    if method == 'linear':
        interp_func = mcalc.interp
    elif method == 'log':
        interp_func = mcalc.log_interp
    else:
        raise ValueError('Unknown option for intepolation: {0}'.format(str(method)))

    #-- determine output dims
    dims_in = variables[0].dims
    if dims_in == (u'time',u'lev',u'lat',u'lon'):
        dims_out = (u'time',levdim,u'lat',u'lon')
        interp_axis = 1

    elif dims_in == (u'lev',u'lat',u'lon'):
        dims_out = (levdim,u'lat',u'lon')
        interp_axis = 0

    elif dims_in == (u'time',u'lev'):
        dims_out = (u'time',levdim)
        interp_axis = 1

    else:
        raise ValueError('Bad dimensions for intepolation: {0}'.format(str(dims_in)))

    #-- check that all vars with level dimension have the same dimensions
    if not all(da.dims == dims_in for da in variables if 'lev' in da.dims):
        raise ValueError('All arrays must have the same dims')

    #-- loop over vars and interpolate
    data_on_new_lev = {}
    coords = {}
    for da in variables:

        if 'lev' not in da.dims: continue
        print('interpolating '+da.name)
        data_on_new_lev[da.name] = xr.DataArray(interp_func(new_levels.values,
                                                            coord_field.values,
                                                            da.values,
                                                            axis=interp_axis),
                                                dims=dims_out,
                                                attrs=da.attrs)

    #-- put back variables that don't have a level dimension to interpolate
    data_on_new_lev.update({da.name:da for da in variables if 'lev' not in da.dims})

    #-- get coordinates
    coords = {c:variables[0][c] for c in variables[0].coords if c != 'lev'}
    coords['plev'] = new_levels

    return xr.Dataset(data_vars= data_on_new_lev,coords=coords)


if __name__ == '__main__':

    from workflow.argpass import pickleparse

    #-- set defaults
    control = pickleparse(default={'file_in':None,
                                   'file_out':None,
                                   'file_in_vertical_coord':None,
                                   'coord_field_name':None,
                                   'remap_variables':None,
                                   'isel': {}},
                          description='Remap vertical coordinate',
                          required_parameters=['file_in','file_out',
                                               'coord_field_name',
                                               'remap_variables'])

    file_in = control['file_in']
    file_out = control['file_out']
    if 'file_in_vertical_coord' in control:
        file_in_vertical_coord = control['file_in_vertical_coord']
    else:
        file_in_vertical_coord = control['file_in']
    coord_field_name = control['coord_field_name']
    remap_variables = control['remap_variables']
    isel = control['isel']

    xrods = {'decode_coords':False,'decode_times':False}
    coord_field = xr.open_dataset(file_in_vertical_coord,**xrods)[coord_field_name]
    ds = xr.open_dataset(file_in,**xrods)

    if isel:
        isel_kwargs = {}
        for k,v in isel.items():
            if isinstance(v,dict):
                isel_kwargs[k] = slice(v['start'],v['stop'])
            else:
                isel_kwargs[k] = v

        coord_field = coord_field.isel(**isel_kwargs)
        ds = ds.isel(**isel_kwargs)

    dso = remap_vertical_coord(coord_field,*[ds[v] for v in remap_variables])
    print(dso)
    dso.to_netcdf(file_out,unlimited_dims='time')
