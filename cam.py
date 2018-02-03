#! /usr/bin/env python
import os
import numpy as np
import xarray as xr
import esm_tools as esm

na = np.newaxis

g   = 9.80616    # m/s^2
mwair = 28.966

re   = 6.37122e06;         # m
deg2rad  = np.pi / 180.0;  # degree/radian
con  = re * deg2rad;       # m

#-------------------------------------------------------------------------------
#--- FUNCTION
#-------------------------------------------------------------------------------

def potential_temperature(pressure,temperature,return_type='DataArray'):

    if 'units' in pressure.attrs:
        if pressure.attrs['units'] in ['hPa','mbar']:
            pot_temp_ref_press = 1000.
        elif pressure.attrs['units'] in ['Pa']:
            pot_temp_ref_press = 10000.
    else:
        pressure_order = np.floor(np.log10(np.nanmean(pressure.max())))
        if pressure_order == 4.:
            pot_temp_ref_press = 10000.
        elif pressure_order == 3.:
            pot_temp_ref_press = 1000.
        else:
            raise ValueError('Cannot determine units of pressure:\n {0}'.format(str(pressure)))

    ptemp = temperature * (pot_temp_ref_press/pressure)**0.286
    ptemp.attrs['long_name'] = 'Potential temperature'
    if 'units' in temperature.attrs:
        ptemp.attrs['units'] = temperature.attrs['units']

    if return_type == 'Dataset':
        ptemp = xr.Dataset({'theta':ptemp})

    return ptemp

#-------------------------------------------------------------------------------
#--- FUNCTION
#-------------------------------------------------------------------------------

def pres_hybrid(ds,layer_center=True,layer_interface=False):
    '''Calculate pressure at the hybrid levels.

    Parameters
    ----------

    ds : xarray Dataset
       Dataset must contain P0,PS,and hybrid coefficients hya[m,i] hyb[m,i]
    cell_center : logical,optional,default = True
       compute pressure on cell centers
    cell_interface : logical,optional,default = False
       compute pressure on cell interfaces
    '''
    if not layer_center and not layer_interface:
        return ds

    if layer_center:
        hya = 'hyam'
        hyb = 'hybm'
        esm.require_variables(ds,['P0','PS',hya,hyb])
        ds['Pm'] = (ds.P0 * ds[hya] + ds.PS * ds[hyb]) * 0.01 # kg m/m^2/s^2 = Pa,*0.01 convert to hPa
        if 'time' in ds.Pm.dims:
            ds['Pm'] = ds.Pm.transpose('time','lev','lat','lon')

        ds['Pm'].attrs['long_name'] = 'Pressure (layer center)'
        ds['Pm'].attrs['units'] = 'hPa'
    if layer_interface:
        hya = 'hyai'
        hyb = 'hybi'
        esm.require_variables(ds,['P0','PS',hya,hyb])
        ds['Pi'] = (ds.P0 * ds[hya] + ds.PS * ds[hyb]) * 0.01 # kg m/m^2/s^2 = Pa
        if 'time' in ds.Pm.dims:
            ds['Pi'] = ds.Pi.transpose('time','ilev','lat','lon')

        ds['Pi'].attrs['long_name'] = 'Pressure (interface)'
        ds['Pi'].attrs['units'] = 'hPa'

    attrs = ds.PS.attrs
    ds['PS'] = ds.PS * 0.01
    attrs['units'] = 'hPa'
    ds.PS.attrs = attrs

    return ds

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
                   204.7640686,77.2164535]),dims=('lev'))

    pressure_levels = xr.DataArray(
        np.array([50.,100.,150.,200.,250.,290.,330.,370.,
                  400.,420.,480.,500.,520.,540.,560.,580.,
                  600.,620.,630.,640.,650.,660.,670.,680.,
                  690.,700.,710.,720.,730.,740.,750.,760.,
                  780.,790.,800.,810.,820.,830.,840.,850.,
                  860.,870.,880.,890.,900.,910.,920.,930.,
                  940.,950.,960.,970.,980.,990.,1000.,1005.,
                  1010.,1015.,1020.]),dims=('lev'))

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
    if method not in ['linear','log']:
        raise ValueError('Unknown option for intepolation: {0}'.format(str(method)))
    if method == 'linear':
        interp_func = mcalc.interp
    elif method == 'log':
        interp_func = mcalc.log_interp

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
    data_on_plev = {}
    coords = {}
    for da in variables:
        print('interpolating '+da.name)
        if 'lev' not in da.dims: continue

        data_on_plev[da.name] = xr.apply_ufunc(interp_func,new_levels,
                                               coord_field,da,
                                               kwargs={'axis':interp_axis},
                                               dask='parallelized',
                                               output_dtypes=[np.float32])
        data_on_plev[da.name] = data_on_plev[da.name].rename({'lev':levdim})
        #data_on_plev[da.name] = data_on_plev[da.name].transpose(*dims_out)

    #-- put back variables that don't have a level dimension to interpolate
    data_on_plev.update({da.name:da for da in variables if 'lev' not in da.dims})
    dso = xr.Dataset(data_vars=data_on_plev)#,coords=coords)
    #print data_on_plev[da.name]
    dso[levdim].values = new_levels.values
    #data_on_plev[da.name][levdim] = new_levels
    return dso

#-------------------------------------------------------------------------------
#--- FUNCTION
#-------------------------------------------------------------------------------

def interp_columns(points_lon,points_lat,grid_lon,grid_lat,*variables):

    from scipy.interpolate import griddata

    grid_lon.values = np.where(grid_lon.values>180.,grid_lon.values-360.,grid_lon.values)

    x,y = np.meshgrid(grid_lon,grid_lat)

    points_lon.values = np.where(points_lon.values>180.,points_lon.values-360.,points_lon.values)
    coords = {c:points_lon[c] for c in points_lon.coords}
    dims = points_lon.dims
    npoints = len(points_lon)

    data_vars = {}
    for var in variables:
        print(var.name)
        if var.dims == (u'lev',u'lat',u'lon'):
            nk = var.shape[0]
            data_vars[var.name] = xr.DataArray(np.ones((npoints,nk))*np.nan,
                                               dims=dims+('lev',))
            for k in range(nk):
                data_vars[var.name].values[:,k] = griddata(points=(x.ravel(),y.ravel()),
                                                           values=var.values[k,:,:].ravel(),
                                                           xi=(points_lon.values,points_lat.values),
                                                           method='linear')
        elif var.dims == (u'lat',u'lon'):
            data_vars[var.name] = xr.DataArray(np.ones((npoints))*np.nan,
                                               dims=dims)
            data_vars[var.name].values = griddata(points=(x.ravel(),y.ravel()),
                                                  values=var.values[:,:].ravel(),
                                                  xi=(points_lon.values,points_lat.values),
                                                  method='linear')
        else:
            raise ValueError('Bad dimensions for intepolation: {0}'.format(str(var.dims)))


    dso = xr.Dataset(data_vars=data_vars,coords=coords)
    return dso

#-------------------------------------------------------------------------------
#--- FUNCTION
#-------------------------------------------------------------------------------

def interp_columns_esmf(points_lon,points_lat,grid_lon,grid_lat,*variables):

    import ESMF
    import esmf_tools
    grid = esmf_tools.create_grid(grid_lon.values,grid_lat.values)
    locstream = esmf_tools.create_locstream_spherical(points_lon.values,points_lat.values)

    coords = {c:points_lon[c] for c in points_lon.coords}
    dims = points_lon.dims
    npoints = len(points_lon)

    data_vars = {}
    srcfield = ESMF.Field(grid,name='srcfield')
    dstfield = ESMF.Field(locstream,name='dstfield')


    regrid = ESMF.Regrid(srcfield,dstfield,
                         regrid_method=ESMF.RegridMethod.BILINEAR,
                         unmapped_action=ESMF.UnmappedAction.ERROR)


    attrs = {}
    for var in variables:
        dstfield.data[...] = np.nan

        if var.dims == (u'lev',u'lat',u'lon'):
            nk = var.shape[0]
            data_vars[var.name] = xr.DataArray(np.ones((npoints,nk))*np.nan,
                                               dims=dims+('lev',),
                                               attrs = var.attrs)
            for k in range(nk):
                srcfield.data[...] = var.values[k,:,:].T
                dstfield = regrid(srcfield,dstfield,zero_region=ESMF.Region.SELECT)
                data_vars[var.name].values[:,k] = dstfield.data

        elif var.dims == (u'lat',u'lon'):
            srcfield.data[...] = var.values[:,:].T
            data_vars[var.name] = xr.DataArray(np.ones((npoints))*np.nan,
                                               dims=dims,
                                               attrs = var.attrs)

            dstfield = regrid(srcfield,dstfield,zero_region=ESMF.Region.SELECT)
            data_vars[var.name].values[:] = dstfield.data
        else:
            raise ValueError('Bad dimensions for intepolation: {0}'.format(str(var.dims)))

    dso = xr.Dataset(data_vars=data_vars,coords=coords)
    return dso

#-------------------------------------------------------------------------------
#--- FUNCTION
#-------------------------------------------------------------------------------

def interp_within_column(points_z,model_z,*model_var_columns,**kwargs):

    #from scipy.interpolate import interp1d

    ncolumns = len(points_z)

    model_z_surface = kwargs.pop('model_z_surface',None)

    data_vars = {}
    for var in model_var_columns:
        data_vars[var.name] = xr.DataArray(np.ones((ncolumns))*np.nan,
                                           dims=points_z.dims,
                                           attrs = var.attrs)

        #-- if the array is descending (i.e. Altitude,not pressure)
        #   reverse it
        stride = 1
        if model_z[0,0]>model_z[0,-1]:
            stride = -1

        if var.ndim == 2:
            for c in range(ncolumns):
                z = model_z.values[c,::stride]
                y = var.values[c,::stride]
                data_vars[var.name].values[c] = np.interp(x=points_z.values[c],
                                                          xp=z,
                                                          fp=y,
                                                          left=y[0],
                                                          right=y[-1])
        else:
            data_vars[var.name].values = var.values

    coords = {c:points_z[c] for c in points_z.coords}
    dso = xr.Dataset(data_vars=data_vars,coords=coords)
    return dso

#-------------------------------------------------------------------------------
#--- FUNCTION
#-------------------------------------------------------------------------------

def regional_mean(ds,rmask,**kwargs):
    area = grid_tools.compute_grid_area(ds.lon.values,ds.lat.values)

    wgt = rmask * area
    wgt = wgt/wgt.sum()
    wgt = wgt.compute()

    ravg = (ds * wgt).sum(dim=['lat','lon'])

    return ravg
