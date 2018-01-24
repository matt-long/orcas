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

def pres_hybrid(ds,layer_center=True,layer_interface=False):
    '''Calculate pressure at the hybrid levels.

    Parameters
    ----------

    ds : xarray Dataset
       Dataset must contain P0, PS, and hybrid coefficients hya[m,i] hyb[m,i]
    cell_center : logical, optional, default = True
       compute pressure on cell centers
    cell_interface : logical, optional, default = False
       compute pressure on cell interfaces
    '''
    if not layer_center and not layer_interface:
        return ds

    if layer_center:
        hya = 'hyam'
        hyb = 'hybm'
        esm.require_variables(ds,['P0','PS',hya,hyb])
        ds['Pm'] = (ds['P0'] * ds[hya] + ds['PS'] * ds[hyb]) * 0.01 # kg m/m^2/s^2 = Pa, *0.01 convert to hPa
        if 'time' in ds.Pm.dims:
            ds['Pm'] = ds.Pm.transpose('time','lev','lat','lon')

        ds['Pm'].attrs['long_name'] = 'Pressure (layer center)'
        ds['Pm'].attrs['units'] = 'hPa'
    if layer_interface:
        hya = 'hyai'
        hyb = 'hybi'
        esm.require_variables(ds,['P0','PS',hya,hyb])
        ds['Pi'] = (ds['P0'] * ds[hya] + ds['PS'] * ds[hyb]) * 0.01 # kg m/m^2/s^2 = Pa
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

def interp_to_pressure(pressure,*variables,**kwargs):
    '''Interpolate to common pressure levels.

    Parameters
    ----------

    pressure : xarray DataArray
        4d pressure field
    *variables : DataArrays
    **kwargs : options
       method = {'linear','log'}
       pressure_levels = np.array
    '''

    import metpy.calc as mcalc

    #-- parse arguments
    pressure_levels = kwargs.pop('pressure_levels',
                                 np.array([50.,100., 150., 200., 250.,
                                           290., 330., 370., 400., 420.,
                                           480.,500., 520., 540., 560.,
                                           580., 600., 620., 630., 640.,
                                           650.,660., 670., 680., 690.,
                                           700., 710., 720., 730., 740.,
                                           750.,760., 780., 790., 800.,
                                           810., 820., 830., 840., 850.,
                                           860.,870., 880., 890., 900.,
                                           910., 920., 930., 940., 950.,
                                           960.,970., 980., 990.,1000.,
                                           1005.,1010.,1015.,1020.]))
    method = kwargs.pop('method','log')

    dims_in = variables[0].dims
    if dims_in == (u'time', u'lev', u'lat', u'lon'):
        dims_out = (u'time', u'plev', u'lat', u'lon')
        interp_axis = 1
    elif dims_in == (u'lev', u'lat', u'lon'):
        dims_out = (u'plev', u'lat', u'lon')
        interp_axis = 0
    elif dims_in == (u'time',u'lev'):
        dims_out = (u'time',u'plev')
        interp_axis = 1
    else:
        raise ValueError('Bad dimensions for intepolation: {0}'.format(str(dims_in)))

    if not all(DataArray.dims == dims_in for DataArray in variables):
        raise ValueError('All arrays must have the same dims')


    values = [DataArray.values for DataArray in variables]

    if method == 'linear':
        data_on_plev = mcalc.interp(pressure_levels,
                                    pressure.values,
                                    *values, axis=interp_axis)
    elif method == 'log':
        data_on_plev = mcalc.log_interp(pressure_levels,
                                        pressure.values,
                                        *values, axis=interp_axis)
    else:
        raise ValueError('Unknown option for intepolation: {0}'.format(str(method)))

    if isinstance(data_on_plev,list):
        data_on_plev = {da.name:xr.DataArray(arr,dims=dims_out)
                        for da,arr in zip(variables,data_on_plev)}
    else:
        data_on_plev = {variables[0].name:xr.DataArray(data_on_plev,dims=dims_out)}

    coords = {c:variables[0][c] for c in variables[0].coords if c != 'lev'}
    coords['plev'] = xr.DataArray(pressure_levels,dims=('plev'))

    dso = xr.Dataset(data_vars=data_on_plev,coords=coords)

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
        if var.dims == (u'lev', u'lat', u'lon'):
            nk = var.shape[0]
            data_vars[var.name] = xr.DataArray(np.ones((npoints,nk))*np.nan,
                                               dims=dims+('lev',))
            for k in range(nk):
                data_vars[var.name].values[:,k] = griddata(points=(x.ravel(),y.ravel()),
                                                           values=var.values[k,:,:].ravel(),
                                                           xi=(points_lon.values,points_lat.values),
                                                           method='linear')
        elif var.dims == (u'lat', u'lon'):
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
    grid = esmf_tools.create_grid(grid_lon.values, grid_lat.values)
    locstream = esmf_tools.create_locstream_spherical(points_lon.values,points_lat.values)

    coords = {c:points_lon[c] for c in points_lon.coords}
    dims = points_lon.dims
    npoints = len(points_lon)

    data_vars = {}
    srcfield = ESMF.Field(grid,name='srcfield')
    dstfield = ESMF.Field(locstream,name='dstfield')


    regrid = ESMF.Regrid(srcfield, dstfield,
                         regrid_method=ESMF.RegridMethod.BILINEAR,
                         unmapped_action=ESMF.UnmappedAction.ERROR)


    attrs = {}
    for var in variables:
        dstfield.data[...] = np.nan

        if var.dims == (u'lev', u'lat', u'lon'):
            nk = var.shape[0]
            data_vars[var.name] = xr.DataArray(np.ones((npoints,nk))*np.nan,
                                               dims=dims+('lev',),
                                               attrs = var.attrs)
            for k in range(nk):
                srcfield.data[...] = var.values[k,:,:].T
                dstfield = regrid(srcfield, dstfield, zero_region=ESMF.Region.SELECT)
                data_vars[var.name].values[:,k] = dstfield.data

        elif var.dims == (u'lat', u'lon'):
            srcfield.data[...] = var.values[:,:].T
            data_vars[var.name] = xr.DataArray(np.ones((npoints))*np.nan,
                                               dims=dims,
                                               attrs = var.attrs)

            dstfield = regrid(srcfield, dstfield, zero_region=ESMF.Region.SELECT)
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

        #-- if the array is descending (i.e. Altitude, not pressure)
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
