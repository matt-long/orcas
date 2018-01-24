#! /usr/bin/env python
from config_calc import *

dout = os.path.join(diro['out'],'flux_products')

if not os.path.exists(dout):
    call(['mkdir','-p',dout])

files = {'GK01G' : os.path.join(dout,'GarciaKeeling2001_seasonal_Gruber2001_ann-clim-O2_flux-1.125x1.125.nc'),
         'GK01R' : os.path.join(dout,'GarciaKeeling2001_seasonal_ResplandyEtAl2016_ann-clim-O2_flux-1.125x1.125.nc'),
         'TAK09' : os.path.join(dout,'Takahashi2009-clim-CO2_flux-4x5.nc'),
         'LGB17' : os.path.join(dout,'LandschutzerEtAl2017-clim-CO2_flux-1x1.nc')}

src_grid = {'GK01G' : {'grid_type': '1.125x1.125','dlat': 1.125,'dlon': 1.125,
                       'left_lon_corner' : -180.},
            'GK01R' : {'grid_type': '1.125x1.125','dlat': 1.125,'dlon': 1.125,
                                   'left_lon_corner' : -180.},
            'TAK09' : {'grid_type': '4x5','dlat': 4.,'dlon': 5.,
                       'left_lon_corner' : 0.},
            'LGB17' : {'grid_type': '1x1','dlat': 1.,'dlon': 1.,
                       'left_lon_corner' : -180.}}

def files_res(key,grid='native'):
    file_in = files[key]
    if grid == 'native':
        return file_in
    else:
        return file_in.replace('.nc','_to_'+grid+'.nc')

EorW = lambda x: 'W' if x < 0 else 'E'
lon_str = lambda lon: '%d%s'%(abs(lon),EorW(lon))

if __name__ == '__main__':

    import animate_calc
    from regrid import regrid

    clobber = False

#-------------------------------------------------------------------------------
#-- run notebooks to compile fluxes at native resolution
#-------------------------------------------------------------------------------

    notebook = ['gruber_etal_2001_o2_ann_flux.ipynb',
                'resplandy_etal_2016_o2_ann_flux.ipynb',
                'garcia_keeling_2001_o2_fluxes.ipynb',
                'landschutzer_2017.ipynb',
                'takahashi_2009.ipynb']

    if not os.path.exists(files['GK01G']):
        for i,nb in enumerate(notebook):
            print('\n'.join(['-'*80,'[%d] running %s'%(i,nb),'-'*80]))
            out = animate_calc.exec_nb(nb)
            print

#-------------------------------------------------------------------------------
#-- make grid files for each grid
#-------------------------------------------------------------------------------

    src_grids = []
    for d in src_grid.values():
        left_lon_corner_str = lon_str(d['left_lon_corner'])
        grid = '_'.join(['latlon',d['grid_type'],left_lon_corner_str])
        src_grids.append(grid)
        fname = regrid.grid_file(grid)
        if not os.path.exists(fname) or clobber:
            print('grid_name = %s'%grid)
            d.update({'grid_out_fname' : fname})
            ok = regrid.gen_latlon_grid_file(**d)

#-------------------------------------------------------------------------------
#-- remap to destination grid_file
#-------------------------------------------------------------------------------
    rectilinear_grids = [{'grid_name': 'T62',
                          'latlon_file' : '/glade/p/cesm/cseg/inputdata/atm/datm7/NYF/nyf.ncep.T62.050923.nc'},
                          {'grid_name': 'f09',
                          'latlon_file' : '/glade/p/work/mclong/grids/f09_f09.nc'}]

    print('-'*40)
    print('generating rectilinear grids')
    for d in rectilinear_grids:
        grid = '_'.join(['rectilinear',d['grid_name']])
        fname = regrid.grid_file(grid)
        if not os.path.exists(fname) or clobber:
            print('grid_name = %s'%grid)
            d.update({'grid_out_fname' : fname})
            ok = regrid.gen_rectilinear_grid_file(**d)

#-------------------------------------------------------------------------------
#-- remap to destination grid_file
#-------------------------------------------------------------------------------

    dst_grids = ['latlon_1x1_180W','rectilinear_f09','POP_gx1v6']

    interp_method = 'conserve'
    prefill_opt = 'zeros'
    postfill_opt = 'none'

    #-- compute weight files
    for src in src_grids:
        for dst in dst_grids:
            for interp_method in ['bilinear','conserve']:
                wgtFile = regrid.wgt_file(src,dst,interp_method)
                srcGridFile = regrid.grid_file(src)
                dstGridFile = regrid.grid_file(dst)

                if not os.path.exists(wgtFile) or clobber:
                    ok = regrid.gen_weight_file(wgtFile = wgtFile,
                                                srcGridFile = srcGridFile,
                                                dstGridFile = dstGridFile,
                                                InterpMethod = interp_method)

    #-- loop over flux products and remap
    for product,file_in in files.items():

        #-- determine the variables
        with xr.open_dataset(file_in) as ds:
            variables = [str(v) for v in ds.variables
                         if all(c in ds[v].dims for c in ['lat','lon'])]

            time_coordname = {v:'none' for v in variables}
            time_coordname.update({v:'time' for v in variables if 'time' in ds[v].dims})

        d = src_grid[product]
        left_lon_corner_str = lon_str(d['left_lon_corner'])
        src = '_'.join(['latlon',d['grid_type'],left_lon_corner_str])

        #-- loop over destination grids
        for dst in dst_grids:
            wgtFile = regrid.wgt_file(src,dst,interp_method)
            if not os.path.exists(wgtFile):
                print('missing weight file: %s'%wgtFile)
                exit(1)

            file_out = files_res(product,dst.split('_')[1])

            if os.path.exists(file_out) and not clobber:
                continue

            #-- loop over variables
            print(file_in)
            print('\t-->%s'%file_out)
            outfile_opt = 'create'
            for varname_in in variables:
                varname_out = varname_in
                print('\t\t%s'%varname_in)
                ok = regrid.regrid_var(wgtFile = wgtFile,
                                       fname_in = file_in,
                                       varname_in = varname_in,
                                       time_coordname = time_coordname[varname_in],
                                       depth_coordname = 'none',
                                       vert_grid_file = 'none',
                                       fname_out = file_out,
                                       varname_out = varname_out,
                                       src_grid = src,
                                       dst_grid = dst,
                                       postfill_opt = postfill_opt,
                                       prefill_opt = prefill_opt,
                                       outfile_opt = outfile_opt)
                outfile_opt = 'append'
