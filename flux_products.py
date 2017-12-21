#! /usr/bin/env python
from config_calc import *

droot = os.path.join(dataroot,'pco2-ldeo')
ds = xr.open_dataset(os.path.join(droot,'ldeo_monthly_clim_v2009_c20150807.nc'))

dout = os.path.join(diro['out'],'flux_products')

if not os.path.exists(dout):
    call(['mkdir','-p',dout])

files = {'GK01G' : os.path.join(dout,'GarciaKeeling2001_seasonal_Gruber2001_ann-clim-O2_flux-1.125x1.125.nc'),
         'GK01R' : os.path.join(dout,'GarciaKeeling2001_seasonal_ResplandyEtAl2016_ann-clim-O2_flux-1.125x1.125.nc'),
         'TAK09' : os.path.join(dout,'Takahashi2009-clim-CO2_flux-5x4.nc'),
         'LGB17' : os.path.join(dout,'LandschutzerEtAl2017-clim-CO2_flux-1x1.nc')}

if __name__ == '__main__':

    import run_notebooks

    notebook = ['gruber_etal_2001_o2_ann_flux.ipynb',
                'resplandy_etal_2016_o2_ann_flux.ipynb',
                'garcia_keeling_2001_o2_fluxes.ipynb',
                'landschutzer_2017.ipynb',
                'takahashi_2009.ipynb']

    for i,nb in enumerate(notebook):
        print('\n'.join(['-'*80,'[%d] running %s'%(i,nb),'-'*80]))
        out = run_notebooks.exec_nb(nb)
        print
