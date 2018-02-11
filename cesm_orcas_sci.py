#! /usr/bin/env python
import os
from subprocess import call
from config_calc import dataroot2

#------------------------------------------------------------
#--- constants
#------------------------------------------------------------
mwair = 28.966
mw = {'O2'  : 32.0,
      'CO2' : 44.0,
      'N2'  : 28.0}

#------------------------------------------------------------
#--- case specifications
#------------------------------------------------------------
case_definitions = {
    'bgeos5.B20TRC5CN.f09_g16.BPRD_orcas_sci.003': {'mdl_name':'b13geos5',
                                                    'datestr' : '20070101-20160229'},
    'b13.B20TRC5CN.f09_g16.BPRD_orcas_sci.001' : {'mdl_name':'b13',
                                                  'datestr' : '20070101-20160229'},
    'bgeos5.B20TRC5CN.f09_g16.BPRD_orcas_sci.004' : {'mdl_name':'b13geos5',
                                                     'datestr' : '20070101-20160229'},
    'bgeos5.B20TRC5CN.f09_g16.BPRD_orcas_sci.004a' : {'mdl_name':'b13geos5',
                                                     'datestr' : '20070101-20160229'}}
for case in case_definitions.keys():
    case_definitions[case]['droot'] = os.path.join(dataroot2,'hpss-mirror',case)

#------------------------------------------------------------
#--- FUNCTION
#------------------------------------------------------------
def convert_units(species,offset,units):
    convert = lambda x: x

    if species == 'O2':
        norm_by_chi = 1. / 0.2095
    else:
        norm_by_chi = 1.

    if units == 'kg/kg':
        convert = lambda x: (x - offset) * 1.0e6 * mwair / mw[species] * norm_by_chi
    elif units == 'mol/mol':
        convert = lambda x: (x - offset) * 1.0e6 * norm_by_chi

    units = 'ppmv'
    if species == 'O2':
        units = 'per meg'
    if species == 'dye':
        units = 'fraction'
    if species == 'age':
        units = 'years'

    return convert,units

#------------------------------------------------------------
#--- constituents simulated thru srf_emis mechanism (vmr, mol/mol):
#-----------------------------------------------------------
def trace_gas_tracers(case):
    offset = {}

    #------------------------------------------------------------
    #--- case
    #------------------------------------------------------------
    if case in ['bgeos5.B20TRC5CN.f09_g16.BPRD_orcas_sci.003',
                'b13.B20TRC5CN.f09_g16.BPRD_orcas_sci.001']:

        '''
        first rounds of ORCAS science cases, run late October, 2016.
        these cases included CarbonTracker fluxes, but there was a unit-conversion
        error: the fluxes were supplied in mol/m^2/s (or maybe year), but interpreted
        as molecules/cm^2/s....so only a factor of about 10^23 off.  The other fluxes
        were fine.
        '''

        # initialized concentration (the offset) in mmr
        mmr = {'CO2' : 360.0e-6*44./28.966, # molCO2/molAIR * gCO2/molCO2 * molAIR/gAIR = kg/kg
               'O2'  : 360.0e-6*32./28.966}

        o2info = {'species': 'O2',
                  'units' : 'mol/mol'}
        co2info = {'species': 'CO2',
                  'units' : 'mol/mol'}

        srf_emis = {'O2_GKA' : o2info.copy(),     # Garcia & Keeling
                    'aO2_GKA' : o2info.copy(),    # Garcia & Keeling (abiotic)
                    'CO2_L14C' : co2info.copy(),  # Landschutzer
                    'CO2_T09' : co2info.copy()}   # Takahashi

        names = {'O2_GKA' : 'GK2001',
                 'aO2_GKA' : 'GK2001 (abiotic)',
                 'CO2_L14C' : 'Landschutzer (2014)',
                 'CO2_T09' :  'Takahashi (2009)'}

        dyeinfo = {'species':'dye','units':'mol/mol'}
        srf_emis.update({k:dyeinfo.copy() for k in ['IDL_S%03d'%i for i in range(1,6)]})

        dyeinfo = {'species':'age','units':'mol/mol'}
        srf_emis.update({k:dyeinfo.copy() for k in ['IDL_T%03d'%i for i in range(0,6)]})

        idl_names = {'IDL_S001' : 'Land',
                     'IDL_S002' : 'Ocean 90N-30S',
                     'IDL_S003' : 'Ocean 30S-45S',
                     'IDL_S004' : 'Ocean 45S-60S',
                     'IDL_S005' : 'Ocean 60S-90S',
                     'IDL_T000' : 'Age of air'}

        for tracer,info in srf_emis.items():
            srf_emis[tracer]['long_name'] = tracer
            if 'IDL' in tracer:
                offset[tracer] = 0.
                if tracer in idl_names:
                    srf_emis[tracer]['long_name'] = idl_names[tracer]
            else:
                if tracer in names:
                    srf_emis[tracer]['long_name'] = names[tracer]
                offset[tracer] = mmr[info['species']]/mw[info['species']]*mwair

        #--- constituents specified thru tracer_gas or co2_cycle
        o2info = {'species': 'O2',
                  'units' : 'kg/kg'}
        co2info = {'species': 'CO2',
                  'units' : 'kg/kg'}

        cpl_emis = {'CO2' : co2info.copy(),
                    'CO2_FFF' : co2info.copy(),
                    'CO2_LND' : co2info.copy(),
                    'CO2_OCN' : co2info.copy(),
                    'O2_OCN' : o2info.copy()}

        for tracer,info in cpl_emis.items():
            cpl_emis[tracer]['long_name'] = tracer
            if tracer == 'CO2':
                offset[tracer] = 0.
            elif tracer == 'O2_OCN':
                offset[tracer] = 100.0e-6
            else:
                offset[tracer] = mmr[info['species']]
            cpl_emis[tracer]['long_name'] = tracer

        tracers = srf_emis.copy()
        tracers.update(cpl_emis.copy())

    #------------------------------------------------------------
    #--- case
    #------------------------------------------------------------
    elif case in ['bgeos5.B20TRC5CN.f09_g16.BPRD_orcas_sci.004',
                  'bgeos5.B20TRC5CN.f09_g16.BPRD_orcas_sci.004a',
                  'bmerra.B20TRC5CN.f09_g16.BPRD_orcas_sci.001']:
        ''' round two of ORCAS science cases.
        '''

        # initialized concentration (the offset) in mmr
        mmr = {'CO2' : 360.0e-6*44./28.966, # molCO2/molAIR * gCO2/molCO2 * molAIR/gAIR = kg/kg
               'O2'  : 360.0e-6*32./28.966}

        o2info = {'species': 'O2',
                  'units' : 'mol/mol'}
        co2info = {'species': 'CO2',
                   'units' : 'mol/mol'}

        srf_emis = {'O2_GKA' : o2info.copy(),     # Garcia & Keeling
                    'aO2_GKA' : o2info.copy(),    # Garcia & Keeling (abiotic)
                    'CO2_L14C' : co2info.copy(),  # Landschutzer
                    'CO2_T09' : co2info.copy(),   # Takahashi
                    'CO2_C15O' : co2info.copy(),  # Carbontracker 2015 (ocean)
                    'CO2_C15L' : co2info.copy(),  # Carbontracker 2015 (land)
                    'CO2_C15F' : co2info.copy(),  # Carbontracker 2015 (fossil)
                    'CO2_C15T' : co2info.copy(),  # Carbontracker 2015 (total)
                    'CO2_CROO' : co2info.copy(),  # Carbontracker NRT Optimized (ocean)
                    'CO2_CROL' : co2info.copy(),  # Carbontracker NRT Optimized (land)
                    'CO2_CROF' : co2info.copy(),  # Carbontracker NRT Optimized (fossil)
                    'CO2_CROT' : co2info.copy(),  # Carbontracker NRT Optimized (total)
                    'CO2_CRPO' : co2info.copy()}  # Carbontracker NRT Priors (ocean)

        names = {'O2_GKA' : 'GK2001',
                 'aO2_GKA' : 'GK2001 (abiotic)',
                 'CO2_L14C' : 'Landschutzer (2014)',
                 'CO2_T09' :  'Takahashi (2009)',
                 'CO2_C15O' : 'CT2015 (ocean)',
                 'CO2_C15L' : 'CT2015 (land)',
                 'CO2_C15F' : 'CT2015 (FF)',
                 'CO2_C15T' : 'CT2015 (total)',
                 'CO2_CROO' : 'CT-NRT (ocean)',
                 'CO2_CROL' : 'CT-NRT (land)',
                 'CO2_CROF' : 'CT-NRT (FF)',
                 'CO2_CROT' : 'CT-NRT (total)',
                 'CO2_CRPO' : 'CT-NRT (ocean,prior)'}

        q = 0
        for i,scale_by in enumerate([0.5,0.25,-0.25,-0.5,-0.75,-1.5]):
            tmpd = co2info.copy()
            tmpd['long_name'] = 'Takahashi (Dec x %+.0f%%)'%(scale_by*100)
            srf_emis.update({'CO2_T09'+chr(q+97) : tmpd.copy()})
            q += 1

            tmpd['long_name'] = 'Takahashi (Jan x %+.0f%%)'%(scale_by*100)
            srf_emis.update({'CO2_T09'+chr(q+97) : tmpd.copy()})
            q += 1


        dyeinfo = {'species':'dye','units':'mol/mol'}
        srf_emis.update({k:dyeinfo.copy() for k in ['IDL_S%03d'%i for i in range(1,7)]})

        dyeinfo = {'species':'age','units':'mol/mol'}
        srf_emis.update({k:dyeinfo.copy() for k in ['IDL_T%03d'%i for i in range(0,7)]})

        idl_names = {'IDL_S001' : 'Land (except Antarctic)',
                     'IDL_S002' : 'Anarctica',
                     'IDL_S003' : 'Ocean 90N-30S',
                     'IDL_S004' : 'Ocean 30S-45S',
                     'IDL_S005' : 'Ocean 45S-60S',
                     'IDL_S006' : 'Ocean 60S-90S',
                     'IDL_T000' : 'Age of air',
                     'IDL_T001' : 'Age of air (Land, except Antarctic)',
                     'IDL_T002' : 'Age of air (Anarctica)',
                     'IDL_T003' : 'Age of air (Ocean 90N-30S)',
                     'IDL_T004' : 'Age of air (Ocean 30S-45S)',
                     'IDL_T005' : 'Age of air (Ocean 45S-60S)',
                     'IDL_T006' : 'Age of air (Ocean 60S-90S)'}

        for tracer,info in srf_emis.items():
            if not 'long_name' in srf_emis[tracer]:
                srf_emis[tracer]['long_name'] = tracer

            if 'IDL' in tracer:
                offset[tracer] = 0.
                if tracer in idl_names:
                    srf_emis[tracer]['long_name'] = idl_names[tracer]
            else:
                if tracer in names:
                    srf_emis[tracer]['long_name'] = names[tracer]
                offset[tracer] = mmr[info['species']]/mw[info['species']]*mwair

        #--- constituents specified thru tracer_gas or co2_cycle
        o2info = {'species': 'O2',
                  'units' : 'kg/kg'}
        co2info = {'species': 'CO2',
                  'units' : 'kg/kg'}

        cpl_emis = {'CO2' : co2info.copy(),
                    'CO2_FFF' : co2info.copy(),
                    'CO2_LND' : co2info.copy(),
                    'CO2_OCN' : co2info.copy(),
                    'O2_OCN' : o2info.copy()}

        for tracer,info in cpl_emis.items():
            cpl_emis[tracer]['long_name'] = tracer
            if tracer == 'CO2':
                offset[tracer] = 0.
            elif tracer == 'O2_OCN':
                offset[tracer] = 100.0e-6
            else:
                offset[tracer] = mmr[info['species']]


        tracers = srf_emis.copy()
        tracers.update(cpl_emis.copy())

    elif case == 'b.e112.B1850LENS.f09_g16.abio-gas.005':
        ''' cesm 1.1 abiotic gas simulation
        '''

        # initialized concentration (the offset) in mmr
        # molCONSTITUENT/molAIR * gCONSTITUENT/mol * molAIR/gAIR = kg/kg
        mmr = {'N2'  : 100.0e-6*28./28.966}

        o2info = {'species': 'O2',
                  'units' : 'mol/mol'}
        co2info = {'species': 'CO2',
                  'units' : 'mol/mol'}

        #--- constituents specified thru tracer_gas or co2_cycle
        o2info = {'species': 'O2',
                  'units' : 'kg/kg'}
        co2info = {'species': 'CO2',
                  'units' : 'kg/kg'}

        cpl_emis = {'aN2_OCN' : {'species':'N2','units':'kg/kg'}}

        for tracer,info in cpl_emis.items():
            cpl_emis[tracer]['long_name'] = tracer
            offset[tracer] = mmr[info['species']]

            tracers = cpl_emis.copy()

    else:
        raise ValueError('Unknown case: {0}'.format(case))

    for tracer,info in tracers.items():
        tracers[tracer]['convert'],tracers[tracer]['units'] = \
            convert_units(info['species'],offset[tracer],info['units'])
    return tracers

#-------------------------------------------------------------------------------
#-- function
#-------------------------------------------------------------------------------

def open_casedata(case,component,stream,variables,
                  transformed=''):

    import xarray as xr

    case_def = case_definitions[case]
    subdir = 'proc/tseries/daily'

    if transformed:
        subdir = subdir+'_'+transformed

    ds = {}
    for v in variables:
        f = os.path.join(case_def['droot'],component,subdir,'.'.join([case,stream,v,case_def['datestr'],'nc']))
        ds = xr.merge((ds,xr.open_dataset(f,drop_variables=['time_written',
                                                            'date_written',
                                                            'date','datesec',
                                                            'mdt','nbsec','nsbase',
                                                            'nsteph','ntrk','ntrm',
                                                            'ntrn'])))

    return ds

#-------------------------------------------------------------------------------
#-- function
#-------------------------------------------------------------------------------

def convert_dataset(ds,case):
    tracer_info = trace_gas_tracers(case)
    dso = ds.copy()

    for v in ds.variables:

        if v in tracer_info:
            dso[v] = tracer_info[v]['convert'](ds[v])
            dso[v].attrs['units'] = tracer_info[v]['units']
            dso[v].attrs['long_name'] = tracer_info[v]['long_name']

        elif 'SFCO2' in v and ds[v].attrs['units'] == 'kg/m2/s':
            dso[v] = ds[v] * 1000. / 44. * 86400. * 365.
            dso[v].attrs['units'] = 'mol m$^{-2}$ yr$^{-1}$'
            dso[v].attrs['long_name'] = v+' surface flux'

        elif 'SFO2' in v and ds[v].attrs['units'] == 'kg/m2/s':
            dso[v] = ds[v] * 1000. / 32. * 86400. * 365.
            dso[v].attrs['units'] = 'mol m$^{-2}$ yr$^{-1}$'
            dso[v].attrs['long_name'] = v+' surface flux'

    return dso

#-------------------------------------------------------------------------------
#--- PREPROCESS
#-------------------------------------------------------------------------------
if __name__ == '__main__':

    from workflow import chunktime as ct
    from workflow import task_manager as tm
    from workflow.argpass import picklepass

    clobber = False

    case = 'bgeos5.B20TRC5CN.f09_g16.BPRD_orcas_sci.004a'
    component = 'atm'
    stream = 'cam.h0'
    case_def = case_definitions[case]

    chunk_size = 100

    #---------------------------------------------------------------------------
    #-- add derived variables
    #---------------------------------------------------------------------------

    diri = os.path.join(case_def['droot'],component,'proc/tseries/daily')
    diro = diri

    fi = os.path.join(diri,'.'.join([case,stream,'{varname}',case_def['datestr'],'nc']))
    fo = os.path.join(diro,'.'.join([case,stream,'{varname}',case_def['datestr'],'nc']))

    script = os.path.abspath('./transform_dataset.py')

    derived_vars = {'theta' : {'file_in' : [fi.format(varname='T'), fi.format(varname='PS')],
                               'func' : 'compute_potential_temperature'},
                    'Pm' :  {'file_in' : [fi.format(varname='PS')],
                             'func' : 'compute_pressure'}
                    }

    for vnew,specs in derived_vars.items():

        control = {'file_in':specs['file_in'],
                   'file_out':fo.format(varname=vnew),
                   'function':specs['func']}

        #jid = tm.submit([script,picklepass(control)],
        #                memory='50GB')

        jid = ct.apply(script=script,
                       kwargs = control,
                       chunk_size = chunk_size,
                       clobber=clobber,
                       cleanup=True,
                       submit_kwargs_i={'memory':'60GB','constraint':'geyser'},
                       submit_kwargs_cat={'memory':'100GB'})
    tm.wait()

    #---------------------------------------------------------------------------
    #-- remap to new vertical coordinate
    #---------------------------------------------------------------------------

    diro = os.path.join(case_def['droot'],component,'proc/tseries/daily_z3')
    if not os.path.exists(diro):
        call(['mkdir','-p',diro])

    fi = os.path.join(diri,'.'.join([case,stream,'{varname}',case_def['datestr'],'nc']))
    fo = os.path.join(diro,'.'.join([case,stream,'{varname}',case_def['datestr'],'nc']))

    file_in_z3 = fi.format(varname='Z3')

    variables = ['Z3','Q','U','V','Pm','theta']+[k for k in trace_gas_tracers(case)]

    script = os.path.abspath('./calc_remap_vertical_coord.py')

    for v in variables:
        file_in = fi.format(varname=v)
        file_out = fo.format(varname=v)

        control = {'file_in':file_in,
                   'file_out':file_out,
                   'file_in_vertical_coord':file_in_z3,
                   'remap_variables':[v],
                   'coord_field_name':'Z3'}

        jid = ct.apply(script=script,
                       kwargs = control,
                       chunk_size = chunk_size,
                       clobber=clobber,
                       cleanup=True,
                       submit_kwargs_i={'memory':'100GB','constraint':'geyser'},
                       submit_kwargs_cat={'memory':'300GB'})
    tm.wait()

    #---------------------------------------------------------------------------
    #-- transform: regional mean
    #---------------------------------------------------------------------------

    diri = os.path.join(case_def['droot'],component,'proc/tseries/daily')
    diri_z3 = os.path.join(case_def['droot'],component,'proc/tseries/daily_z3')

    diro = os.path.join(case_def['droot'],component,'proc/tseries/daily_so_ocean_mean')
    if not os.path.exists(diro):
        call(['mkdir','-p',diro])

    fiz3 = os.path.join(diri_z3,'.'.join([case,stream,'{varname}',case_def['datestr'],'nc']))
    fi = os.path.join(diri,'.'.join([case,stream,'{varname}',case_def['datestr'],'nc']))
    fo = os.path.join(diro,'.'.join([case,stream,'{varname}',case_def['datestr'],'nc']))

    script = os.path.abspath('./transform_dataset.py')

    variables = ['Z3','Q','U','V','Pm','theta']+[k for k in trace_gas_tracers(case)]
    variables += ['SF'+k for k in trace_gas_tracers(case) if 'IDL' not in k]

    for v in variables:
        if 'SF' in v:
            file_in = fi.format(varname=v)
        else:
            file_in = fiz3.format(varname=v)

        file_out = fo.format(varname=v)

        control = {'file_in':file_in,
                   'file_out':file_out,
                   'function':'so_ocean_mean',
                   'kwargs' : {'varlist':[v]}}
        jid = ct.apply(script=script,
                       kwargs = control,
                       chunk_size = chunk_size,
                       clobber=clobber,
                       cleanup=True,
                       submit_kwargs_i={'memory':'60GB','constraint':'geyser'},
                       submit_kwargs_cat={'memory':'100GB'})
    tm.wait()

    #---------------------------------------------------------------------------
    #-- transform: 80W
    #---------------------------------------------------------------------------

    diri = os.path.join(case_def['droot'],component,'proc/tseries/daily')
    diro = os.path.join(case_def['droot'],component,'proc/tseries/daily_80W')
    if not os.path.exists(diro):
        call(['mkdir','-p',diro])

    fi = os.path.join(diri,'.'.join([case,stream,'{varname}',case_def['datestr'],'nc']))
    fo = os.path.join(diro,'.'.join([case,stream,'{varname}',case_def['datestr'],'nc']))

    script = os.path.abspath('./transform_dataset.py')

    variables = ['Z3','Q','U','V','Pm','theta']+[k for k in trace_gas_tracers(case)]
    for v in variables:
        file_in = fi.format(varname=v)
        file_out = fo.format(varname=v)

        control = {'file_in':file_in,
                   'file_out':file_out,
                   'function':'80W'}

        jid = ct.apply(script=script,
                       kwargs = control,
                       chunk_size = chunk_size,
                       clobber=clobber,
                       cleanup=True,
                       submit_kwargs_i={'memory':'60GB','constraint':'geyser'},
                       submit_kwargs_cat={'memory':'100GB'})
    tm.wait()
