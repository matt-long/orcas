#! /usr/bin/env python
import os
from subprocess import call
from datetime import datetime, timedelta
import tempfile
from glob import glob

from workflow import task_manager as tm

os.putenv('PYTHONUNBUFFERED','no')

case_info = {'bmerra.e12.B20TRC5CN.f09_g16.BPRD.002' :
             {'date_range' : [datetime(1979,1,1),datetime(2003,12,31)],
              'droot' : '/glade/scratch/mclong/hpss-mirror',
              'startup_delay' : True },
             #
             'bmerra.e12.B20TRC5CN.f09_g16.BPRD.003' :
             {'date_range' : [datetime(2004,1,1),datetime(2015,12,31)],
              'droot' : '/glade/scratch/mclong/hpss-mirror',
              'startup_delay' : False},
             #
             'b13.B20TRC5CN.f09_g16.BPRD_orcas_sci.001':
                 {'date_range' : [datetime(2007,1,1),datetime(2016,2,29)],
                  'droot' : '/glade/scratch/mclong/archive',
                  'startup_delay' : False},
             #
             'bgeos5.B20TRC5CN.f09_g16.BPRD_orcas_sci.003':
                 {'date_range' : [datetime(2007,1,1),datetime(2016,2,29)],
                  'droot' : '/glade/scratch/mclong/archive',
                  'startup_delay' : False},
             #
             'bgeos5.B20TRC5CN.f09_g16.BPRD_orcas_sci.004':
                 {'date_range' : [datetime(2007,1,1),datetime(2016,2,29)],
                  'droot' : '/glade/scratch/mclong/hpss-mirror',
                  'startup_delay' : False},
             #
             'bgeos5.B20TRC5CN.f09_g16.BPRD_orcas_sci.004a':
                 {'date_range' : [datetime(2007,1,1),datetime(2016,2,29)],
                  'droot' : '/glade/scratch/mclong/archive',
                  'startup_delay' : False},
             #
             'bmerra.B20TRC5CN.f09_g16.BPRD_orcas_sci.001':
                 {'date_range' : [datetime(2007,1,1),datetime(2016,2,29)],
                  'droot' : '/glade/scratch/mclong/hpss-mirror',
                  'startup_delay' : False}
             }


#---------------------------------------------------------
#--- function
#---------------------------------------------------------
def filedate_list(date1,date2,
                  datelist_type='nday1_monfile',
                  startup_delay = True,
                  calendar = 'noleap'):
    '''
    make a list of date strings appropriate for files at a particular
    output/averaging frequency

    date1: first date in sequence
    date2: last date in sequence
    '''

    if not any([calendar == c for c in ['gregorian','noleap']]):
        print 'calendar option not recognized: '+calendar
        return []

    nday = (date2 - date1).days + 1
    oneday = timedelta(days=1)

    #--- make list of days
    if datelist_type == 'month':
        datelist = ['%04d-%02d'%(y,m)
                    for y in range(date1.year,date2.year+1)
                    for m in range(1,13)
                    if date1 <= datetime(y,m,1) <= date2]

    elif datelist_type == 'nday1_monfile':
        datelist = ['%04d-%02d-01'%(y,m)
                    for y in range(date1.year,date2.year+1)
                    for m in range(1,13)
                    if date1 <= datetime(y,m,1) <= date2]
        if startup_delay:
            datelist[0] = '%04d-%02d-02'%(date1.year,1)

    elif any([datelist_type == d for d in ['nday1','nday1_hour','nhour1_day']]):
        if calendar == 'noleap':
            datelist = [(date1 + timedelta(days=n)).strftime('%Y-%m-%d')
                        for n in range(0, nday)
                        if not ((date1 + timedelta(days=n)).month == 2
                                and (date1 + timedelta(days=n)).day == 29)]

        elif calendar == 'gregorian':
            datelist = [(date1 + timedelta(days=n)).strftime('%Y-%m-%d')
                        for n in range(0, nday)]

        if any([datelist_type == d for d in ['nday1_hour','nhour1_day']]):
            datelist = ['%s-%05d'%(d,0) for d in datelist]


    return datelist

#---------------------------------------------------------
#--- function
#---------------------------------------------------------
def make_timeseries(case,
                    component_stream,
                    output_dir = '',
                    output_case = '',
                    component_dict={'pop':'ocn','cam':'atm'},
                    calendar = 'noleap',
                    clobber = False):
    lJID = []
    err = False
    if not output_case:
        output_case = case[0]

    if not output_dir:
        output_dir = case_info[output_case]['droot']

    for component_name,stream_dict in component_stream.items():

        component_dir = component_dict[component_name]
        print '='*80
        print component_dir
        print '='*80


        for stream,datelist_type in stream_dict.items():

            if any([datelist_type == d for d in
                    ['nday1_monfile','nday1','nday1_hour']]):
                freq = 'daily'
            elif datelist_type == 'month':
                freq = 'monthly'
            else:
                print datelist_type+' unknown'
                exit(1)

            print '.'*80
            print stream+' --> '+freq
            print '.'*80

            opth = '/'.join([output_dir,output_case,component_dir,'proc','tseries',freq])
            if not os.path.isdir(opth): call(['mkdir','-p',opth])

            #-- list the files to concatenate
            dpthlist = []
            hfillist = []
            (fid,tmpfile) = tempfile.mkstemp('','filelist.')

            fid = open(tmpfile,'w')
            for i in range(0,len(case)):

                droot = case_info[case[i]]['droot']
                date_range = case_info[case[i]]['date_range']
                startup_delay = case_info[case[i]]['startup_delay']

                datelist = filedate_list(date_range[0],date_range[1],
                                         datelist_type=datelist_type,
                                         startup_delay=startup_delay,
                                         calendar = calendar)

                ######################################################
                #### BEGIN KLUDGE
                if case[i] == 'bmerra.e12.B20TRC5CN.f09_g16.BPRD.002':
                    if component_name == 'pop' and 'nday1' in stream:
                        datelist.extend(['1997-08-14','1997-12-31'])
                        datelist.sort()
                #### END KLUDGE
                ######################################################

                if i == 0:
                    output_datestr = datelist[0].replace('-','').replace('00000','')

                if i == len(case)-1:
                    output_datestr = output_datestr+'-'+\
                        datelist[-1].replace('-','').replace('00000','').replace('1201','1231')

                dpth = '/'.join([droot,case[i],component_dir,'hist'])
                for d in datelist:
                    hfil = '.'.join([case[i],component_name,stream,d,'nc'])
                    if not os.path.exists(dpth+'/'+hfil):
                        print 'file missing: '+dpth+'/'+hfil
                        return False
                    dpthlist.append(dpth)
                    hfillist.append(hfil)
                    fid.write('%s\n'%(dpth+'/'+hfil))

                hfil_ls = glob(dpth+'/'+'.'.join([case[i],component_name,stream,'????-*','nc']))
                hfil_ls.sort()

                if len(hfil_ls) != len(datelist):
                    for i,hf in enumerate(hfil_ls):
                        d = filedate(hf)
                        if not d in datelist:
                            if not (d < datelist[0] or d > datelist[-1]):
                                print 'unexpected file: '+os.path.basename(hf)
                                exit()
                                err = True

            fid.close()

            tvar,gvar = nc_var_list(dpthlist[0]+'/'+hfillist[0])

            #-- loop over time-dependent variables
            for v in tvar:
                #-- construct ouput filename
                file_out = opth+'/'+'.'.join([output_case,component_name,stream,v,output_datestr,'nc'])

                strvar = ','.join(gvar)+','+v
                if not os.path.exists(file_out) or clobber:
                    cmd = ['cat',tmpfile,'|','ncrcat','-O','-v',strvar,'-o',file_out,';']
                    cmd.extend(['ncks', '-O', '-4', '-L', '1',file_out,file_out,';'])
                    jid = tm.submit(cmd,memory='100GB',modules=['nco'])
                    lJID.append(jid)

    ok = True
    if lJID:
        ok = tm.wait(lJID)
    return ok


#----------------------------------------------------------------
#---- FUNCTION
#----------------------------------------------------------------
def filedate(name):
   import re
   date =  re.findall(r'[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]-00000',name)
   if not date:
       date = re.findall(r'[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]',name)
   if not date:
       date = re.findall(r'[0-9][0-9][0-9][0-9]-[0-9][0-9]',name)
   return date[0]

#------------------------------------
#--- FUNCTION
#------------------------------------
def hpss_put(ldir,hdir):
    cmd = ['hsi','"cd '+hdir+'; '+'cput -h -R '+ldir+'"']
    jid = tm.submit(cmd)


#------------------------------------
#--- FUNCTION
#------------------------------------
def nc_var_list(fname):
    import netCDF4 as nc4
    nc = nc4.Dataset(fname,'r')
    var = nc.variables

    var_list_time = []
    var_list_grid = []
    for v in var:
        if 'time' in nc.variables[v].dimensions \
                and nc.variables[v].dtype == 'float32':
            var_list_time.append(v)
        else:
            var_list_grid.append(v)
    return var_list_time,var_list_grid

#------------------------------------
#--- main
#------------------------------------
if __name__ == '__main__':

    #case = ('bmerra.e12.B20TRC5CN.f09_g16.BPRD.002','bmerra.e12.B20TRC5CN.f09_g16.BPRD.003')
    #stream =  {'cam' : {'h1' : 'nday1_hour',
    #                    'h0' : 'month'}}}

    #CASE = [('b13.B20TRC5CN.f09_g16.BPRD_orcas_sci.001',),
    #        ('bgeos5.B20TRC5CN.f09_g16.BPRD_orcas_sci.003',)]

    CASE = [('bgeos5.B20TRC5CN.f09_g16.BPRD_orcas_sci.004a',)]

    #CASE = [('bmerra.B20TRC5CN.f09_g16.BPRD_orcas_sci.001',)]
    STREAM =  {'cam' : {'h1' : 'nday1_hour',
                        'h0' : 'nday1_hour'},
               'pop' : {'h.ecosys.nday1' : 'nday1_monfile',
                        'h' : 'month',
                        'h.nday1' : 'nday1_monfile'}}

    stream = {k:v for k,v in STREAM.items() if k != 'pop'}

    for case in CASE:
        output_dir = '/glade/scratch/mclong/hpss-mirror'
        output_case = case[0]

        if len(case) > 1:
            if case[0] == 'bmerra.e12.B20TRC5CN.f09_g16.BPRD.002' and \
               case[1] == 'bmerra.e12.B20TRC5CN.f09_g16.BPRD.003':
                output_case = 'bmerra.e12.B20TRC5CN.f09_g16.BPRD.2p3'
            else:
                print 'unknown case concatenation.'
                exit(1)

        print 'working on timeseries:'
        print '\t%s'%case
        print '\t-->%s'%output_case

        ok = make_timeseries(case,stream,
                             output_dir=output_dir,
                             output_case=output_case,
                             calendar = 'gregorian')
        if not ok:
            print 'error'
            exit(1)

        ok = hpss_put(output_dir+'/'+output_case,'/home/mclong/csm/')
