#-- os interaction
import os
import sys
from glob import glob
from subprocess import call
import socket

#-- analysis
from datetime import datetime
import xarray as xr
import numpy as np

hostname = socket.gethostname()

#-- path additions
path_tools = ['./easy']
for p in path_tools:
    sys.path.insert(0,os.path.abspath(os.path.expanduser(p)))

hostname = socket.gethostname()
if any(s in hostname for s in ['cheyenne','yslogin','geyser','caldera']):
    scratch = '/glade/scratch/'+os.environ['USER']
    dataroot = '/glade/p/work/'+os.environ['USER']
else:
    print('hostname not found')

#-- directories
calc_name = 'orcas'
diro = {}
diro['out'] = os.path.join(scratch,'calcs',calc_name)
diro['fig'] =  'fig'
diro['log'] = './logs'

for pth in diro.values():
    if not os.path.exists(pth):
        call(['mkdir','-p',pth])
