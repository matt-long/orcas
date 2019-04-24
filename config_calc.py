#-- os interaction
import os
from subprocess import check_call
import socket
import re

def _get_host():
    """Function to determine which base class to use.
    """
    filter1 = r'^cheyenne'
    filter2 = r'r\d(i)\d(n)\d*'
    cheyenne_filter = re.compile('|'.join([filter1, filter2]))
    dav_filter = re.compile(r'^casper')
    mac_filter = re.compile(r'^alpenhorn')

    hostname = socket.gethostname()

    host_on_cheyenne = cheyenne_filter.search(hostname)
    host_on_dav = dav_filter.search(hostname)
    host_on_mac = mac_filter.search(hostname)

    if host_on_cheyenne:
        return 'cheyenne'

    elif host_on_dav:
        return 'dav'

    elif host_on_mac:
        return 'mac'
    else:
        raise ValueError('cannot determine host')

USER = os.environ['USER']
calc_name = 'orcas'

host = _get_host()
if host in ['cheyenne', 'dav']:
    scratch = f'/glade/scratch/{USER}'
    dataroot = f'/glade/work/{USER}'
    dataout = '/glade/p/eol/stephens/longcoll/mclong_calcs'
    dataroot2 = f'/glade/p/eol/stephens/longcoll'

elif host == 'mac':
    scratch = '/Users/mclong/scratch'
    dataroot = '/Users/mclong/data'
    dataout = f'/Users/mclong/data/calcs/{calc_name}'
    dataroot2 = '/Users/mclong/data/calcs/orcas'

#-- directories
diro = {}
diro['work'] = f'{scratch}/calcs/{calc_name}'
diro['out'] = dataout
diro['fig'] = 'fig'
diro['log'] = 'logs'

for pth in diro.values():
    if not os.path.exists(pth):
        check_call(['mkdir','-p',pth])
