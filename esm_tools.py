#! /usr/bin/env python
import os
import sys
import xarray as xr
import numpy as np
from netcdftime import utime
import time
import tempfile
import json

#----------------------------------------------------------------
#-- function
#----------------------------------------------------------------

def require_variables(ds,req_var):
    missing_var_error = False
    for v in req_var:
        if v not in ds:
            print('ERROR: Missing required variable: %s'%v)
            missing_var_error = True
    if missing_var_error:
        sys.exit(1)
