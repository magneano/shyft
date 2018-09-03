# This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
# See file COPYING for more details **/
from shyft import api
try:
    from cftime import utime
except ModuleNotFoundError as mnf:
    from netcdftime import utime

import numpy as np

""" These are the current supported regular time-step intervals """
delta_t_dic = {'days': 24*3600, 'hours': 3600, 'minutes': 60,
               'seconds': 1}


def convert_netcdf_time(time_spec, t):
    """
    Converts supplied numpy array to  shyft time given netcdf time_spec.
    Throws exception if time-unit is not supported, i.e. not part of delta_t_dic
    as specified in this file.

    Parameters
    ----------
        time_spec: string
           from netcdef  like 'hours since 1970-01-01 00:00:00'
        t: numpy array
    Returns
    -------
        numpy array type int64 with new shyft time units (seconds since 1970utc)
    """
    u = utime(time_spec)
    t_origin = api.Calendar(int(u.tzoffset)).time(u.origin.year, u.origin.month, u.origin.day, u.origin.hour, u.origin.minute, u.origin.second).seconds
    delta_t = delta_t_dic[u.units]
    return t_origin + delta_t * t[:].astype(np.int64)
