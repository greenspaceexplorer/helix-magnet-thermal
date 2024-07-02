# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 10:32:09 2016

@author: green
"""

from datetime import datetime as dt
from datetime import timedelta as td
import time

def ConvertDate(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

in_name = "crestFlight2011_gpsInfo.csv"
out_name = "HELIX_20191101.txt"
position_in = open(in_name,"r")
position_out = open(out_name,"w")

start_date = dt(year=2019,month=11,day=1)

for line in position_in:
    coord = line.split(",")
    delta_t = td(days=float(coord[0]))
    new_date = start_date + delta_t
    coord[0] = ConvertDate(new_date)
    print("{} E K{} {} {}".format(coord[0],coord[3],coord[1]\
        ,coord[2]),file=position_out)
        
position_in.close()
position_out.close()

print("Converted coordinates in {} to WMM format in {}"\
    .format(in_name,out_name))

    