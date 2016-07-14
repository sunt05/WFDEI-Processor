# preload packages
import numpy as np
import os
import sys
from datetime import datetime, timedelta, date
from scipy import interpolate
import time
import pandas as pd
from netCDF4 import Dataset



dir=os.path.dirname(os.path.realpath(__file__))
# print(dir)

# determing grid index according to coordinates
def lon_lat_grid(lat, lon):
    lon_deci = lon - int(lon)
    lat_deci = lat - int(lat)

    if lon >= 0:
        if 0 <= lon_deci < 0.5:
            lon = int(lon) + 0.25
        else:
            lon = int(lon) + 0.75
    else:
        if -0.5 < lon_deci <= 0:
            lon = -(-int(lon) + 0.25)
        else:
            lon = -(-int(lon) + 0.75)

    if lat >= 0:
        if 0 <= lat_deci < 0.5:
            lat = int(lat) + 0.25
        else:
            lat = int(lat) + 0.75
    else:
        if -0.5 < lat_deci <= 0:
            lat = -(-int(lat) + 0.25)
        else:
            lat = -(-int(lat) + 0.75)

    return lat, lon


# Get City Index: WATCH
def WATCH_get_city_index(lat, lon):
    nc = Dataset(os.path.join(dir,"WFD-land-lat-long-z.nc"))
    for i in range(0, 67420):
        if nc.variables['Latitude'][i] == lat and nc.variables['Longitude'][i] == lon:
            index = i
            break
    return index


# Get City Index: WFDEI
def WFDEI_get_city_index(lat, lon):
    with open(os.path.join(dir,'WFDEI-land-long-lat-height.txt')) as f:
        ls = [line.split() for line in f]
    for i in range(7, len(ls)):
        if float(ls[i][0]) == lon and float(ls[i][1]) == lat:
            return int(ls[i][4]), int(ls[i][3])
            break


# generate WFDEI filename
def path_WFDEI(directory, var, year, month):
    if var == "Rainf":
        path = os.path.join(directory, var + '_WFDEI_CRU')
        fn = var + '_WFDEI_CRU_' + str(year) + "%02.d" % month + ".nc"
    else:
        path = os.path.join(directory, var + '_WFDEI')
        fn = var + '_WFDEI_' + str(year) + "%02.d" % month + ".nc"

    path = os.path.join(path, fn)
    return path


# import WFDEI data
def read_WFDEI(directory, var, year_start, year_end, xlat, xlon):
    print('reading in ' + var + ':')
    rawdata = []
    for year in range(year_start, year_end + 1):
        print('     working on ' + str(year) + '...')
        for month in range(1, 13):
            # determine file name to read in
            fn = path_WFDEI(directory, var, year, month)
            # print(fn)
            # get WFDEI dataset:
            nc = Dataset(fn)

            # read in data:
            for i in range(0, len(nc.variables[var][:, xlat, xlon])):
                # determine date time string
                date = str(year) + "%02.d" % month + \
                    "%02.d" % (i / 8 + 1) + "%02.d" % (i % 8 * 3)

                # note the staring index in WFDEI data is 1 whereas in python is 0
                # so xlat-1 and xlon-1 are needed.
                rawdata.append(
                    (date, nc.variables[var][:, xlat - 1, xlon - 1][i]))
    # convert to time series
    ts_data = pd.DataFrame(rawdata, columns=['time', var])
    ts_data["time"] = pd.to_datetime(ts_data['time'], format="%Y%m%d%H")
    print(var + ' successfully imported!')
    return ts_data

# load and rearrange WFDEI data for interpolation


def load_WFDEI_3h(WFDEI_path, year_start, year_end, lat, lon):
    print('*************** WFDEI Data Loader *************** ')
    print('start year: ' + str(year_start))
    print('end year: ' + str(year_end))
    print('input coordinates: (' + '%.2f' % lat + ', ' + '%.2f' % lon + ')')
    # convert user input coordinates into grid index in WFDEI dataset
    glat, glon = lon_lat_grid(lat, lon)
    xlat, xlon = WFDEI_get_city_index(glat, glon)

    # WFDEI variables
    var_list = ["SWdown", "LWdown", "Rainf", "Tair", "PSurf", "Wind", "Qair"]
    # import WFDEI raw data
    data_raw = {var: read_WFDEI(
        WFDEI_path, var, year_start, year_end, xlat, xlon) for var in var_list}
    # join all variable time series into one DataFrame
    try:
        ts_data_raw = [k for k in data_raw.itervalues()]
    except AttributeError as e:
        ts_data_raw = [k for k in data_raw.values()]

    ts_data_raw = [xx.set_index('time') for xx in ts_data_raw]
    data_raw_3h = pd.concat(ts_data_raw, axis=1, join='outer')

    print('*************** WFDEI Data Successfully Loaded *************** ')
    return data_raw_3h
