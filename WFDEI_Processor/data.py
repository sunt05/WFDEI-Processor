import numpy as np
from netCDF4 import Dataset
import os
import datetime
from matplotlib.dates import date2num
import Tair

# Get City Index
def WATCH_get_city_index(city_latitude, city_longtitude):
    nc = Dataset("WFD-land-lat-long-z.nc")
    for i in range(0, 67420):
        if nc.variables['Latitude'][i] == city_latitude and nc.variables['Longitude'][i] == city_longtitude:
            index = i
            break
    return index

# WFDEI Get City Index
def WFDEI_get_city_index(city_latitude, city_longtitude):
    with open('WFDEI-land-long-lat-height.txt') as f:
        ls = [line.split() for line in f]
    for i in range(7, len(ls)):
        if float(ls[i][0]) == city_longtitude and float(ls[i][1]) == city_latitude:
            return int(ls[i][4]), int(ls[i][3])
            break

# Get Tair
def get_data(filepath, index, var):
    Tair = []
    for list in os.listdir(filepath):
        filename = filepath + list
        nc = Dataset(filename)
        days = len(nc.variables[var][:, index])
        for i in range(0, days):
            date = list[-9: -3]
            if i < 9:
                date += '0' + str(i + 1)
            else:
                date += str(i + 1)
            Tair.append((date, nc.variables[var][:, index][i]))

    Tair = (np.array(Tair)).astype(np.float)
    return Tair

# Threshold


def get_threshold(data, percent):
    return np.percentile(data[:, 1], percent)

# get year period


def get_year_period(data):
    return str(data[0][0])[:6], str(data[-1][0])[:6]

# longitude latitude grid


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

# float date to datetime


def date_time(Tair):
    date = []
    for i in Tair:
        date.append(datetime.datetime.strptime(str(int(i[0])), "%Y%m%d"))
    return np.array(date)

# statistics about heatwaves [year, occurrence, cumulative, max, min]


def HW_statistics(HW_index, dates):
    dates_HW = []
    for i in range(0, len(HW_index)):
        dates_HW.append(dates[HW_index[i][0]: HW_index[i][1] + 1])
    dates_HW = np.array(dates_HW)
    # f = open("../London_HW(1961-1990).txt", 'w')
    # print >> f, dates_HW

    statis = []
    year_last = dates_HW[0][0].year - 1
    for i in dates_HW:
        year = i[0].year
        length = len(i)
        if year != year_last:
            for i in range(1, year - year_last):
                statis.append([year_last + i, 0, 0, 0, 0])
            array_t = [year, 1, length, length, length]
            statis.append(array_t)
            year_last = year
        else:
            statis[-1][1] += 1
            statis[-1][2] += length
            statis[-1][3] = max(statis[-1][3], length)
            statis[-1][4] = min(statis[-1][4], length)

    return np.array(statis)
