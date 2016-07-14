##########################################################################
# WFDEI Interpolator
##########################################################################
# Purpose:
# 1. interpolate WFDEI data at hourly scale;
# 2. prepare input data for SUEWS.
##########################################################################
# Authors:
# Lingbo Xue, Uni. of Reading, L.Xue@student.reading.ac.uk
# Ting Sun, Uni. of Reading, Ting.Sun@reading.ac.uk
##########################################################################
# History:
# 20160610 LX: initial version.
# 20160618 TS: refactored with pandas.
# 20160619 TS & LX: solar related problems fixed.
# 20160706 TS & LX: solar parts replaced with Pysolar.
# 20160707 TS: relative humidity modified to be consistent with SUEWS.
# 20160707 TS: interactive input implemented.
# 20160713 TS: spectrum-based rain model is included.
# 20160713 TS: time step can be specified.
##########################################################################
# To do:
# 1. adapt this code for WATCH data so as to make this code a generic
# interpolaor.
# 2. incorporate heat wave finder
##########################################################################

# preload packages
import numpy as np
import os
import sys
from datetime import datetime, timedelta, date
from scipy import interpolate
import time
import pandas as pd
from Pysolar import solar
from scipy import signal
from RainModel.rainHelper import downscaleTimeSeries
from RainModel.RainModels import LogNormalRainIntensity, WeibullClusterExtent
import q2RH
import WATCH_Utility


# correct SW data based on energy conservation according to raw data:
def correct_SW(SW_proc, SW_raw, lat, lon):
    # copy SW_proc
    SW_crt = SW_proc.copy()

    # correct SW based on sunrise/sunset
    sol_elev = np.array([solar.GetAltitudeFast(
        lat, lon, t - timedelta(minutes=30)) for t in SW_proc.index])
    sol_elev_reset = np.sin(np.radians(sol_elev.copy()))
    sol_elev_reset[sol_elev > 0] = 1
    sol_elev_reset[sol_elev <= 0] = 0

    # SWdown correction:
    # force nocturnal values to zero:
    SW_crt = sol_elev_reset * SW_proc
    # rescale based on 3-hourly values:
    avg_raw = SW_raw.resample('D').mean()
    avg_crt = SW_crt.resample('D').mean()
    ratio = (avg_raw / avg_crt).reindex(index=SW_crt.index).fillna(method='pad')
    SW_crt = (ratio * SW_crt).fillna(method='pad')

    return SW_crt


# spectrum-based downscaling procedure to process raw rainfall data to
# tstep-min rainfall
def rain_model_AG(rain_raw, tstep_min):
    """

    Parameters
    ----------
    rain_raw : raw rainfall data
    tstep0_min : timestep of rain_raw in minute
    tstep_min : timestep to downscale at in minute

    Returns : tstep-min rainfall data in pandas series with index of time series
    -------

    """
    # get original timestep
    ix = rain_raw.index
    tstep0 = ix[1] - ix[0]
    # convert to min
    tstep0_min = tstep0.total_seconds() / 60

    rainIntensityLogMean = -1.103937465
    rainIntensityLogSd = 1.365991286
    rainExtentScale = 218.93547799
    rainExtentShape = 0.64986935
    rainExponent = -1.01324
    rainWetProb = 0.5236637

    intensityModel = LogNormalRainIntensity(
        rainIntensityLogMean, rainIntensityLogSd)
    extentModel = WeibullClusterExtent(rainExtentScale, rainExtentShape)
#     rainfallModel['powerLawExponent'] = rain_nml['RainParameters']['rainExponent']
    probabilityOfRainInsideCluster = rainWetProb

    # 60 is original tstep (min)
    realz = downscaleTimeSeries(rain_raw, tstep0_min, tstep_min, intensityModel,
                                extentModel, rainExponent, probabilityOfRainInsideCluster)

    # repack to pandas DataFrame
    ix = pd.date_range(rain_raw.index[0], freq='%dMin' %
                       tstep_min, periods=len(realz))
    rain_pd = pd.Series(realz, index=ix, name=rain_raw.name)

    # return tstep-min rainfall data in pandas series with index of time series
    return rain_pd


# process raw rain data with different models to tstep_min:
# the output will be in mm s-1
def process_Rain(rain_raw, tstep_min, opt):
    # get original timestep
    ix = rain_raw.index
    tstep0 = ix[1] - ix[0]
    # make the date range cover whole periods
    ix_new = pd.date_range(ix[0], ix[-1] + tstep0, freq=tstep0)
    # expand ending part
    rain_new = rain_raw.copy().reindex(index=ix_new)

    # even distributio over the desired timestep:
    if opt == 0:
        rain_proc = rain_new.resample(
            '%d min' % tstep_min).fillna(method='bfill')

    # spectrum-based downscaling:
    if opt == 1:
        rain_proc = rain_model_AG(rain_raw, tstep_min)

    # make the date range cover whole periods
    ix_new = pd.date_range(ix[0], ix[-1] + tstep0,
                           freq='%0.2d min' % tstep_min)
    # reindex to conform the processed series to desired coverage
    rain_proc = rain_proc.reindex(index=ix_new)
    rain_proc = rain_proc.fillna(value=0.)
    return rain_proc


# interpolate 3-hourly raw data to hourly results for SUEWS
def process_SUEWS_forcing_1h(data_raw_3h, lat, lon, rain_opt):
    #     print('*************** WFDEI Data Processor *************** ')
    # expand over the whole date range at the scale of 1 h
    ix = pd.date_range(data_raw_3h.index[
        0], data_raw_3h.index[-1] + timedelta(hours=3), freq="H")
    data_raw_1h = data_raw_3h.reindex(index=ix).resample('1h').mean()

    # create space for processed data
    data_proc_1h = data_raw_1h.copy()

    # Take off 30-min so that values represent average over previous hour
    # sol_elev = np.array([Astral.solar_elevation(
    # a, t - timedelta(minutes=30), lat, lon) for t in
    # data_raw_1h["SWdown"].index])

    sol_elev = np.array([solar.GetAltitudeFast(
        lat, lon, t - timedelta(minutes=30)) for t in data_raw_1h["SWdown"].index])

    sol_elev_reset = np.sin(np.radians(sol_elev.copy()))
    sol_elev_reset[sol_elev > 0] = 1
    sol_elev_reset[sol_elev <= 0] = 0

    # normal interpolation for instantaneous variables:
    # i.e., Tair, Wind, Psurf, Qair.
    # these variables have been processed
    var_List_inst = ["Tair", "PSurf", "Wind", "Qair"]
    for var in var_List_inst:
        data_proc_1h[var] = data_raw_1h[var].interpolate(method='polynomial', order=1).rolling(
            window=2, center=False).mean().fillna(method='pad')
        # fill the first three hours with value at 3:00 h
        data_proc_1h[var][0:4] = data_proc_1h[var][3]

    # normal interpolation for variables averaged over previous 3 hours:
    # i.e., SWdown, LWdown, Rainf
    var_List_p3h = ["SWdown", "LWdown", "Rainf"]

    # SWdown, LWdown:
    for var in var_List_p3h[:-1]:
        # convert to 30-min instantaneous values
        x0 = data_raw_1h[var].resample('30min').mean(
        ).interpolate(method='polynomial', order=1)
        # shift to get delta/6, so that values at
        # [t-1,t,t+1]=xt+delta*[1,3,5]/6
        data_proc_1h[var] = x0.shift(-3)[::2]

    # SWdown correction:
    data_proc_1h["SWdown"] = correct_SW(
        data_proc_1h["SWdown"], data_raw_3h["SWdown"], lat, lon)
    # # force nocturnal values to zero:
    # data_proc_1h["SWdown"] = sol_elev_reset * data_proc_1h["SWdown"]
    # # rescale based on 3-hourly values:
    # avg_3h_SWdown = data_raw_3h["SWdown"].resample('D').mean()
    # avg_1h_SWdown = data_proc_1h["SWdown"].resample('D').mean()
    # ratio_SWdown = (avg_3h_SWdown / avg_1h_SWdown).reindex(
    #     index=ix).resample('1h').mean().fillna(method='pad')
    # data_proc_1h["SWdown"] = (
    #     ratio_SWdown * data_proc_1h['SWdown']).fillna(method='pad')

    # Rainf: evenly distribute over the 3-h period
    # data_proc_1h['Rainf'] = data_raw_1h['Rainf'].interpolate(
    #     method='polynomial', order=0).shift(-2).fillna(method='pad', limit=2).fillna(value=0)
    # Rainf:
    # option 0: even distributio over the desired timestep:
    # option 1: spectrum-based downscaling:
    data_proc_1h['Rainf'] = process_Rain(data_raw_3h['Rainf'], 60, rain_opt)

    # export processed data
    header = ["iy", "id", "it", "imin", "qn", "qh", "qe", "qs", "qf", "U", "RH", "Tair", "pres",
              "rain", "kdown", "snow", "ldown", "fcld", "wuh", "xsmd", "lai", "kdiff", "kdir", "wdir"]
    data_out_1h = pd.DataFrame(index=data_proc_1h.index, columns=header)
    var_out_list = ['SWdown', 'LWdown', 'Rainf', 'Tair', 'PSurf', 'Wind']
    var_out_zip = np.array(
        [var_out_list, ['kdown', 'ldown', 'rain', 'Tair', 'pres', 'U']]).T

    # refill starting & ending missing values
    var_fill_list = ['SWdown', 'LWdown', 'Tair', 'PSurf', 'Wind', 'Qair']
    for var in var_fill_list:
        data_proc_1h[var][:4] = data_proc_1h[var][4]
        data_proc_1h[var] = data_proc_1h[var].fillna(method='pad')

    # fill in variables:
    for p in var_out_zip:
        data_out_1h[p[1]] = data_proc_1h[p[0]]

    # RH calculation:
    data_out_1h['RH'] = vq2rh(data_proc_1h['Qair'],
                              data_proc_1h['PSurf'] / 1000,
                              data_proc_1h['Tair'] - 273.15)

    # unit conversion:
    # Tair: K -> degC
    data_out_1h['Tair'] -= 273.15
    # rainfall: kg m-2 -> mm  60 x 60 s / 1000 kg m-3 * 1000 mm m-1
    data_out_1h['rain'] *= 60 * 60
    # presure: Pa -> kPa
    data_out_1h['pres'] /= 1000

    # process timestamps
    data_out_1h['iy'] = data_proc_1h.index.year
    data_out_1h['id'] = data_proc_1h.index.dayofyear
    data_out_1h['it'] = data_proc_1h.index.hour
    data_out_1h['imin'] = data_proc_1h.index.minute

    # replace nan with -999
    data_out_1h = data_out_1h.fillna(value=-999)

    print('*************** WFDEI Data Successfully Processed *************** ')

    return data_out_1h


# interpolate 3-hourly raw data to t-minute results for SUEWS
def process_SUEWS_forcing_tmin(data_raw_3h, lat, lon, tstep_min, rain_opt):
    #     print('*************** WFDEI Data Processor *************** ')
    tstep = timedelta(minutes=tstep_min)
    tstart_raw, t0_raw, tend_raw = data_raw_3h.index[[0, 1, -1]]

    # expand over the whole date range at the scale of 1 h
    ix = pd.date_range(tstart_raw, tend_raw + timedelta(hours=3), freq=tstep)
    data_raw_tmin = data_raw_3h.reindex(index=ix).resample(tstep).mean()

    # create space for processed data
    data_proc_tmin = data_raw_tmin.copy()

    # Take off 30-min so that values represent average over previous hour
    # sol_elev = np.array([Astral.solar_elevation(
    # a, t - timedelta(minutes=30), lat, lon) for t in
    # data_raw_tmin["SWdown"].index])

    # sol_elev = np.array([solar.GetAltitudeFast(
    #     lat, lon, t - timedelta(minutes=30)) for t in data_raw_tmin["SWdown"].index])
    #
    # sol_elev_reset = np.sin(np.radians(sol_elev.copy()))
    # sol_elev_reset[sol_elev > 0] = 1
    # sol_elev_reset[sol_elev <= 0] = 0

    # normal interpolation for instantaneous variables:
    # i.e., Tair, Wind, Psurf, Qair.
    # these variables have been processed
    var_List_inst = ["Tair", "PSurf", "Wind", "Qair"]
    for var in var_List_inst:
        data_proc_tmin[var] = data_raw_tmin[var].interpolate(method='polynomial', order=1).rolling(
            window=2, center=False).mean().fillna(method='pad')
        # fill the first three hours with value at 3:00 h
        data_proc_tmin[var][0:4] = data_proc_tmin[var][3]

    # normal interpolation for variables averaged over previous 3 hours:
    # i.e., SWdown, LWdown, Rainf
    var_List_p3h = ["SWdown", "LWdown", "Rainf"]

    # SWdown, LWdown:
    for var in var_List_p3h[:-1]:
        # convert to instantaneous values: put the values back to tstep / 2
        # positions
        x0 = data_raw_tmin[var].resample(
            timedelta(hours=3 / 2)).mean().shift(-1)
        # resample to a finer timestep
        x0 = x0.resample(timedelta(minutes=1)).mean()
        # interpolate to fill nan values
        x0 = x0.interpolate(method='polynomial', order=1)
        # resample to desired timestep
        x0 = x0.resample(tstep).mean()
        # average over the previous periods
        x0 = x0.rolling(window=2, center=False).mean().fillna(method='pad')
        # shift to get delta/6, so that values at
        data_proc_tmin[var] = x0

    # SWdown correction:
    data_proc_tmin["SWdown"] = correct_SW(
        data_proc_tmin["SWdown"], data_raw_3h["SWdown"], lat, lon)
    # # force nocturnal values to zero:
    # data_proc_tmin["SWdown"] = sol_elev_reset * data_proc_tmin["SWdown"]
    # # rescale based on 3-hourly values:
    # avg_3h_SWdown = data_raw_3h["SWdown"].resample('D').mean()
    # avg_tmin_SWdown = data_proc_tmin["SWdown"].resample('D').mean()
    # ratio_SWdown = (avg_3h_SWdown / avg_tmin_SWdown).reindex(
    #     index=ix).resample('1h').mean().fillna(method='pad')
    # data_proc_tmin["SWdown"] = (
    #     ratio_SWdown * data_proc_tmin['SWdown']).fillna(method='pad')

    # Rainf:
    # option 0: even distributio over the desired timestep;
    # option 1: spectrum-based downscaling.
    data_proc_tmin['Rainf'] = process_Rain(
        data_raw_3h['Rainf'], tstep_min, rain_opt)

    # export processed data
    header = ["iy", "id", "it", "imin", "qn", "qh", "qe", "qs", "qf", "U", "RH", "Tair", "pres",
              "rain", "kdown", "snow", "ldown", "fcld", "wuh", "xsmd", "lai", "kdiff", "kdir", "wdir"]
    data_out_tmin = pd.DataFrame(index=data_proc_tmin.index, columns=header)
    var_out_list = ['SWdown', 'LWdown', 'Rainf', 'Tair', 'PSurf', 'Wind']
    var_out_zip = np.array(
        [var_out_list, ['kdown', 'ldown', 'rain', 'Tair', 'pres', 'U']]).T

    # refill starting & ending missing values
    var_fill_list = ['SWdown', 'LWdown', 'Tair', 'PSurf', 'Wind', 'Qair']
    for var in var_fill_list:
        data_proc_tmin[var][:t0_raw] = data_proc_tmin[var][t0_raw]
        data_proc_tmin[var] = data_proc_tmin[var].fillna(method='pad')

    # fill in variables:
    for p in var_out_zip:
        data_out_tmin[p[1]] = data_proc_tmin[p[0]]

    # RH calculation:
    data_out_tmin['RH'] = vq2rh(data_proc_tmin['Qair'],
                                data_proc_tmin['PSurf'] / 1000,
                                data_proc_tmin['Tair'] - 273.15)

    # unit conversion:
    # Tair: K -> degC
    data_out_tmin['Tair'] -= 273.15
    # rainfall: kg m-2 -> mm  60 x 60 s / 1000 kg m-3 * 1000 mm m-1
    data_out_tmin['rain'] *= 60 * 60
    # presure: Pa -> kPa
    data_out_tmin['pres'] /= 1000

    # process timestamps
    data_out_tmin['iy'] = data_proc_tmin.index.year
    data_out_tmin['id'] = data_proc_tmin.index.dayofyear
    data_out_tmin['it'] = data_proc_tmin.index.hour
    data_out_tmin['imin'] = data_proc_tmin.index.minute

    # replace nan with -999
    data_out_tmin = data_out_tmin.fillna(value=-999)

    print('*************** WFDEI Data Successfully Processed *************** ')

    return data_out_tmin


def write_SUEWS_forcing_1h(WFDEI_path, output_path, year_start, year_end, lat, lon, rain_opt):
    # load raw 3-hourly data
    data_raw_3h = load_WFDEI_3h(input_path, year_start, year_end, lat, lon)

    # process raw data to hourly forcings for SUEWS
    data_out_1h = process_SUEWS_forcing_1h(data_raw_3h, lat, lon, rain_opt)

    # output files of each year
    print('output files:')
    for year in range(year_start, year_end + 1):
        data_out_1h_year = data_out_1h[
            datetime(year, 1, 1) + timedelta(minutes=60):datetime(year + 1, 1, 1)]
        file_output_year = os.path.expanduser(
            os.path.join(output_path, 'WFDEI_' + str(year) + '.txt'))
        data_out_1h_year.to_csv(file_output_year, sep=" ",
                                index=False, float_format='%.4f')
        print(file_output_year)

    print('********* WFDEI Data Processing Successfully Finished *********')


def write_SUEWS_forcing_tmin(WFDEI_path, output_path, year_start, year_end, lat, lon, tstep_min, rain_opt):
    # load raw 3-hourly data
    data_raw_3h = load_WFDEI_3h(input_path, year_start, year_end, lat, lon)

    # process raw data to hourly forcings for SUEWS
    data_out_tmin = process_SUEWS_forcing_tmin(
        data_raw_3h, lat, lon, tstep_min, rain_opt)

    # output files of each year
    print('output files:')
    for year in range(year_start, year_end + 1):
        data_out_tmin_year = data_out_tmin[
            datetime(year, 1, 1) + timedelta(minutes=tstep_min):datetime(year + 1, 1, 1)]
        file_output_year = os.path.expanduser(
            os.path.join(output_path, 'WFDEI_' + str(year) + '_' + str(tstep_min) + 'min.txt'))
        data_out_tmin_year.to_csv(file_output_year, sep=" ",
                                  index=False, float_format='%.4f')
        print(file_output_year)

    print('********* WFDEI Data Processing Successfully Finished *********')


##########################################################################
# running section:
# provide parameters here
##########################################################################
# read in data from WFDEI files

input_path = '/Users/sunt05/Documents/Data/WFDEI/'
# input_path = '/Volumes/DATA-TS/WFDEI/'

output_path = '~/Downloads/'

year_start, year_end = 2012, 2012

lat, lon = 51.51, -0.12  # London
# lat, lon = 51.58, -1.8  # Swindon

tstep_min = 60  # minutes

# 0 for unven distribution, 1 for spectrum-based method
rain_opt = 1

start = time.time()
# write_SUEWS_forcing_1h(input_path, output_path,
#                        year_start, year_end, lat, lon,rain_opt)

write_SUEWS_forcing_tmin(input_path, output_path,
                         year_start, year_end, lat, lon, tstep_min, rain_opt)

end = time.time()
print('time used in processing:' + '%.2f' % (end - start) + ' s')


# # interatively get input parameters
# while True:
#     # WFDEI path:
#     while True:
#         input_path = raw_input("Please input the path for WFDEI data: ")
#         input_path = os.path.realpath(os.path.expanduser(input_path))
#         print input_path
#         if os.path.lexists(input_path):
#             break
#         else:
#             print "No such directory. Try again..."
#
#     # output path:
#     while True:
#         output_path = raw_input("Please input the path for output: ")
#         output_path = os.path.realpath(os.path.expanduser(output_path))
#         print output_path
#         if os.path.lexists(output_path):
#             break
#         else:
#             print "No such directory. Try again..."
#
#     # year range:
#     while True:
#         year_start = int(raw_input(
#             "Please input the start year (YYYY): "))
#         year_end = int(raw_input(
#             "Please input the end year (YYYY): "))
#         print(1979 <= year_start <= year_end <= 2014)
#         if 1979 <= year_start <= year_end <= 2014:
#             break
#         else:
#             print "Please input valid years. Try again..."
#
#     # coordinates:
#     while True:
#         lat = float(raw_input(
#             "Please input the latitude (in deg): "))
#         lon = float(raw_input(
#             "Please input the longitude (in deg): "))
#         print(-90 < lat < 90 and -180 < lon < 180)
#         if -90 < lat < 90 and -180 < lon < 180:
#             break
#         else:
#             print "Please input valid coordinates. Try again..."
#
#     start = time.time()
#     write_SUEWS_forcing_1h(input_path, output_path,
#                            year_start, year_end, lat, lon)
#     end = time.time()
#     print('time used in processing:' + '%.2f' % (end - start) + ' s')
#
#     t = raw_input('Do you want to quit? Y/N')
#     if t == 'Y' or t == 'y':
#         ftp.quit()
#         break
