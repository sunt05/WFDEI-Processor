import numpy as np

# calculate saturation vapour pressure [Pa]
def evpsat(T_degC, p_hPa):
    ''' get saturation pressure [Pa] for a given air temperature [degC] and pressure [hPa]'''
    from numpy import log10

    T_degC = max(T_degC, 0.001, key=lambda x: abs(x))
    p_kPa = p_hPa / 10

    if T_degC > 0:
        # For T > 0 degC (f corrects for the fact that we are not dealing with
        # pure water)
        e_mb_pos = 6.1121 * np.exp(((18.678 - T_degC / 234.5)
                                    * T_degC) / (T_degC + 257.14))
        f_pos = 1.00072 + p_kPa * (3.2e-6 + 5.9e-10 * T_degC ** 2)
        esat_hPa_pos = e_mb_pos * f_pos
        esat_hPa = esat_hPa_pos
    else:
        # For T <= 0 degC
        e_mb_neg = 6.1115 * np.exp(((23.036 + T_degC / 333.7)
                                    * T_degC) / (T_degC + 279.82))
        f_neg = 1.00022 + p_kPa * (3.83e-6 + 6.4e-10 * T_degC ** 2)
        esat_hPa_neg = e_mb_neg * f_neg
        esat_hPa = esat_hPa_neg

    esat_Pa = esat_hPa * 100
    return esat_Pa


# calculate actual vapour pressure [Pa]
def eact(qv, p_kPa):
    # Specific gas constant for dry air  (Rv = R/RMM(dry air)) / J K^-1 kg^-1
    Rd = 287.04
    # Specific gas contstant for water vapour (Rv = R/RMM(H20)) / J K^-1 kg^-1
    Rv = 461.5
    # calculate actual vapour pressure [Pa]
    evp_Pa = Rv * qv * p_kPa * 1000 / (Rd + qv * (Rv - Rd))
    return evp_Pa


# convert specific humidity [kg/kg] to relative humidity [%]
def q2rh(qv, p_kPa, T_degC):
    eact_Pa = eact(qv, p_kPa)
    esat_Pa = evpsat(T_degC, p_kPa * 10)
    rh_pct = 100 * eact_Pa / esat_Pa

    return rh_pct


# vectorize q2rh
def vq2rh(qv, p_kPa, T_degC):
    eact_Pa = eact(qv, p_kPa)
    vevpsat = np.vectorize(evpsat)
    esat_Pa = vevpsat(T_degC, p_kPa * 10)
    rh_pct = 100 * eact_Pa / esat_Pa

    return rh_pct

# functions for calculating RH
# from package of meteo
Mw = 18.0160  # molecular weight of water
Md = 28.9660  # molecular weight of dry air


def esat(T):
    ''' get sateration pressure (units [Pa]) for a given air temperature (units [K])'''
    from numpy import log10
    TK = 273.15
    e1 = 101325.0
    logTTK = log10(T / TK)
    x1 = 10.79586 * (1 - TK / T)
    x2 = 5.02808 * logTTK
    x3 = 1.50474 * 1e-4 * (1. - 10**(-8.29692 * (T / TK - 1)))
    x4 = 0.42873 * 1e-3 * (10**(4.76955 * (1 - TK / T)) - 1) - 2.2195983
    xx = x1 - x2 + x3 + x4
    esat = e1 * 10 ** xx
    return esat


def sh2mixr(qv):
    '''conversion from specific humidity (units [kg/kg]) to mixing ratio (units also [kg/kg])'''
    return qv / (1. - qv)


def mixr2rh(mixr, p, T):
    '''purpose: conversion mixing ratio to relative humidity [kg/kg] (not tested)'''
    return mixr * p / ((mixr + Mw / Md) * esat(T))


def sh2rh(qv, p, T):
    '''conversion from specific humidity (units [kg/kg]) to relative humidity in percentage'''
    return mixr2rh(sh2mixr(qv), p, T)
