# Helper methods for RainfallGenerator.py
import numpy as np
import numpy.fft as fft

def generateNoise(beta, nt):
    # Generates fractal noise (metagaussian) given beta (spectral power exponent, nt (length of series, must be 2^n (int n))
    if beta is None:
        raise TypeError('Spectral power exponent is not set')       
    
    st = np.abs(beta)
    kt = np.concatenate((np.arange(0, nt/2+1), np.arange(-nt/2+1, 0)))
    kt = kt**2
    kt[0]=0.000001
    kt=kt**(-st/4)
    kt = kt/np.sqrt(np.sum(np.abs(kt.flat)**2)) * nt
    
    # Metagauss
    ph = np.random.randn(nt)
    ph = fft.fftn(ph)
    ph = ph/np.abs(ph)
    ph[0] = 0
    ph = kt*ph
    metaGauss = fft.ifftn(ph).real # Meta-gaussian time series a la RainFARM
    # Force the metaGaussian noise to be unit variance and zero-centred
    metaGauss = (metaGauss - np.mean(metaGauss))/np.std(metaGauss)
    return metaGauss

def rleTotals(vec):
    # Return start & length non-zero runs in an array
    vals = vec
    vec = vec > 0
    d = np.hstack(([0], vec, [0]))
    difs = np.diff(d)
    run_starts, = np.where(difs > 0)
    run_ends, = np.where(difs < 0)
    results = {}
    results['start'] = run_starts
    results['length'] = run_ends-run_starts
    results['total'] = [sum(vals[run_starts[i]:run_ends[i]]) for i in np.arange(run_starts.shape[0])]
    return results
    
   
def realizeForEvent(rain,  minLength, maxLength, timeStepsPerInputBin):
    # Realize each rain event in the coarse time series
    # Checks there's no gap within the event longer than timeStepsPerInputBin
    rain.setDurationLimits(minLength, maxLength)
    validGapsInSeries = False
    while validGapsInSeries == False:
        rain.realizeTimeSeries()
        a = rain.getTimeSeries()
        gapLengths = rleTotals(a == 0.0)
        validity = np.array(gapLengths['length']) < timeStepsPerInputBin
        validGapsInSeries = sum(validity) == validity.shape[0]
    
    return a
    
    
def generateTestCoarseData():

    # Generates hourly data as a coarse met time series for input to SuewsDataProcessing to5min method
    # Cols: YYYY, DOY, H, M, <21 cols of random numbers>
    # result[:,13] is rainfall accumulation so has intentional gaps in it
    
    # Generate whole year of days
    Y = 2015
    DOYS = np.arange(0, 365)+1 # DOY 1 = 1st jan
    H = np.arange(0,24)
    nMin = 1# Minutes per hour (1 = hourly data)
    M = np.arange(0,nMin)
    nRandomCols = 21 # Number of columns of data to generate
    
    # Make a day of minutes
    dayOfMins =np.tile(M, 24)
    # ... and the hour to which they belong
    hoursEachMin = np.tile(H, nMin)
    # Make a year of minutes
    yearOfMins = np.tile(dayOfMins, 365)
    yearOfHours = np.tile(hoursEachMin, 365)
    yearOfYears = np.tile(Y, yearOfMins.shape[0])
    yearOfDOYS = np.repeat(DOYS, hoursEachMin.shape[0])
    randomColumn = np.random.uniform(0,1, (yearOfHours.shape[0], nRandomCols))
    # Put it all together
    result = np.empty([yearOfMins.shape[0], nRandomCols+4])
    result[:,0] = yearOfYears
    result[:,1] = yearOfDOYS
    result[:,2] = yearOfHours
    result[:,3] = yearOfMins
    result[:,4:nRandomCols+4] = randomColumn
    # Manipulate rain to be mostly below a threshold so it's a bit bursty
    result[result[:,13] < 0.8, 13] = 0 # Only allow 20% of data through

    return result
    

def downscaleTimeSeries(series, res_in, res_out, IntensityModel, TemporalExtentModel, power_expt, wet_fraction):
    from RainGenerator import RainGenerator
    # Downscales a 1D vector series of rainfall accumulation
    # Inputs: series:               np.array of rain accumulations. Assumption: accumulation at Ti is accumulation between T(i-1) and T(i)
    #         res_in:               Resolution of input time series (minutes)
    #         res_out:              Desired resolution (minutes)
    #         IntensityModel:       Statistical model of rainfall intensity with a .transform(data) method that takes normally distributed input data, and 
    #                               .calculateRainThreshold(wet_fraction) model that works out where to place a rain/no rain threshold. 
    #         TemporalExtentModel:  Statistical model of rain period duration (pulses of rain) with .realize(lowerLimit, upperLimit) method. 
    #         power_expt:           Rain power spectrum exponent (power ~ freq^power_expt)
    #         wet_fraction:         The mean amount of time during a rain cluster that it's really raining (because intermittency exists)
    if sum(series<0) > 0:
        raise ValueError('Rainfall time series cannot contain negative data')
    
    if res_out < 1:
        raise ValueError('This technique is only reasonable for timescales down to ~1 min')
        
    if divmod(res_in, res_out)[1] > 0:
        raise ValueError('res_in must be an integer multiple of res_out: %f and %f given' % (res_in, res_out))

    timeConversion = res_in/res_out
    # Preallocate the time series 
    newTS = np.zeros(series.shape[0]*timeConversion)
     
    # Oversample original time series to get rainfall proportion in each coarse time bin
    osample = series[np.repeat(np.arange(0,series.shape[0]),timeConversion)]
        
    # Step through series and identify runs of non-zero rain accumulation
    rainOccurring = rleTotals(series)
    rainDurations = rainOccurring['length']*res_in 
    # Instantiate rainfall generator 
    rain = RainGenerator(IntensityModel, TemporalExtentModel) 
    rain.setTimeRes(res_out) 
    rain.setSpectralPowerExpt(power_expt)
    rain.setRainFraction(wet_fraction)
    
    # Realize a time series for each rain event (doesn't conserve mass, just duration)
    realz = [realizeForEvent(rain, rainDurations[i]-res_in+res_out, rainDurations[i], timeConversion) for i in range(rainDurations.shape[0])]

    # Apply  per-coarse-time-bin mass conservation.
    for i, ts in enumerate(realz): # Step through each realized event
        ots = ts
        for j, mass in enumerate(np.arange(0,rainOccurring['length'][i])):
            # Each coarse time bin has a rainfall depth that must be conserved
            massIndex = int(rainOccurring['start'][i]+j)
            mass = series[massIndex]

            # Find relevant part of realized time series
            startIndex = int(j*timeConversion) # Start of coarse bin in new time series
            endIndex = startIndex + min(ts.shape[0]-startIndex, int(timeConversion))
            if startIndex == endIndex:
                replaceRange = startIndex
            
            else:
                replaceRange = np.arange(startIndex, endIndex)
            
            ts[replaceRange] = ts[replaceRange]/sum(ts[replaceRange]) * mass
        # Populate new time series
        startPoint = int(rainOccurring['start'][i]*timeConversion)
        insertRange = np.arange(startPoint, startPoint+ts.shape[0])         
        newTS[insertRange] = ts

    return newTS
