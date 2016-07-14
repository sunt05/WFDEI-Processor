import numpy as np
from rainHelper import generateNoise
class RainGenerator:
    # A rain cluster characterised by an Intensity model, Extent Model
    def __init__(self, Intensity, Extent):
        # Initialize Rain Generator. Inputs:
        # A rain Intensity object with a .transform() method (mm/hr)
        # A cluster extent object with a .realize() method (min)
        # lowerExtent and upperExtent are the time range (minutes) to aim for
        self.Intensity = Intensity  # Model of intensity dist
        self.Extent = Extent        # Model of extent dist
        self.durationLimits = None
        self.rainFraction = None
        self.timeRes = None
        self.extent = None
        self.beta = None
        self.metaGauss = None
        self.rainRate = None
        self.extentBins = None
           
    def setDurationLimits(self, lowerExtent, upperExtent):
        self.durationLimits = [lowerExtent, upperExtent]
        # Reset any previously realized values of event temporal extent
        self.extent = None
        self.extentBins = None
        
    def setSpectralPowerExpt(self, beta):
        self.beta = None
        if (beta > 0) | (beta < -3):
            raise ValueError('Spectral power exponent must be between 0 and -2: %f given' % (p1))
        self.beta = beta
        
    def setRainFraction(self, p1):   
        # Set mean fraction of cluster that contains rain (0 to 1)
        self.rainFraction = None
        if (p1 < 0) | (p1 > 1):
            raise ValueError('Rain fraction must range from 0 to 1. %f given' % (p1))
        self.rainFraction = p1

    def setTimeRes(self, timeRes):
        # Define time resolution in minutes
        self.timeRes = timeRes

        
    def getExtent(self):          
        return self.extent
    
    def realizeExtent(self):
        # Realize overall duration of rain event including intermittency
        self.extent = None
        self.extentBins = None
        
        if self.durationLimits is None:
            raise TypeError('Duration range must be set')
        
        # Durations shorter than time res not allowed
        lowerLimit = max(self.durationLimits[0], self.timeRes)
        self.extent = self.Extent.realize(lowerLimit, self.durationLimits[1])       # Exact extent in minutes
        self.extent = np.ceil(self.extent/self.timeRes)*self.timeRes                # Round to time resolution of new time series
        self.extentBins = int(self.extent/self.timeRes)                             # Work out length of vector in realized time series
        
    def realizeTimeSeries(self):
        self.metaGauss = None
        self.rainRate = None
        
        if self.extent is None:
            self.realizeExtent() # realize if it hasn't been called yet
       
        # Boundary condition: First and final time step of realized series must be wet.
        # Force this and assume equal intensity during both bins for short series
        if self.extentBins <= 1:
            self.rainRate = np.ones(self.extentBins)
            return()
        
        # Otherwise realize time series stochastically
        zeroThreshold = self.Intensity.calculateRainThreshold(self.rainFraction)

        # Take time series length to next greatest 2^n
        nt = 2**(np.ceil(np.log(self.extentBins)/np.log(2)))
        
        mask = np.zeros(nt)
        mask[0:self.extentBins] = np.ones(self.extentBins)
        
        validExtent = False
        while not validExtent:
            # Realize repeatedly until boundary conditions are met
            self.metaGauss = None
            self.rainRate = None
            metaGauss = generateNoise(self.beta, nt)
            
            # Cut back to nt length from its base-2 version
            self.metaGauss = metaGauss[mask > 0] 
            
            # Transform to rain intensity distrib and set no-rain threshold
            self.rainRate = self.Intensity.transform(self.metaGauss)
            self.rainRate[self.rainRate < zeroThreshold] = 0
            # Boundary condition: Is generated extent == the time aimed for?
            validExtent = ((np.ptp(np.where(self.rainRate > 0))+1) >= self.extentBins)
            
            
    def getTimeSeries(self):
        if self.rainRate is None:
            self.realizeTimeSeries()
            
            
        return self.rainRate
        

    
