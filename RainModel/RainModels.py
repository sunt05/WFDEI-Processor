import numpy as np
from scipy.stats import norm

class LogNormalRainIntensity:
    # An instantaneous rainfall intensity distribution model
    def __init__(self, logmean, logsd):
        self.logMean = logmean
        self.logSd = logsd
    
    def transform(self, normalData):
        # Transform normally distributed data to log-normal of the given type
        bbb = np.exp(self.logMean + normalData*self.logSd)
        return bbb
    
    def calculateRainThreshold(self, p1):
        # Calculate quantile of p1 (0 to 1) based on parameters of log-normal
        return np.exp(norm.isf(p1, loc=self.logMean, scale=self.logSd))
        
class WeibullClusterExtent:
    # An model to realize cluster extent (temporal)
    def __init__(self, scale, shape):
        self.scale = scale
        self.shape = shape
        
    def realize(self, lowerLimit, upperLimit):
        # Sample from Weibull to get overall length of cluster
        # Parameters specified in minutes; output is in minutes
        valid = False 
        while not valid:
            extent = self.scale * np.random.weibull(self.shape, 1) # Minutes
            # Round to time resolution 
            valid = (extent >= lowerLimit) & (extent <= upperLimit)
        
        return extent