from WFDEI_Processor import *
from WFDEI_Processor.HW_Utility import HWFinder

import pandas as pd

data_out_1h=pd.read_csv('~/Downloads/WFDEI_2012_60min.txt',sep=' ')
HWFinder(data_out_1h['Tair'],123,122)
