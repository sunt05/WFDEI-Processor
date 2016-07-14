##########################################################################
# running section:
# provide parameters here
##########################################################################
from WFDEI_Processor.WFDEI_Processor import *
import time


# read in data from WFDEI files

input_path = '/Users/sunt05/Documents/Data/WFDEI/'
# input_path = '/Volumes/DATA-TS/WFDEI/'

output_path = '~/Downloads/'

year_start, year_end = 2012, 2012

lat, lon = 51.51, -0.12  # London
# lat, lon = 51.58, -1.8  # Swindon

tstep_min = 60  # minutes

# 0 for unven distribution, 1 for spectrum-based method
rain_opt = 0

start = time.time()
# write_SUEWS_forcing_1h(input_path, output_path,
#                        year_start, year_end, lat, lon,rain_opt)

write_SUEWS_forcing_tmin(input_path, output_path,
                         year_start, year_end, lat, lon, tstep_min, rain_opt)

end = time.time()
print('time used in processing:' + '%.2f' % (end - start) + ' s')
