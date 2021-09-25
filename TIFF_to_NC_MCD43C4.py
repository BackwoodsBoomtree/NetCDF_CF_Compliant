#-------------------------------------------------------------------------------
# Name:        Timeseries TIFF to netCDF
#
# Purpose:     Bundles tiff files into nc for MODIS MCD43C4.
#
# Author:      Russell Doughty, PhD
#
# Created:     September 25th, 2021
#-------------------------------------------------------------------------------

import os
from netCDF4 import Dataset,date2num
from datetime import datetime, timedelta
from osgeo import gdal
import numpy as np
