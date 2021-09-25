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

# Input arguments for 8-day data
input_dir   = '/mnt/g/MCD43C4/tif/8-day/0.25/LSWI' # input folder contains TIFFs from all years
output_dir  = '/mnt/g/MCD43C4/nc/8-day/0.25/LSWI'
output_name = 'MCD43C4.A'
vi          = 'LSWI'
time        = '8-day' # must be 'Daily', '8-day', or 'Monthly
res         = '0.25'
year_list   = list(range(2018, 2020 + 1)) # Start and end year


def tiff_to_netcdf(input_dir, output_dir, output_name, vi, time, res, year):
    
    raster_list = sorted([os.path.join(input_dir, l) for l in os.listdir(input_dir) if l.endswith('.tif')]) # All rasters in input_dir
    
    # Check whether the specified path exists or not
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
             
        raster_list = raster_list
    
    if len(raster_list) != 46 and year != 2000 and year != 2001: 
        quit('ERROR: There should be 46 input rasters for the year ' + str(year) + ', but there are only ' + str(len(raster_list)) + '.')

    # Printing for visual quality control
    print('Packing up these ' + str(len(raster_list)) + ' rasters into the .nc file for the year ' + str(year) + ':')
    print("\n".join(raster_list))

    # Create NetCDF file
    out_name = ".".join([output_name, str(year), vi, time, res, 'nc'])
    out_file    = os.path.join(output_dir, out_name)
    nco         = Dataset(out_file, 'w', clobber = True, format = "NETCDF4")
    
    # Meta-data
    nco.Conventions = 'CF-1.9'
    nco.title       = " ".join(['MCD43C4', vi, str(year), str(res + '-degrees')])
    nco.institution = 'University of Oklahoma College of Atmospheric and Geographic Sciences'
    nco.source      = 'Original data found at: https://e4ftl01.cr.usgs.gov/MOTA/MCD43C4.006/'
    nco.references  = 'Schaaf, C., Wang, Z. (2015). MCD43C4 MODIS/Terra+Aqua BRDF/Albedo Nadir BRDF-Adjusted Ref Daily L3 Global 0.05Deg CMG V006 [Data set]. NASA EOSDIS Land Processes DAAC. Accessed 2021-09-25 from https://doi.org/10.5067/MODIS/MCD43C4.006'
    nco.history     = """[2021-09-25] File created."""
    nco.comment     = ('This data has been produced from the original 0.05-degree daily MCD43C4 files by Dr. Russell Doughty. ' + 
                        'A strict snow filter of 0 was used to exclude data with any percentage of snow per Walther et al. (2016). ' +
                        'Data flagged as 4 or 5 were excluded from the data as there was very little data for tropics when 3 was excluded. ' +
                        'See the Readme.md at https://github.com/GeoCarb-OU/MCD43C4_VIs for more info on how this data was processed. ' +
                        'Please cite Schaaf, C., Wang, Z. (2015) if you use this data in your project.')
    
    # Get shape, extent, coordinates from one of the input rasters
    ds         = gdal.Open(raster_list[0])
    a          = ds.ReadAsArray()
    nlat, nlon = np.shape(a)

    b        = ds.GetGeoTransform()  # bbox, interval
    lon_list = np.arange(nlon) * b[1] + b[0]
    lat_list = np.arange(nlat) * b[5] + b[3]

    # Create dimensions, variables and attributes
    nco.createDimension('lon', nlon)
    nco.createDimension('lat', nlat)
    nco.createDimension('time', None)
    
    # Time
    time               = nco.createVariable('time', 'f8', ('time',))
    time.standard_name = 'time'
    time.calendar      = 'gregorian'
    time.units         = 'days since %s-01-01' % str(year)
    dates              = [datetime(year, 1, 1) + n * timedelta(days = 8) for n in range(0, 46)] # list of dates
    
    # Lon
    lon               = nco.createVariable('lon', 'f4', ('lon',))
    lon.standard_name = 'longitude'
    lon.units         = 'degrees_east'
    
    # Lat
    lat               = nco.createVariable('lat', 'f4', ('lat',))
    lat.standard_name = 'latitude'
    lat.units         = 'degrees_north'

    # CRS - These values pulled from gdalinfo after converting GeoTiffs to WGS84
    crs                             = nco.createVariable('crs', 'i4')
    crs.long_name                   = 'World Geodetic System 1984 (WGS84)'
    crs.grid_mapping_name           = 'latitude_longitude'
    crs.longitude_of_prime_meridian = 0.0
    crs.semi_major_axis             = 6378137.0
    crs.inverse_flattening          = 298.257223563
    crs:crs_wtk                     = 'GEOGCRS["WGS 84", DATUM["World Geodetic System 1984", ELLIPSOID["WGS 84",6378137,298.257223563, LENGTHUNIT["metre",1]]], PRIMEM["Greenwich",0, ANGLEUNIT["degree",0.0174532925199433]], CS[ellipsoidal,2], AXIS["geodetic latitude (Lat)",north, ORDER[1], ANGLEUNIT["degree",0.0174532925199433]], AXIS["geodetic longitude (Lon)",east, ORDER[2], ANGLEUNIT["degree",0.0174532925199433]], ID["EPSG",4326]]'
    
    # Variable
    var = nco.createVariable(vi, 'i4', ('time', 'lat', 'lon'), zlib = True, fill_value = -9999)
    var.long_name    = " ".join(['MCD43C4', vi, str(year), str(res + '-degrees')])
    var.units        = 'Index'
    var.scale_factor = 0.0001
    var.add_offset   = 0.0
    var.grid_mapping = 'crs'
    var.set_auto_maskandscale(False)

    # Write lon,lat
    lon[:] = lon_list
    lat[:] = lat_list

    # Fill missing DOYs with dummy rasters for 2000 and 2001
    if year == 2000:
        dummy             = gdal.Open(raster_list[0])
        dummy             = dummy.ReadAsArray()
        dummy[dummy > -1] = -9999
        
        for i in range(len(raster_list) + 7):
            time[i] = date2num(dates[i], units = time.units, calendar = time.calendar) # netCDF4 function that translates datetime object into proper format for .nc
            if i > 6:
                layer = gdal.Open(raster_list[i - 7])
                array = layer.ReadAsArray()  # data
                var[i, :, :] = array
            else:
                var[i, :, :] = dummy
        nco.close()
        print('I have created the file: %s\n' % out_file)

    elif year == 2001:
        dummy             = gdal.Open(raster_list[0])
        dummy             = dummy.ReadAsArray()
        dummy[dummy > -1] = -9999            
        
        for i in range(len(raster_list)+1):
            time[i] = date2num(dates[i], units = time.units, calendar = time.calendar) # netCDF4 function that translates datetime object into proper format for .nc
            if i == 22:
                GPP[i, :, :] = dummy
            elif i < 22:
                layer = gdal.Open(raster_list[i])
                array = layer.ReadAsArray()  # data
                GPP[i, :, :] = array
            elif i > 22:
                layer = gdal.Open(raster_list[i - 1])
                array = layer.ReadAsArray()  # data
                GPP[i, :, :] = array
        nco.close()
        print('I have created the file: %s\n' % out_file)                    
    
    else:
        # For all years after 2001. Step through each raster for the year, writing time and data to NetCDF
        for i in range(len(raster_list)):
            time[i] = date2num(dates[i], units = time.units, calendar = time.calendar) # netCDF4 function that translates datetime object into proper format for .nc
            layer = gdal.Open(raster_list[i])
            array = layer.ReadAsArray()  # data
            var[i, :, :] = array
        nco.close()
        print('I have created the file: %s\n' % out_file)

sub_dirs      = [f.path for f in os.scandir(input_dir) if f.is_dir()] # Get subdirectories of input_dir
filtered_dirs = []

for i in range(len(year_list)): # filter dirs by the year list
    if str(year_list[i]) in sub_dirs[i]:
        filtered_dirs.append(sub_dirs[i])
        
for j in range(len(filtered_dirs)):
    tiff_to_netcdf(filtered_dirs[j], output_dir, output_name, vi, time, res, year_list[j])
    


