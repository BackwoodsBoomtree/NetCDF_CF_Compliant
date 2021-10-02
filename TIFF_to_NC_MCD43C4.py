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
from dateutil.relativedelta import relativedelta
from osgeo import gdal
import numpy as np

# Input arguments 
input_dir   = '/mnt/g/MCD43C4/tif/Daily/0.05'
output_dir  = '/mnt/g/MCD43C4/nc/Daily/0.05'
output_name = 'MCD43C4.A'
vi_list     = 'EVI', 'NDVI', 'NIRv', 'LSWI'
interval    = 'Daily' # must be 'Daily', '8-day', or 'Monthly
res         = 0.05
extent      = -180, 180, -90, 90
year_list   = list(range(2018, 2020 + 1)) # Start and end year

def create_nc_obj(raster_list, output_name, vi, interval, res, extent, year):
    # Create NetCDF file    
    nco = Dataset(output_name, 'w', clobber = True, format = "NETCDF4")
    
    # Meta-data  
    nco.LongName        = " ".join(['MCD43C4', interval, vi, str(year), str(str(res) + '-degrees')])
    nco.ShortName       = "_".join(['MCD43C4', interval, vi, str(year), str(str(res) + '-degrees')])
    nco.GranuleID       = os.path.basename(output_name)
    nco.VersionID       = '1.0'
    nco.Format          = 'NetCDF4'
    nco.Conventions     = 'CF-1.9'
    nco.ProcessingLevel = 'Level 4'
    nco.Source          = 'MCD43C4 0.05-degree Daily product: https://e4ftl01.cr.usgs.gov/MOTA/MCD43C4.006/'
    nco.ProcessingCenter     = 'University of Oklahoma College of Atmospheric and Geographic Sciences'
    nco.ContactPersonName    = 'Doughty, Russell'
    nco.ContactPersonRole    = 'Research Scientist'
    nco.ContactPersonEmail   = 'russell.doughty@ou.edu'
    nco.ContactPersonAddress = "GeoCarb Mission, University of Oklahoma, Norman, OK, 73019, USA";
    nco.SouthernBoundingCoordinate = str(extent[2])
    nco.NorthernBoundingCoordinate = str(extent[3])
    nco.WesternBoundingCoordinate  = str(extent[0])
    nco.EasternBoundingCoordinate  = str(extent[1])
    nco.LatitudeResolution  = str(res)
    nco.LongitudeResolution = str(res)
    nco.ProductionDateTime  = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    nco.comment     = ('A strict snow filter of 0 was used to exclude data with any percentage of snow per Walther et al. (2016). ' +
                        'Data with a QC flag of 4 or 5 were excluded as there was very little data for tropics when 3 was excluded. ' +
                        'See the Readme.md at https://github.com/GeoCarb-OU/MCD43C4_VIs for more info on how this data was processed.')

    lon_list   = np.arange(extent[0], extent[1], res)
    lat_list   = np.flip(np.arange(extent[2] + res, extent[3] + res, res))
    
    # Create dimensions, variables and attributes
    nco.createDimension('lon', len(lon_list))
    nco.createDimension('lat', len(lat_list))
    nco.createDimension('time', None)
    
    # Time
    time               = nco.createVariable('time', 'f8', ('time',))
    time.standard_name = 'time'
    time.calendar      = 'gregorian'
    time.units         = 'days since %s-01-01' % str(year)
    if interval == 'Daily':
        dates              = [datetime(year, 1, 1) + n * timedelta(days = 1) for n in range(len(raster_list))] # list of dates    
    if interval == '8-day':
        dates              = [datetime(year, 1, 1) + n * timedelta(days = 8) for n in range(0, 46)] # list of dates
    if interval == 'Monthly':
        dates              = [datetime(year, 1, 1) + relativedelta(month = n) for n in range(1, 13)] # list of dates    
    
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
    crs.crs_wtk                     = 'GEOGCRS["WGS 84", DATUM["World Geodetic System 1984", ELLIPSOID["WGS 84",6378137,298.257223563, LENGTHUNIT["metre",1]]], PRIMEM["Greenwich",0, ANGLEUNIT["degree",0.0174532925199433]], CS[ellipsoidal,2], AXIS["geodetic latitude (Lat)",north, ORDER[1], ANGLEUNIT["degree",0.0174532925199433]], AXIS["geodetic longitude (Lon)",east, ORDER[2], ANGLEUNIT["degree",0.0174532925199433]], ID["EPSG",4326]]'
    
    # Variable
    var = nco.createVariable(vi, 'i4', ('time', 'lat', 'lon'), zlib = True, fill_value = -9999)
    var.long_name    = " ".join(['MCD43C4', vi, str(year), str(str(res) + '-degrees')])
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
        print('I have created the file: %s\n' % output_name)

    elif year == 2001:
        dummy             = gdal.Open(raster_list[0])
        dummy             = dummy.ReadAsArray()
        dummy[dummy > -1] = -9999            
        
        for i in range(len(raster_list)+1):
            time[i] = date2num(dates[i], units = time.units, calendar = time.calendar) # netCDF4 function that translates datetime object into proper format for .nc
            if i == 22:
                var[i, :, :] = dummy
            elif i < 22:
                layer = gdal.Open(raster_list[i])
                array = layer.ReadAsArray()  # data
                var[i, :, :] = array
            elif i > 22:
                layer = gdal.Open(raster_list[i - 1])
                array = layer.ReadAsArray()  # data
                var[i, :, :] = array
        nco.close()
        print('I have created the file: %s\n' % output_name)                    
    
    else:
        # For all years after 2001. Step through each raster for the year, writing time and data to NetCDF
        for i in range(len(raster_list)):
            time[i] = date2num(dates[i], units = time.units, calendar = time.calendar) # netCDF4 function that translates datetime object into proper format for .nc
            layer = gdal.Open(raster_list[i])
            array = layer.ReadAsArray()  # data
            var[i, :, :] = array
        nco.close()
        print('I have created the file: %s\n' % output_name)
        

def tiff_to_netcdf(input_dir, output_dir, output_name, vi, interval, res, year):
    
    raster_list = sorted([os.path.join(input_dir, l) for l in os.listdir(input_dir) if l.endswith('.tif')]) # All rasters in input_dir
    
    # Check whether the specified path exists or not
    output_dir = os.path.join(output_dir, vi)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Make sure there are enough files to package up (2000 and 2001 have short MODIS records)
    if interval == 'Daily':
        if len(raster_list) != 365 and len(raster_list) != 366 and year != 2000 and year != 2001: 
            quit('ERROR: There should be 365 or 365 input rasters for the year ' + str(year) + ', but there are only ' + str(len(raster_list)) + '.')                
    elif interval == '8-day':
        if len(raster_list) != 46 and year != 2000 and year != 2001: 
            quit('ERROR: There should be 46 input rasters for the year ' + str(year) + ', but there are only ' + str(len(raster_list)) + '.')
    elif interval == 'Monthly':
        if len(raster_list) != 12 and year != 2000 and year != 2001: 
            quit('ERROR: There should be 12 input rasters for the year ' + str(year) + ', but there are only ' + str(len(raster_list)) + '.')           

    # Printing for visual quality control
    print('Packing up these ' + str(len(raster_list)) + ' rasters into the .nc file for the year ' + str(year) + ' , starting with:')
    print(raster_list[0])

    # 8-day and Monthly TIFs are placed into a single nc file for the year,
    # and daily tifs are placed into monthly nc files
    if interval == '8-day' or interval == 'Monthly':
        out_name    = ".".join([output_name, str(year), vi, interval, str(res), 'nc'])
        out_file    = os.path.join(output_dir, out_name)
        create_nc_obj(raster_list, out_file, vi, interval, res, extent, year)
        
    elif interval == 'Daily':
        doy = 0
        for m in list(range(1, 13)):
            if m == 1 or m == 3 or m ==5 or m == 7 or m == 8 or m == 10 or m == 12:
                m_raster_list = raster_list[doy : (doy + 31)]
                out_name      = ".".join([output_name, str(year), str(m).zfill(2), vi, interval, str(res), 'nc'])
                out_file      = os.path.join(output_dir, out_name)
                doy           = doy + 31
                print(m_raster_list[0], m_raster_list[-1])
                create_nc_obj(m_raster_list, out_file, vi, interval, res, extent, year)
            elif m == 4 or m == 6 or m == 9 or m == 11:
                m_raster_list = raster_list[doy : (doy + 30)]
                out_name      = ".".join([output_name, str(year), str(m).zfill(2), vi, interval, str(res), 'nc'])
                out_file      = os.path.join(output_dir, out_name)
                doy           = doy + 30
                print(m_raster_list[0], m_raster_list[-1])
                create_nc_obj(m_raster_list, out_file, vi, interval, res, extent, year)
            elif m == 2:
                if len(raster_list) == 365:
                    m_raster_list = raster_list[doy : (doy + 28)]
                    out_name      = ".".join([output_name, str(year), str(m).zfill(2), vi, interval, str(res), 'nc'])
                    out_file      = os.path.join(output_dir, out_name)
                    doy           = doy + 28
                    print('Feb was not a leap year for ' + str(year))
                    print(m_raster_list[0], m_raster_list[-1])
                    create_nc_obj(m_raster_list, out_file, vi, interval, res, extent, year)
                if len(raster_list) == 366:
                    m_raster_list = raster_list[doy : (doy + 29)]
                    out_name      = ".".join([output_name, str(year), str(m).zfill(2), vi, interval, str(res), 'nc'])
                    out_file      = os.path.join(output_dir, out_name)
                    doy           = doy + 29
                    print('Feb was a leap year for ' + str(year))
                    print(m_raster_list[0], m_raster_list[-1])                    
                    create_nc_obj(m_raster_list, out_file, vi, interval, res, extent, year)


for v in range(len(vi_list)):          
    i_dir = os.path.join(input_dir, vi_list[v])
    sub_dirs      = [f.path for f in os.scandir(i_dir) if f.is_dir()] # Get subdirectories of input_dir
    filtered_dirs = []

    for i in range(len(year_list)): # filter dirs by the year list
        if str(year_list[i]) in sub_dirs[i]:
            filtered_dirs.append(sub_dirs[i])

    for j in range(len(filtered_dirs)):
        tiff_to_netcdf(filtered_dirs[j], output_dir, output_name, vi_list[v], interval, res, year_list[j])