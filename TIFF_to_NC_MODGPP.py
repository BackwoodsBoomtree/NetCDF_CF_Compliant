#-------------------------------------------------------------------------------
# Name:        Timeseries TIFF to netCDF
#
# Purpose:     Bundles tiff files into nc for MODIS GPP.
#
# Author:      Russell Doughty, PhD
#
# Created:     October 23, 2020
#-------------------------------------------------------------------------------

import os
from netCDF4 import Dataset,date2num
from datetime import datetime, timedelta
from osgeo import gdal
import numpy as np

# Input arguments for 8-day data
input_dir = 'C:/Russell/Projects/Geometry/Data/oco3/8-day' # input folder contains TIFFs from all years
output_dir = 'C:/Russell/Projects/Geometry/Data/oco3/8-day_NC'

years = range(2015,2018+1) # Start and end year

def tiff_to_netcdf(input_dir, output_dir, years):
    
    raster_list = sorted([os.path.join(input_dir, l) for l in os.listdir(input_dir) if l.endswith('.tif')]) # All rasters in input_dir

    dummy = gdal.Open(raster_list[0])
    dummy = dummy.ReadAsArray()
    dummy[dummy > -1] = -9999

    for year in range(len(years)):
        
        raster_list_year = [r for r in raster_list if '.'+str(years[year]) in r] # List of rasters for the year
        
        if len(raster_list_year) != 46 and years[year] != 2000 and years[year] != 2001: 
            quit('ERROR: There should be 46 input rasters for the year ' + str(years[year]) + ', but there are only ' + str(len(raster_list_year)) + '.')
 
        # Printing for visual quality control
        print('Packing up these '+str(len(raster_list_year))+' rasters into the .nc file for the year '+str(years[year])+':')
        print("\n".join(raster_list_year))

        # Create NetCDF file
        out_file = os.path.join(output_dir,'GPP.MOD17A2H.v006.'+str(years[year])+'.8-day.CMG_0.05.nc')
        
        nco = Dataset(out_file, 'w', clobber=True, format="NETCDF4")
        
        # Meta-data
        nco.Conventions = 'CF-1.7'
        nco.title = 'MODIS Gross Primary Production'
        nco.institution = 'California Institute of Technology'
        nco.source = 'Original data found at: https://e4ftl01.cr.usgs.gov/MOLT/MOD17A2H.006/'
        nco.references = 'Running, S., Mu, Q., Zhao, M. (2015). MOD17A2H MODIS/Terra Gross Primary Productivity 8-Day L4 Global 500m SIN Grid V006 [Data set]. NASA EOSDIS Land Processes DAAC. Accessed 2020-10-23 from https://doi.org/10.5067/MODIS/MOD17A2H.006'
        nco.history = """[2020-10-23] File created."""
        nco.comment = 'This data has been aggregated from the original MOD17A2H 8-day, 500 m data by Russell Doughty. Please cite the references if you use this data in your project.'
        
        # Get shape, extent, coordinates from one of the input rasters
        ds = gdal.Open(raster_list[0])
    
        a = ds.ReadAsArray()
        nlat, nlon = np.shape(a)
    
        b = ds.GetGeoTransform()  # bbox, interval
        lon_list = np.arange(nlon)*b[1]+b[0]
        lat_list = np.arange(nlat)*b[5]+b[3]
    
        # Create dimensions, variables and attributes
        nco.createDimension('lon', nlon)
        nco.createDimension('lat', nlat)
        nco.createDimension('time', None)
        
        # Time
        time = nco.createVariable('time', 'f8', ('time',))
        time.standard_name = 'time'
        time.calendar = "gregorian"
        time.units = 'days since %s-01-01 00:00:00' % str(years[year])
        dates = [datetime(years[year],1,1)+n*timedelta(days=8) for n in range(0, 46)] # list of dates
        
        # Lon
        lon = nco.createVariable('lon', 'f4', ('lon',))
        lon.standard_name = 'longitude'
        lon.units = 'degrees_east'
        
        # Lat
        lat = nco.createVariable('lat', 'f4', ('lat',))
        lat.standard_name = 'latitude'
        lat.units = 'degrees_north'
    
        # CRS - These values pulled from gdalinfo after converting GeoTiffs to WGS84
        crs = nco.createVariable('crs', 'i4')
        crs.long_name = 'World Geodetic System 1984 (WGS84)'
        crs.grid_mapping_name = 'latitude_longitude'
        crs.longitude_of_prime_meridian = 0.0
        crs.semi_major_axis = 6378137.0
        crs.inverse_flattening = 298.257223563
        crs.spatial_ref = 'GEOGCS["WGS 84",DATUM["World Geodetic System 1984",ELLIPSOID["WGS 84",6378137,298.257223563,LENGTHUNIT["metre",1]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],CS[ellipsoidal,2],AXIS["geodetic latitude (Lat)",north,ORDER[1],ANGLEUNIT["degree",0.0174532925199433]],AXIS["geodetic longitude (Lon)",east,ORDER[2],ANGLEUNIT["degree",0.0174532925199433]],USAGE[SCOPE["unknown"],AREA["World"],BBOX[-90,-180,90,180]],ID["EPSG",4326]]'

        # GPP
        GPP = nco.createVariable('GPP', float,  ('time', 'lat', 'lon'),
                                  zlib=True, chunksizes=None, least_significant_digit=3, fill_value=-9999)
        GPP.long_name = 'Gross Primary Production'
        GPP.units = 'g C m-2 day-2'
        GPP.scale_factor = 0.01
        GPP.add_offset = 0.0
        GPP.grid_mapping = 'crs'
        GPP.set_auto_maskandscale(False)
    
        # Write lon,lat
        lon[:] = lon_list
        lat[:] = lat_list

        # Fill missing DOYs with dummy rasters for 2000 and 2001
        if years[year] == 2000:
            for i in range(len(raster_list_year) + 7):
                time[i] = date2num(dates[i],units=time.units,calendar=time.calendar) # netCDF4 function that translates datetime object into proper format for .nc
                if i > 6:
                    layer = gdal.Open(raster_list_year[i - 7])
                    array = layer.ReadAsArray()  # data
                    GPP[i, :, :] = array
                else:
                    GPP[i, :, :] = dummy
            nco.close()
            print('I have created the file: %s\n' % out_file)

        elif years[year] == 2001:
            for i in range(len(raster_list_year)+1):
                time[i] = date2num(dates[i],units=time.units,calendar=time.calendar) # netCDF4 function that translates datetime object into proper format for .nc
                if i == 22:
                    GPP[i, :, :] = dummy
                elif i < 22:
                    layer = gdal.Open(raster_list_year[i])
                    array = layer.ReadAsArray()  # data
                    GPP[i, :, :] = array
                elif i > 22:
                    layer = gdal.Open(raster_list_year[i - 1])
                    array = layer.ReadAsArray()  # data
                    GPP[i, :, :] = array
            nco.close()
            print('I have created the file: %s\n' % out_file)                    
       
        else:
            # For all years after 2001. Step through each raster for the year, writing time and data to NetCDF
            for i in range(len(raster_list_year)):
                time[i] = date2num(dates[i],units=time.units,calendar=time.calendar) # netCDF4 function that translates datetime object into proper format for .nc
                layer = gdal.Open(raster_list_year[i])
                array = layer.ReadAsArray()  # data
                GPP[i, :, :] = array
            nco.close()
            print('I have created the file: %s\n' % out_file)

tiff_to_netcdf(input_dir, output_dir, years)