#-------------------------------------------------------------------------------
# Name:        Timeseries TIFF to netCDF
#
# Purpose:     Creates annual netCDF4 files containing timeseries data for a
#              single year or range of years.
#
# Author:      Russell Doughty, PhD
#
# Created:     April 10th, 2020
#-------------------------------------------------------------------------------


import os
from netCDF4 import Dataset,date2num
from datetime import datetime, timedelta
from osgeo import gdal
import numpy as np

# Input arguments for 8-day data
input_dir = '/data/ifs/VPM/product/v20/merge/CMG/8-day_new' # input folder contains TIFFs from all years
output_dir = '/data/ifs/VPM/product/v20/merge/netcdf/CMG_0.05_8-day'
tempo='8-day' # Valid strings are '8-day' or 'monthly'

## Input arguments for monthly data
#input_dir = '/data/ifs/VPM/product/v20/merge/CMG/monthly_new' # input folder contains TIFFs from all years
#output_dir = '/data/ifs/VPM/product/v20/merge/netcdf/CMG_0.05_monthly'
#tempo='monthly' # Valid strings are '8-day' or 'monthly'

years=range(2000,2019+1) # Start and end year


def tiff_to_netcdf(input_dir, output_dir, years, tempo):
    
    # Only monthly and 8-day accepted
    if tempo != 'monthly':
        if tempo != '8-day':
            quit('ERROR: temp variable should be either monthly or 8-day')      
    
    raster_list = sorted([os.path.join(input_dir, l) for l in os.listdir(input_dir) if l.endswith('.tif')]) # All rasters in input_dir

    for year in range(len(years)):
        
        raster_list_year = [r for r in raster_list if '.'+str(years[year]) in r] # List of rasters for the year

        # Ensure number of input rasters is correct
        if tempo == 'monthly' and len(raster_list_year) != 12:
            quit('ERROR: There should be 12 input rasters for the year '+str(years[year]+', but there are only '+ len(raster_list_year)+'.'))
        if tempo == '8-day' and len(raster_list_year) != 46:
            quit('ERROR: There should be 46 input rasters for the year '+str(years[year]+', but there are only '+ len(raster_list_year)+'.'))
 
        # Printing for visual quality control
        print('Packing up these '+str(len(raster_list_year))+' rasters into the .nc file for the year '+str(years[year])+':')
        print("\n".join(raster_list_year))

        # Create NetCDF file
        if tempo == 'monthly':
            out_file = os.path.join(output_dir,'GPP.VPM.v20.'+str(years[year])+'.monthly.CMG_0.05.nc')
        else:
            out_file = os.path.join(output_dir,'GPP.VPM.v20.'+str(years[year])+'.8-day.CMG_0.05.nc')
        
        nco = Dataset(out_file, 'w', clobber=True, format="NETCDF4")
        
        # Meta-data
        nco.Conventions = 'CF-1.7'
        nco.title = 'Vegetation Photosynthesis Model (VPM)'
        nco.institution = 'Department of Microbiology and Plant Biology, University of Oklahoma, Norman, Oklahoma'
        nco.source = 'http://eomf.ou.edu'
        nco.references = 'Zhang, Y., Xiao, X., Wu, X., Zhou, S., Zhang, G., Qin, Y. and Dong, J., 2017. A global moderate resolution dataset of gross primary production of vegetation for 2000-2016. Scientific data, 4, p.170165.'
        nco.history = """
        [2020-04-10] Added units to lat and lon for CF compliance. Also set least significant digit for GPP to 3 decimal places to reduce file size.
        [2020-04-08] Scale factor has been corrected to 1.
        [2020-04-07] CRS now set to WGS84. netCDF files now provided on annual basis.
        [2020-04-04] First creation of VPM in netCDF format.
        """
        nco.comment = 'This product has been aggregated from the original 8-day, 500 m data. Please cite the references if you use this data in your project. Contact Dr. Xiangming Xiao at xiangming.xiao@ou.edu with your questions.'
        
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

        if tempo == 'monthly':
            time.units = 'days since %s-01-01 00:00:00' % str(years[year])
            dates = [datetime(years[year],1+n,1) for n in range(len(raster_list_year))] # list of dates
        elif tempo == '8-day':
            time.units = 'days since %s-01-01 00:00:00' % str(years[year])
            dates = [datetime(years[year],1,1)+n*timedelta(days=8) for n in range(len(raster_list_year))] # list of dates
        
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
        if tempo == 'monthly':
            GPP.units = 'g C m-2 month-2'
        elif tempo == '8-day':
            GPP.units = 'g C m-2 day-2'
        GPP.scale_factor = 1.0
        GPP.add_offset = 0.0
        GPP.grid_mapping = 'crs'
        GPP.set_auto_maskandscale(False)
    
        # Write lon,lat
        lon[:] = lon_list
        lat[:] = lat_list
       
        # Step through each raster for the year, writing time and data to NetCDF
        for i in range(len(raster_list_year)):
    
            time[i] = date2num(dates[i],units=time.units,calendar=time.calendar) # netCDF4 function that translates datetime object into proper format for .nc
            layer = gdal.Open(raster_list_year[i])
            array = layer.ReadAsArray()  # data
            GPP[i, :, :] = array
    
        nco.close()
        
        print('I have created the file: %s\n' % out_file)
    
tiff_to_netcdf(input_dir,output_dir,years,tempo)