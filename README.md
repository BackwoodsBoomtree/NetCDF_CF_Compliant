# NetCDF CF Compliant Files

The codes in this repo can be used to make CF compliant NetCDF files from directories of tif files.

## Notes

* The portions of the code that account for missing files in 2000 and 2001 likely need to be updated when running for those years.
* The code for packaging the VPM and MOD GPP files likely needs to be updated when they are needed. Use MCD43C4.py for reference when updating.
* Pay attention to the size of the resultant nc files. For instance, for 0.05-degree data, package 8-day and monthly tif files by year, and package daily files by month.

## Resources
* https://cfconventions.org/
* https://unidata.github.io/netcdf4-python/

## Author

Russell Doughty, PhD
russell.doughty@ou.edu
