from netCDF4 import Dataset
import numpy as np

data = precipDriver("rice", 20)
latitude   = np.arange(-90,90,180/len(data))
longitude  = np.arange(-180,180,360/len(data)) 
print(len(latitude))
fake_data  = np.zeros((len(latitude),len(longitude)))  # generate a blank numpy array

####################################################################################################
#### SAVE ALL RESULTS TO A NETCDF

f1 = Dataset('./example.nc','w',format='NETCDF4_CLASSIC')

lat                             = f1.createDimension('lat',len(latitude))
lon                             = f1.createDimension('lon',len(longitude))

lat                             = f1.createVariable('lat',np.double,('lat',))
lon                             = f1.createVariable('lon',np.double,('lon',))
some_data                       = f1.createVariable('some_data',np.float32,('lat','lon'))

# Global Attributes
f1.description                  = 'This are global attributes for your file.'

# Variable Attributes
lat.units                       = 'degrees_north'
lat.long_name                   = 'latitude'
lat.comment                     = 'center of grid cell'
lon.units                       = 'degrees_east'
lon.long_name                   = 'longitude'
lon.comment                     = 'center of grid cell'
some_data.units                 = 'it is good to describe your data'
some_data.long_name             = 'replace long_name with whatever. that is another option'

f1.variables['lat'][:]          = latitude[:]
f1.variables['lon'][:]          = longitude[:]
f1.variables['some_data'][:,:]  = fake_data[:,:]

f1.close()
