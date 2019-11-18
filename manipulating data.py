# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 13:16:34 2019

@author: Jared
"""
#%% Importing Modules
import os
os.environ['PROJ_LIB'] = r'C:\Users\Jared\Anaconda3\pkgs\proj4-5.2.0-ha925a31_1\Library\share'
import matplotlib
from pylab import *
from mpl_toolkits.basemap import Basemap #Basemap isn't working
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

'''
I need to take in a NetCDF file, 
subtract 1 from 2 (2 is the control group)
average the last twenty years of precipitation and temperature,
then output a new netCDF file with the same dimensions

Also need to take in a netcdf file for ozone,
average months 4,5,6,7 of the file, 
subtract 3 from 4 (4 is the control group)
aerage those months
the output a new netCDF file 

'''
#%% calculate average for precipitation and temp.
def Average(my_cdf_file): #takes location of CDF file on computer (string), outputs CDF file
    fh = Dataset(my_cdf_file, mode='r')
    temps1 = fh.variables['GFDL_ann_avg_temp1'][:]
    temps2 = fh.variables['GFDL_ann_avg_temp2'][:]
    diff = temps1[24:44,:,:] - temps2[24:44,:,:]
    mean_diff = np.average(diff[:,:],axis=0)
    return CDF_Outputter(mean_diff)

#%% calculate average for Ozone
def OzoneAverage(my_cdf_ozone):
    oh = Dataset(my_cdf_ozone, mode='r')
    M7_3 = oh.variables['M7_Run3'][:]
    M7_4 = oh.variables['M7_Run4'][:]
    M12_3 = oh.variables['M12_Run3'][:]
    M12_4 = oh.variables['M12_Run4'][:]
    
    M7_diff = M7_3[4:8,:,:] - M7_4[4:8,:,:]
    M12_diff = M12_3[4:8,:,:] - M12_4[4:8,:,:]
    
    M7_mean = np.average(M7_diff[:,:],axis=0)
    M12_mean = np.average(M12_diff[:,:],axis=0)
    return OzoneCDF(M7_mean, M12_mean)
    
#%% Ozone CDF mapper (only difference is that there's an extra dimension in Ozone)
def OzoneCDF(M7_mean, M12_mean):
    if(len(M7_mean) != len(M12_mean)):
        print("Error! M7 and M12 files different lengths!!!")
        
    latitude   = np.arange(-90,90,180/len(M7_mean))
    longitude  = np.arange(-180,180,180/len(M7_mean))
    f1 = Dataset('CESM2_final_ozone_average.nc','w',format='NETCDF4_CLASSIC')
    
    lat                             = f1.createDimension('lat',len(latitude))
    lon                             = f1.createDimension('lon',len(longitude))
    
    lat                             = f1.createVariable('lat',np.double,('lat',))
    lon                             = f1.createVariable('lon',np.double,('lon',))
    
    M7                              = f1.createVariable('M7', np.float32,('lat','lon'))
    M12                             = f1.createVariable('M12', np.float32,('lat','lon'))
    
    f1.description                  = 'In file: mean ozone exposure in PPB for growing season months.'
    
    lat.units                       = 'degrees_north'
    lat.long_name                   = 'latitude'
    lat.comment                     = 'center of grid cell'
    lon.units                       = 'degrees_east'
    lon.long_name                   = 'longitude'
    lon.comment                     = 'center of grid cell'
    M7.units                        = 'PPB'
    M7.long_name                    = 'M7'
    M12.units                       = 'PPB'
    M12.long_name                   = 'M12'
    
    f1.variables['lat'][:]          = latitude[:]
    f1.variables['lon'][:]          = longitude[:]
    f1.variables['M7'][:,:]         = M7_mean[:,:]
    f1.variables['M12'][:,:]        = M12_mean[:,:]
    
    f1.close()
#%% output net cdf file
def CDF_Outputter(mean_diff): #takes in numpy array, outputs netCDF file with lat and lon variables
    latitude   = np.arange(-90,90,180/len(mean_diff))
    longitude  = np.arange(-180,180,180/len(mean_diff))
    f1 = Dataset('GFDL_final_temp_average.nc','w',format='NETCDF4_CLASSIC')
    lat                             = f1.createDimension('lat',len(latitude))
    lon                             = f1.createDimension('lon',len(longitude))
    
    lat                             = f1.createVariable('lat',np.double,('lat',))
    lon                             = f1.createVariable('lon',np.double,('lon',))
    temp_diff                       = f1.createVariable('temp_diff',np.float32,('lat','lon'))
    
    f1.description                  = 'In file: average temp difference between control group and methane simulation.'
    
    lat.units                       = 'degrees_north'
    lat.long_name                   = 'latitude'
    lat.comment                     = 'center of grid cell'
    lon.units                       = 'degrees_east'
    lon.long_name                   = 'longitude'
    lon.comment                     = 'center of grid cell'
    mean_diff.units                 = 'degrees K'
    mean_diff.long_name             = 'mean_diff'
    
    f1.variables['lat'][:]          = latitude[:]
    f1.variables['lon'][:]          = longitude[:]
    f1.variables['temp_diff'][:,:]  = mean_diff[:,:]
    
    f1.close()
    

#%% Main method
if __name__ == "__main__":
    CESM_temp = "C:/Users/Jared/Projects_for_Shindell/Climate_Models/4jared/GFDL_spatial_temp_diff_05x05.nc"
    #Average(CESM_temp)
    CESM_precip = "C:/Users/Jared/Projects_for_Shindell/Climate_Models/4jared/CESM2_spatial_precip_diff_05x05.nc"
    #Average(CESM_precip)
    CESM_ozone = "C:/Users/Jared/Projects_for_Shindell/Climate_Models/4jared/CESM2_agriozone.nc"
    OzoneAverage(CESM_ozone)
    
