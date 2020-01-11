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
    temps1 = fh.variables['GFDL_ann_avg_precip1'][:]
    temps2 = fh.variables['GFDL_ann_avg_precip2'][:]
    diff = temps1[10:44,:,:] - temps2[10:44,:,:]
    mean_diff = np.average(diff[:,:],axis=0)
    return CDF_Outputter(mean_diff)

#%% calculate average for Ozone
def OzoneAverage(my_cdf_ozone):
    oh = Dataset(my_cdf_ozone, mode='r')
    M7_3 = oh.variables['M7_Run3'][:]
    M7_4 = oh.variables['M7_Run4'][:]
    #print(type(M7_3))
    M12_3 = oh.variables['M12_Run3'][:]
    M12_4 = oh.variables['M12_Run4'][:]
    
    
    M7_3_mean = M7_3[4:8,:,:]
    M7_4_mean = M7_4[4:8,:,:]
    
    M7_3_avg = np.average(M7_3_mean,axis=0)
    M7_4_avg = np.average(M7_4_mean,axis=0)
    
    M12_3_mean = M12_3[4:8,:,:]
    M12_4_mean = M12_4[4:8,:,:]
    
    M12_3_avg = np.average(M12_3_mean,axis=0)
    M12_4_avg = np.average(M12_4_mean,axis=0)
    
    
    return OzoneCDF(M7_3_avg, M7_4_avg, M12_3_avg, M12_4_avg)

#generate for run three, run for. Process. Then subtract to get difference. 

#%% Ozone CDF mapper (only difference is that there's an extra dimension in Ozone)
def OzoneCDF(M7_3_mean, M7_4_mean, M12_3_mean, M12_4_mean):
    #this assumes that all things in the array will be the same length
        
    latitude   = np.arange(-90,90,180/len(M7_3_mean))
    longitude  = np.arange(-180,180,180/len(M7_3_mean))
    f1 = Dataset('GFDL_final_ozone_allRuns.nc','w',format='NETCDF4_CLASSIC')
    
    lat                             = f1.createDimension('lat',len(latitude))
    lon                             = f1.createDimension('lon',len(longitude))
    
    lat                             = f1.createVariable('lat',np.double,('lat',))
    lon                             = f1.createVariable('lon',np.double,('lon',))
    
    M7_3                              = f1.createVariable('M7_3', np.float32,('lat','lon'))
    M7_4                              = f1.createVariable('M7_4', np.float32,('lat','lon'))

    M12_3                             = f1.createVariable('M12_3', np.float32,('lat','lon'))
    M12_4                             = f1.createVariable('M12_4', np.float32,('lat','lon'))
    
    f1.description                  = 'In file: mean ozone exposure in PPB for growing season months.'
    
    lat.units                       = 'degrees_north'
    lat.long_name                   = 'latitude'
    lat.comment                     = 'center of grid cell'
    lon.units                       = 'degrees_east'
    lon.long_name                   = 'longitude'
    lon.comment                     = 'center of grid cell'
    M7_3.units                      = 'PPB'
    M7_4.units                      = 'PPB'
    M7_3.long_name                  = 'M7_run3'
    M7_4.long_name                  = 'M7_run4'
    M12_3.units                     = 'PPB'
    M12_4.units                     = 'PPB'
    M12_3.long_name                 = 'M12_run3'
    M12_4.long_name                 = 'M12run_4'
    
    f1.variables['lat'][:]          = latitude[:]
    f1.variables['lon'][:]          = longitude[:]
    f1.variables['M7_3'][:,:]       = M7_3_mean[:,:]
    f1.variables['M7_4'][:,:]       = M7_4_mean[:,:]
    f1.variables['M12_3'][:,:]      = M12_3_mean[:,:]
    f1.variables['M12_4'][:,:]      = M12_4_mean[:,:]

    
    f1.close()
#%% output net cdf file
def CDF_Outputter(mean_diff): #takes in numpy array, outputs netCDF file with lat and lon variables
    latitude   = np.arange(-90,90,180/len(mean_diff))
    longitude  = np.arange(-180,180,180/len(mean_diff))
    f1 = Dataset('GFDL_final_precip_average_new.nc','w',format='NETCDF4_CLASSIC')
    lat                             = f1.createDimension('lat',len(latitude))
    lon                             = f1.createDimension('lon',len(longitude))
    
    lat                             = f1.createVariable('lat',np.double,('lat',))
    lon                             = f1.createVariable('lon',np.double,('lon',))
    precip_diff                       = f1.createVariable('precip_diff',np.float32,('lat','lon'))
    
    f1.description                  = 'In file: average precip difference between control group and methane simulation.'
    
    lat.units                       = 'degrees_north'
    lat.long_name                   = 'latitude'
    lat.comment                     = 'center of grid cell'
    lon.units                       = 'degrees_east'
    lon.long_name                   = 'longitude'
    lon.comment                     = 'center of grid cell'
    precip_diff.units                 = 'degrees K'
    precip_diff.long_name             = 'mean_diff'
    
    f1.variables['lat'][:]          = latitude[:]
    f1.variables['lon'][:]          = longitude[:]
    f1.variables['precip_diff'][:,:]  = mean_diff[:,:]
    
    f1.close()
    
#how come I'm not using that extra variabel
#%% Main method
if __name__ == "__main__":
    GFDL_temp = "C:/Users/Jared/Projects_for_Shindell/Climate_Models/4jared/GFDL_spatial_temp_diff_05x05.nc"
    #Average(GFDL_temp)
    GFDL_precip = "C:/Users/Jared/Projects_for_Shindell/Climate_Models/4jared/GFDL_spatial_precip_diff_05x05.nc"
    #Average(GFDL_precip)
    GFDL_ozone = "C:/Users/Jared/Projects_for_Shindell/Climate_Models/4jared/GFDL_agriozone.nc"
    OzoneAverage(GFDL_ozone)
    
