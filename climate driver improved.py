# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 17:41:46 2019

@author: Jared

Note: this code performs the same functions as climate driver.py. But this code has each individual
contributor to d_ag split up into seperate functions. It also outputs them in a format which mirrors the inputs
in spatial_temp_diff.py

Since Drew is most interested in Ozone, I'm going to start with ozone, understanding that this same process
can be easily replicated for the other contributors, once all is set up
"""
#%% Load in modules
from netCDF4 import Dataset
import numpy as np
import math
#%% Load in Data

my_cdf_file = "C:/Users/Jared/Projects_for_Shindell/Climate_Models/CESM_final_temp_average.nc"
fh = Dataset(my_cdf_file, mode='r')
#print(fh)
my_ozone_file = "C:/Users/Jared/Projects_for_Shindell/Climate_Models/CESM2_final_ozone_allRuns.nc"
oh = Dataset(my_ozone_file, mode='r')
#print(oh)
O_lats = oh.variables['lat'][:] #this is just a list of latitude degrees binwidth 0.5 degrees
O_lons = oh.variables['lon'][:] #longitude degrees, binwidth 0.5 degrees
#print(O_lats.shape)

M7_3 = oh.variables['M7_3'][:] #M7[#][#][#] <-- M7 ozone value at a certain month, latitude, longitude
M7_4 = oh.variables['M7_4'][:]
M12_3 = oh.variables['M12_3'][:]
M12_4 = oh.variables['M12_4'][:]

#%% Extracting Data

lons = fh.variables['lon'][:]
lats = fh.variables['lat'][:]
temps = fh.variables['temp_diff'][:]
crops = {"wheat": [0.024, 0.138, 137, 2.34, 25], "maize": [0.024, 0.034, 124, 2.83, 20], "rice": [0.032, 0.020, 202, 2.47, 25]} 
#for each value, index 0 represents extratropics, index 1 represents tropics
#index 2 represents alpha value for ozone
#index 3 represents beta value for ozone
#index 4 represents run constant
#soybeans have been withheld
#%% initializing tropics/extratropics. 
#the idea behind these is that one will always be zero. I use them to add temperature change for tropical/extratropical regions
tropics = []
for ycoord in lats:
    if abs(ycoord) <= 30:
        tropics.append(1)
    else:
        tropics.append(0)
extratropics = []
for thing in tropics:
    if thing == 1:
        extratropics.append(0)
    else:
        extratropics.append(1)
#%% initializing precip_map
my_precip_file = "C:/Users/Jared/Projects_for_Shindell/Climate_Models/CESM_final_precip_average.nc"
precip_data = Dataset(my_precip_file, mode = 'r')
precip = precip_data.variables['temp_diff'][:] #yes, I know this says temp diff. It's really precip. Settle down, settle down.

#%% Ozone
def ozoneDriver(crop, year, CO2_conc):
    if crop == "maize":
        run3 = M7_3
        run4 = M7_4
    else: 
        run3 = M12_3
        run4 = M12_4
    runs = [run3, run4]
    three = np.array([])
    four = np.array([])
    runCounter = 0
    for run in runs:
        d_ozone = []
        y_counter = 0
        for ycoord in lats:
            x_counter = 0
            lat_degree_ozone = []
            for xcoord in lons:
                ozone_conc = float(run[y_counter][x_counter])
                alpha = float(crops[crop][2])
                beta = float(crops[crop][3])
                run_constant = float(crops[crop][4])
               
                numerator = 1 / e**((ozone_conc / alpha)**beta)
                denominator = 1 / e**((run_constant / alpha)**beta)
                ozone_change = (1 - (numerator / denominator)) * 100
                lat_degree_ozone.append(ozone_change.real)#still don't know why this is complex in the first place. 
                #the variability between datapoints seems right. But the 
                x_counter += 1
            d_ozone.append(lat_degree_ozone)
            y_counter += 1
        if(runCounter == 0):
            three = np.asarray(d_ozone)
            runCounter += 1
        if(runCounter == 1):
            four = np.asarray(d_ozone)
    final = three - four
    CDFMaker(final, 'CESM_ozone_yield_loss_change.nc')
    return final

#%% CO2
def CO2Driver(crop, year, CO2_conc):
    dCO2 = []
    y_counter = 0
    for ycoord in lats:
        x_counter = 0
        lat_degree_CO2 = []
        for xcoord in lons:
            lat_degree_CO2.append(CO2_conc*0.0006) # 0.06% per ppm (from Challinor et al., NCC, 2014) #0.7 ppm is CO2 value
            x_counter += 1
        dCO2.append(lat_degree_CO2)
        y_counter += 1
    myArray = np.asarray(dCO2)
    CDFMaker(myArray, 'CESM_CO2_yield_loss.nc')
    return myArray

#%% Temperature
def tempDriver(crop, year, CO2_conc):
    d_temp = []
    y_counter = 0
    for ycoord in lats:
        x_counter = 0
        lat_degree_temp = []
        for xcoord in lons:
            crop_temp = (crops[crop][1]*tropics[y_counter] + crops[crop][0]*extratropics[y_counter])*temps[year][y_counter][x_counter]
            lat_degree_temp.append(crop_temp)
            x_counter += 1
        d_temp.append(lat_degree_temp)
        y_counter += 1
    return np.asarray(d_temp)

#%% Precipitation
def precipDriver(crop, year, CO2_conc):
    d_precip = []
    y_counter = 0
    for ycoord in lats:
        x_counter = 0
        lat_degree_precip = [] #included for sake of graphing
        for xcoord in lons:
            precip_change = (precip[year][y_counter][x_counter] - precip[0][y_counter][x_counter])/ precip[0][y_counter][x_counter]
            #why am I doing precip change? I suppose that makes sense, but I still want to double check that that's the right thing. 
            lat_degree_precip.append(precip_change*0.0053) #this gives % yield loss by doing %precip change*%yield loss/%precip change
            x_counter += 1
        d_precip.append(lat_degree_precip) 
        y_counter += 1
        
    return np.asarray(d_precip)
#this output is readable by the spatial temp diff and spatial precip diff files. 
#Now, I just don't know what programs you want me to run, or what data you have to compare mine against.
#%% Total Yield Loss
def totalYieldLoss(crop, year, CO2_conc):
    d_ag = []
    y_counter = 0
    if crop == "maize":
        run = M7
    else: 
        run = M12
    for ycoord in lats:
        x_counter = 0
        lat_degree_yield = []
        for xcoord in lons:
            CO2 = CO2_conc*0.0006 # 0.06% per ppm (from Challinor et al., NCC, 2014)
            crop_temp = (crops[crop][1]*tropics[y_counter] + crops[crop][0]*extratropics[y_counter])*temps[year][y_counter][x_counter]
            precip_change = (precip[year][y_counter][x_counter] - precip[0][y_counter][x_counter])/ precip[0][y_counter][x_counter]
            ozone_conc = run[1][y_counter][x_counter]
            ozone_change = 1 - math.exp(-(ozone_conc / crops[crop][2])**crops[crop][3]) / math.exp(-(crops[crop][4] / crops[crop][2])**crops[crop][3]) 
            lat_degree_yield.append(CO2 + crop_temp + precip_change + ozone_change)
            x_counter += 1
        d_ag.append(lat_degree_yield)
        y_counter += 1
    return np.asarray(d_ag)
#%% 
def yieldLoss4Graph(crop, year, CO2_conc):
    d_ag = []
    y_counter = 0
    for ycoord in lats:
        x_counter = 0
        lat_degree_loss = [] #included for sake of graphing
        for xcoord in lons:
            CO2 = CO2_conc*0.0006 # 0.06% per ppm (from Challinor et al., NCC, 2014) //should be crops[crop][0] for extratropics
            crop_temp = (crops[crop][1]*tropics[y_counter] + crops[crop][0]*extratropics[y_counter])*temps[y_counter][x_counter]
            precip_change = (precip[y_counter][x_counter] - precip[y_counter][x_counter])/ precip[y_counter][x_counter]
            x_counter += 1
            lat_degree_loss.append(CO2 + crop_temp + precip_change)
        y_counter += 1
        d_ag.append(lat_degree_loss)
    myArray = (np.asarray(d_ag) * 100) + ozoneDriver(crop, year, CO2_conc)
    CDFMaker(myArray, 'CESM_Final_yield_loss_combined.nc')
        
    return myArray

# I'm worried that this isn't outputting anything correct aside from Ozone. It seems like the temperature values are either two low or not properly reflected
#%% export netCDF file
def CDFMaker(myArray, myCDF):
    latitude   = np.arange(-90,90,180/len(myArray))
    longitude  = np.arange(-180,180,180/len(myArray))
    f1 = Dataset(myCDF,'w',format='NETCDF4_CLASSIC')
    
    lat                             = f1.createDimension('lat',len(latitude))
    lon                             = f1.createDimension('lon',len(longitude))
    
    lat                             = f1.createVariable('lat',np.double,('lat',))
    lon                             = f1.createVariable('lon',np.double,('lon',))
    yield_loss                      = f1.createVariable('yield_loss',np.float32,('lat','lon')) 
    
    f1.description                  = 'In file: total agricultural yield loss due to methane emissions'
    
    lat.units                       = 'degrees_north'
    lat.long_name                   = 'latitude'
    lat.comment                     = 'center of grid cell'
    lon.units                       = 'degrees_east'
    lon.long_name                   = 'longitude'
    lon.comment                     = 'center of grid cell'
    yield_loss.units                = '% loss'
    yield_loss.long_name            = 'yield_loss'
    
    f1.variables['lat'][:]          = latitude[:]
    f1.variables['lon'][:]          = longitude[:]
    f1.variables['yield_loss'][:,:]  = myArray[:,:]
    
    f1.close()
    

#%% Main
if __name__ == "__main__":
    print(yieldLoss4Graph("maize", 20,  1))
    #print(ozoneDriver("maize", 20,  1))
    #myCDF = 'CESM_final_temp_average.nc'
    #print(precipDriver("wheat", 20, CO2_conc = float(input("Enter Change in CO2 Concentration (PPM)")))[100])

    
        
