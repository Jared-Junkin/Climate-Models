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

my_cdf_file = "C:/Users/Jared/Documents/4jared/CESM2_spatial_temp_diff.nc"
fh = Dataset(my_cdf_file, mode='r')
#print(fh)
my_ozone_file = "C:/Users/Jared/Documents/4jared/CESM2_agriozone.nc"
oh = Dataset(my_ozone_file, mode='r')
#print(oh)
O_lats = oh.variables['latitude'][:] #this is just a list of latitude degrees binwidth 0.5 degrees
O_lons = oh.variables['longitude'][:] #longitude degrees, binwidth 0.5 degrees
months = oh.variables['month'][:]
M7 = oh.variables['M7_Run3'][:] #M7[#][#][#] <-- M7 ozone value at a certain month, latitude, longitude
M12 = oh.variables['M12_Run3'][:]
#print(M7[1][1][1])
#print(M12[1][1][1])
#print(type(months[1])) #floats


#%% Extracting Data

lons = fh.variables['longitude'][:]
print(lons)
lats = fh.variables['latitude'][:]
temps = fh.variables['CESM2_ann_avg_temp1'][:]
crops = {"wheat": [0.024, 0.138, 137, 2.34, 25], "maize": [0.024, 0.034, 124, 2.83, 20], "rice": [0.032, 0.020, 202, 2.47, 25]} 

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
my_precip_file = "C:/Users/Jared/Documents/4jared/CESM2_spatial_precip_diff.nc"
precip_data = Dataset(my_precip_file, mode = 'r')
precip = precip_data.variables['CESM2_ann_avg_precip1'][:]
#print(precip[1][100])
#%% Ozone
def ozoneDriver(crop, year, CO2_conc):
    d_ozone = []
    y_counter = 0
    if crop == "maize":
        run = M7
    else: 
        run = M12       
    for ycoord in lats:
        x_counter = 0
        lat_degree_ozone = []
        for xcoord in lons:
            ozone_conc = run[1][y_counter][x_counter]
            ozone_change = 1 - math.exp(-(ozone_conc / crops[crop][2])**crops[crop][3]) / math.exp(-(crops[crop][4] / crops[crop][2])**crops[crop][3]) 
            lat_degree_ozone.append(ozone_change)
            x_counter += 1
        d_ozone.append(lat_degree_ozone)
        y_counter += 1
    return np.asarray(d_ozone)

#%% CO2
def CO2Driver(crop, year, CO2_conc):
    dCO2 = []
    y_counter = 0
    for ycoord in lats:
        x_counter = 0
        lat_degree_CO2 = []
        for xcoord in lons:
            lat_degree_CO2.append(CO2_conc*0.0006) # 0.06% per ppm (from Challinor et al., NCC, 2014)
            x_counter += 1
        dCO2.append(lat_degree_CO2)
        y_counter += 1
    return np.asarray(dCO2)

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
    if crop == "maize":
        run = M7
    else: 
        run = M12
    for ycoord in lats:
        x_counter = 0
        lat_degree_loss = [] #included for sake of graphing
        for xcoord in lons:
            CO2 = CO2_conc*0.0006 # 0.06% per ppm (from Challinor et al., NCC, 2014) //should be crops[crop][0] for extratropics
            crop_temp = (crops[crop][1]*tropics[y_counter] + crops[crop][0]*extratropics[y_counter])*temps[year][y_counter][x_counter]
            precip_change = (precip[year][y_counter][x_counter] - precip[0][y_counter][x_counter])/ precip[0][y_counter][x_counter]
            ozone_conc = run[1][y_counter][x_counter]
            ozone_change = 1 - math.exp(-(ozone_conc / crops[crop][2])**crops[crop][3]) / math.exp(-(crops[crop][4] / crops[crop][2])**crops[crop][3]) 
            #d_ag.append([ycoord, xcoord, CO2 + crop_temp + precip_change + ozone_change])
            x_counter += 1
            lat_degree_loss.append(CO2 + crop_temp + precip_change + ozone_change)
        y_counter += 1
        d_ag.append(lat_degree_loss)
        
    return np.asarray(d_ag)

#%% Main
if __name__ == "__main__":
    #print(yieldLoss4Graph("rice", 20, CO2_conc = float(input("Enter Change in CO2 Concentration (PPM)")))[1])
    print(precipDriver("wheat", 20, CO2_conc = float(input("Enter Change in CO2 Concentration (PPM)")))[100])

    
        
