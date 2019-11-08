from netCDF4 import Dataset
import numpy as np
import math

#my_cdf_file = "C:/Users/Jared/Projects for Shindell/zzz.nc" #original
my_cdf_file = "C:/Users/Jared/Documents/4jared/CESM2_spatial_temp_diff.nc"
fh = Dataset(my_cdf_file, mode='r')
#print(fh)
my_ozone_file = "C:/Users/Jared/Documents/4jared/CESM2_agriozone.nc"
oh = Dataset(my_ozone_file, mode='r')
print(oh)
O_lats = oh.variables['latitude'][:] #this is just a list of latitude degrees binwidth 0.5 degrees
O_lons = oh.variables['longitude'][:] #longitude degrees, binwidth 0.5 degrees
months = oh.variables['month'][:]
M7 = oh.variables['M7_Run3'][:] #M7[#][#][#] <-- M7 ozone value at a certain month, latitude, longitude
M12 = oh.variables['M12_Run3'][:]
print(M7[1][1][1])
print(M12[1][1][1])
#print(type(months[1])) #floats
#%% Extracting Data 
#lons = fh.variables['lon'][:]
#lats = fh.variables['lat'][:]
#temps = fh.variables['TEMPANOMALY'][:].tolist()
lons = fh.variables['longitude'][:]
lats = fh.variables['latitude'][:]
temps = fh.variables['CESM2_ann_avg_temp1'][:]
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
my_precip_file = "C:/Users/Jared/Documents/4jared/CESM2_spatial_precip_diff.nc"
precip_data = Dataset(my_precip_file, mode = 'r')
precip = precip_data.variables['CESM2_ann_avg_precip1'][:]
print(precip[1][100])
#%% main function (calculates losses due to temperature and CO2 changes)
def d_ag(crop, year, CO2_conc = float(input("Enter Change in CO2 Concentration (PPM)"))):
    dCO2 = []
    d_temp = []
    d_precip = [] 
    d_ozone = []
    y_counter = 0
    if crop == "maize":
        run = M7
    else: 
        run = M12
    for ycoord in lats:
        x_counter = 0
        lat_degree_precip = []
        for xcoord in lons:
            dCO2.append([ycoord, xcoord, CO2_conc*0.0006]) # 0.06% per ppm (from Challinor et al., NCC, 2014)
            crop_temp = (crops[crop][1]*tropics[y_counter] + crops[crop][0]*extratropics[y_counter])*temps[year][y_counter][x_counter]
            d_temp.append([ycoord, xcoord, crop_temp])
            precip_change = (precip[year][y_counter][x_counter] - precip[0][y_counter][x_counter])/ precip[0][y_counter][x_counter]
            lat_degree_precip.append([ycoord, xcoord, precip_change*0.0053]) #this gives % yield loss by doing %precip change*%yield loss/%precip change
            ozone_conc = run[1][y_counter][x_counter]
            ozone_change = 1 - math.exp(-(ozone_conc / crops[crop][2])**crops[crop][3]) / math.exp(-(crops[crop][4] / crops[crop][2])**crops[crop][3]) 
            d_ozone.append(ozone_change)
            #note, this doesn't fix the issue of binwidth. It just ignores it.             
          
            #Look at table one in paper to find yield loss given Ozone
            #M7 and M12 are ozone exposure given seven hours with highest level and 12 with highest level.
            ### these will come in a dataset and be given to you by greg
            #A and b are empirical constants. Just put them in. 
            x_counter += 1
        d_precip.append(lat_degree_precip)
        y_counter += 1
    return d_ozone
if __name__ == "__main__":
    print(d_ag("rice", 20)[1])
    
    
    
    
####### Extra Info #######
#%% Dimensions of Various data fromes and names of data
'''
zzz  lat(90), lon(180), temp(45) **no time
lon 
lat
TEMPANOMALY

CESM2 spatial temp --> lat(192), lon(288), years(45)
CESM spatial precip --> same
longitude
latitude
CESM2_ann_avg_temp1

CESM2_agriozone.nc
lat 360 lon 720 12 months
M7_Run3(month,latitude,longitude), 
float32 M12_Run3(month,latitude,longitude), 
float32 M7_Run4(month,latitude,longitude), 
float32 M12_Run4(month,latitude,longitude)
########################################### Note that these variables exist in agriozone and I haven't accessed them yet
GFDL spatial temp --> lat(180), lon(288), years(40)
GFDL spatial precip --> same

GISS spatial precip --> lat(90), lon(144), years(45)
GISS spatial temp --> same

Original file with tempanomaly --> lat(90), lon(180)

'''
#%% fixing temp to remove Nones. Note: this causes code to break when from other source than original nc
'''
for y in range(len(temps)):
     for x in range(len(temps[y])):
         if type(temps[y][x]) != float:
            temps[y][x] = 0.0 
'''