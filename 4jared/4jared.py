from netCDF4 import Dataset
import numpy as np

my_cdf_file = './CESM2_spatial_temp_diff.nc'
fh = Dataset(my_cdf_file, mode='r')

lons = fh.variables['longitude'][:]
lats = fh.variables['latitude'][:]
temps = fh.variables['CESM2_ann_avg_temp1'][:]
crops = {"wheat": [0.024, 0.138], "maize": [0.024, 0.034], "rice": [0.032, 0.020]} #for each value, index 0 represents extratropics, index 1 represents tropics

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

GFDL spatial temp --> lat(180), lon(288), years(40)
GFDL spatial precip --> same

GISS spatial precip --> lat(90), lon(144), years(45)
GISS spatial temp --> same

Original file with tempanomaly --> lat(90), lon(180)

'''

#print(precip) ###what is the difference between ann_avg_precip1 and 2?
#print(len(temps), len(temps[1]), len(temps[1][1]))
#%% fixing temp to remove Nones
for y in range(len(temps)):
     for x in range(len(temps[y])):
         if type(temps[y][x]) != float:
            temps[y][x] = 0.0 #there were a bunch of "None"'s in temps. I replaced them with 0's. What should I do long-term        
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
#%% initializing precip_map (this will later be read in)
my_precip_file = "C:/Users/Jared/Documents/4jared/CESM2_spatial_precip_diff.nc"
precip_data = Dataset(my_precip_file, mode = 'r')
precip = precip_data.variables['CESM2_ann_avg_precip1'][:].tolist()
print(len(precip), len(precip[1]), len(precip[1][1]), type(precip[1][1][1]))
print(precip_data)
#%% main function (calculates losses due to temperature and CO2 changes)
def d_ag(crop, CO2_conc = float(input("Enter Change in CO2 Concentration (PPM)"))):
    dCO2 = []
    d_temp = []
    d_precip = [] 
    y_counter = 0
    for ycoord in lats:
        x_counter = 0
        for xcoord in lons:
            dCO2.append([ycoord, xcoord, CO2_conc*0.006]) # is it 0.006 or 0.0006?
            crop_temp = (crops[crop][1]*tropics[y_counter] + crops[crop][0]*extratropics[y_counter])*temps[y_counter][x_counter]
            #print(tropics[y_counter])
            d_temp.append([ycoord, xcoord, crop_temp])
            #d_precip.append([ycoord, xcoord])
            x_counter+= 1
        y_counter += 1
    return dCO2
if __name__ == "__main__":
    print(d_ag("rice"))