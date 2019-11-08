# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 14:53:01 2019

@author: Jared
"""
#%% importing modules
import numpy as np
import math
import matplotlib.pyplot as plt
#os.environ['PROJ_LIB'] = 'C:/Users/Jared/Anaconda3/share/proj'
#from mpl_toolkits.basemap import Basemap #I can't fucking use basemap on my computer. I don't know why the hell that is.
#%% reading in data
d_ag = yieldLoss4Graph("rice", 20, 1)

#%% Create data
lon = [var[0] for var in d_ag]
lat = [var[1] for var in d_ag]
data = [var[2] for var in d_ag]
print(lat)
plt.scatter(lat, lon, c=data)
plt.show()

#%%

# Create heatmap
heatmap, xedges, yedges = np.histogram2d(x, y, bins=(26,26)) #there's something going on here that sorts it
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]] #also, this is number of occurances. So Idk if this will work
print(xedges)
#print(xedges[0], xedges[-1], yedges[0], yedges[-1])

# Plot heatmap
plt.clf()
plt.title('Pythonspot.com heatmap example')
plt.ylabel('y')
plt.xlabel('x')
plt.imshow(heatmap, extent=extent)
plt.show()
