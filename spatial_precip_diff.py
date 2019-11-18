# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 11:36:53 2019

@author: Jared
"""

import os
os.environ['PROJ_LIB'] = r'C:\Users\Jared\Anaconda3\pkgs\proj4-5.2.0-ha925a31_1\Library\share'
from pylab import *
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

rc('font',family='serif')

####################################################################################################
#### Custom colorbar
colors = [(5,48,97),(33,102,172),(67,147,195),(146,197,222),(209,229,240),(253,219,199),(244,165,130),(214,96,77),(178,24,43),(103,0,31)]

def make_cmap(colors,position=None,bit=False):
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))
    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,len(colors))
    #cmap = mpl.cm.cool
    return cmap

my_cmap = make_cmap(colors,bit=True) #mpl.cm.cool#
####################################################################################################

####################################################################################################
file               = './4jared/CESM2_spatial_precip_diff.nc'
f1                 = Dataset(file, 'r')
CESM2_timesteps    = f1.variables['year'][:]  # time
CESM2_lat          = f1.variables['latitude'][:]
CESM2_lon          = f1.variables['longitude'][:]
CESM2_sim1         = f1.variables['CESM2_ann_avg_precip1'][:]
CESM2_sim2         = f1.variables['CESM2_ann_avg_precip2'][:]
f1.close()


CESM2_final_sum_diff = ozoneDriver("maize", 20, 1)
#print(CESM2_final_sum_diff[144])
'''
What is the purpose of the [:,0:144] 
I can see it changes the output ever so slightly, but why is it important
And this: [10:40,:,:]
'''
#print(CESM2_final_sum_diff[:,0:144][1][1]) #this is 192 x 144. The print statement as is returns a float.
#this is an array of arrays of floats.

CESM2_final_sum_diff2 = precipDriver("maize", 20, 1)

CESM2_final_sum_diff3 = CO2Driver("maize", 20, 1)

#%% Create an empty figure
fig = plt.figure(figsize=(10,2.2),frameon=False)

#%% fit axes to it
ax_CESM2     = fig.add_axes([0.01,0.15,0.28,0.80])
ax_CESM2_2     = fig.add_axes([0.30,0.15,0.28,0.80])
ax_CESM2_3     = fig.add_axes([0.59,0.15,0.28,0.80])

ax_zonal     = fig.add_axes([0.88,0.225,0.11,0.646])
ax_colorbar1 = fig.add_axes([0.07,0.10,0.743,0.06])

'''
I've never seen a graph created like this in numpy. Is this command specific to basemap?
'''

#%% Put world map on graph
m = Basemap(projection='robin',lon_0=0,resolution='c',ax=ax_CESM2)
m.drawmapboundary(fill_color='#BDBDBD')
m.fillcontinents(color='#BDBDBD',lake_color='#BDBDBD',zorder=0)
m.drawcoastlines(color='#000000', linewidth=0.5)
m.drawcountries(color='#000000', linewidth=0.5)

n = Basemap(projection='robin',lon_0=0,resolution='c',ax=ax_CESM2_2)
n.drawmapboundary(fill_color='#BDBDBD')
n.fillcontinents(color='#BDBDBD',lake_color='#BDBDBD',zorder=0)
n.drawcoastlines(color='#000000', linewidth=0.5)
n.drawcountries(color='#000000', linewidth=0.5)
o = Basemap(projection='robin',lon_0=0,resolution='c',ax=ax_CESM2_3)
o.drawmapboundary(fill_color='#BDBDBD')
o.fillcontinents(color='#BDBDBD',lake_color='#BDBDBD',zorder=0)
o.drawcoastlines(color='#000000', linewidth=0.5)
o.drawcountries(color='#000000', linewidth=0.5)

#%% Get levels for colorbar              
#colorbar_levels = [-16, 1, 2, 4, 5, 6, 7, 7.5, 7.75, 8, 8.5, 9, 9.5, 10, 10.5, 15, 20] #how do I make a continuous grpah?
#colorbar_levels = [-32, -16, -8, -4, -2, 0, 2, 4, 8, 16, 32]
norm = mpl.colors.Normalize(vmin = -1, vmax = 1)
colorbar_levels = [-0.1, -0.075, -0.05, -0.025, 0, 0.025, 0.05, 0.075, 0.1]
#norm = mpl.colors.BoundaryNorm(colorbar_levels,my_cmap.N)
#colorbar_levels = mpl.colorbar.ColorbarBase(cmap=cmap,
                                #norm=norm,
                                #orientation='horizontal')

#norm = mpl.colors.Normalize(vmin=5, vmax=16)
'''
this is deceptively important? Where is the data actually mapped onto the graph?
Nothing will be colored without these lines. But the data is still red into the graphs, right?
Is there no default? Why does this return a monotone gray. Why not grayscale, or some other default? 
'''
#%% Make lat and lon into a meshgrid coordinate pair

CESM2_lon   = CESM2_lon[:] + 0.625
CESM2_lat   = CESM2_lat[:] + 0.4712
x,y         = np.meshgrid(CESM2_lon,CESM2_lat)
#print(y[1][1]) returns -88.58..... y[1] returns an array of all 
kx,ky       = n(x,y)


#%% create meshes
'''
pcolor([X, Y,] C, **kwargs)
'''
map1 = m.pcolormesh(kx,ky,CESM2_final_sum_diff[:,:],cmap=mpl.cm.afmhot_r,vmin=-0.05,vmax=0.05,ax=ax_CESM2,norm=norm)

CESM2_2_lon    = CESM2_lon[:] + 0.625
CESM2_2_lat    = CESM2_lat[:] + 0.4712
x,y         = np.meshgrid(CESM2_lon,CESM2_lat)
lx,ly       = m(x,y)
map1 = n.pcolormesh(lx,ly,CESM2_final_sum_diff2[:,:],cmap=mpl.cm.afmhot_r,vmin=-0.05,vmax=0.05,ax=ax_CESM2_2,norm=norm)
CESM2_3_lon    = CESM2_lon[:] + 0.625
CESM2_3_lat    = CESM2_lat[:] + 0.4712
x,y         = np.meshgrid(CESM2_lon,CESM2_lat)
lx,ly       = m(x,y)
map1 = o.pcolormesh(lx,ly,CESM2_final_sum_diff3[:,:],cmap=mpl.cm.afmhot_r,vmin=-0.05,vmax=0.05,ax=ax_CESM2_3,norm=norm)

#%% Make colorbar legend
cbar1 = plt.colorbar(map1,cax=ax_colorbar1,orientation='horizontal',ticks=colorbar_levels)
cbar1.ax.tick_params(labelsize=10)

#%% Make horizontal lines in zonal mean plot
ax_zonal.plot([-1000,1000],[-60,-60],'#808080',alpha=0.9,linestyle=':',linewidth=1)
ax_zonal.plot([-1000,1000],[-30,-30],'#808080',alpha=0.9,linestyle=':',linewidth=1)
ax_zonal.plot([-1000,1000],[0,0],'#808080',alpha=0.9,linestyle=':',linewidth=1)
ax_zonal.plot([-1000,1000],[30,30],'#808080',alpha=0.9,linestyle=':',linewidth=1)
ax_zonal.plot([-1000,1000],[60,60],'#808080',alpha=0.9,linestyle=':',linewidth=1)
ax_zonal.plot([0,0],[-1000,1000],'k',alpha=0.9,linestyle='-',linewidth=1)

#%% Make zonal mean plot
ax_zonal.plot(np.average(CESM2_final_sum_diff[:,:],axis=1),CESM2_lat[:],'b',alpha=0.9,label='total',linestyle='-',linewidth=2)
ax_zonal.plot(np.average(CESM2_final_sum_diff2[:,:],axis=1),CESM2_lat[:],'g',alpha=0.9,label='temp',linestyle='-',linewidth=2)
ax_zonal.plot(np.average(CESM2_final_sum_diff3[:,:],axis=1),CESM2_lat[:],'r',alpha=0.9,label='CO2',linestyle='-',linewidth=2)

#%% Set x and y limits for zonal mean plot. X units are %change and y units are latitude

ax_zonal.legend(loc=2,numpoints=1,fontsize=4,framealpha=1.0)
ax_zonal.tick_params(labelsize=8)
ax_zonal.axis([-80, 80, -90, 90])
ax_zonal.set_xticks([-1,1])
ax_zonal.set_yticks([])
ax_zonal.yaxis.grid(True,linestyle='--',which='major',color='k',alpha=0.4,zorder=0)
ax_zonal.set_title('Zonal Mean [%]',fontsize=8)

#%%% Make titles
ax_CESM2.set_title('Crop Yield Changes due to Ozone',fontsize=9)
ax_CESM2_2.set_title('Crop Yield Changes due to Precipitation',fontsize=9)
ax_CESM2_3.set_title('Crop Yield Changes due to CO2',fontsize=9)

#%% Save output as PDF
plt.savefig('./OzoneDriver_output.pdf')
plt.show()
####################################################################################################