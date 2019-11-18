# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 11:55:32 2019

@author: Jared

Jared's version of graphing code
"""
#%% Importing Modules

import matplotlib
matplotlib.use('Agg')
import os
os.environ['PROJ_LIB'] = r'C:\Users\Jared\Anaconda3\pkgs\proj4-5.2.0-ha925a31_1\Library\share'
from pylab import *
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

rc('font',family='serif')

#%% Get file data
file               = './4jared/CESM2_spatial_precip_diff.nc'
f1                 = Dataset(file, 'r')
CESM2_timesteps    = f1.variables['year'][:]  # time
CESM2_lat          = f1.variables['latitude'][:]
CESM2_lon          = f1.variables['longitude'][:]
CESM2_sim1         = f1.variables['CESM2_ann_avg_precip1'][:]
CESM2_sim2         = f1.variables['CESM2_ann_avg_precip2'][:]
f1.close()

#%% Arrange data for graphing
CESM2_final_sum_diff = yieldLoss4Graph("maize", 20, 0.006)
#%% Create an empty figure
fig = plt.figure(figsize=(10,2.2),frameon=False)

#%% fit axes to it
ax_CESM2     = fig.add_axes([0.01,0.15,0.28,0.80])
ax_zonal     = fig.add_axes([0.88,0.225,0.11,0.646])
ax_colorbar1 = fig.add_axes([0.07,0.10,0.743,0.06])

#%% Put world map on graph
m = Basemap(projection='robin',lon_0=0,resolution='c',ax=ax_CESM2)
m.drawmapboundary(fill_color='#BDBDBD')
m.fillcontinents(color='#BDBDBD',lake_color='#BDBDBD',zorder=0)
m.drawcoastlines(color='#000000', linewidth=0.5)
m.drawcountries(color='#000000', linewidth=0.5)
                
#%% Get levels for colorbar              
colorbar_levels = [-32, -16, -8, -4, -2, 0, 2, 4, 8, 16, 32]
norm = mpl.colors.Normalize(vmin = -40, vmax = 40)

#%% Make lat and lon into a meshgrid coordinate pair

CESM2_lon   = CESM2_lon[:] + 0.625
CESM2_lat   = CESM2_lat[:] + 0.4712
x,y         = np.meshgrid(CESM2_lon,CESM2_lat)
#print(y[1][1]) returns -88.58..... y[1] returns an array of all 
kx,ky       = n(x,y)

#%% create meshes

map1 = m.pcolormesh(kx,ky,CESM2_final_sum_diff[:,:],cmap=mpl.cm.afmhot_r,vmin=0,vmax=16,ax=ax_CESM2,norm=norm)
#I personall prefer seismic, becuase that's what Dr. Shindell used.
#%% Make colorbar legend
cbar1 = plt.colorbar(map1,cax=ax_colorbar1,orientation='horizontal',ticks=colorbar_levels)
cbar1.ax.tick_params(labelsize=10)

#%% Make horizontal lines in zonal mean plot
ax_zonal.plot([-1000,1000],[-60,-60],'#808080',alpha=0.9,linestyle=':',linewidth=1)
ax_zonal.plot([-1000,1000],[-30,-30],'#808080',alpha=0.9,linestyle=':',linewidth=1)
        
#%% Make zonal mean plot
ax_zonal.plot(np.average(CESM2_final_sum_diff[:,:],axis=1),CESM2_lat[:],'b',alpha=0.9,label='CESM2',linestyle='-',linewidth=2)

#%% Set x and y limits for zonal mean plot. X units are %change and y units are latitude

ax_zonal.legend(loc=2,numpoints=1,fontsize=4,framealpha=1.0)
ax_zonal.tick_params(labelsize=8)
ax_zonal.axis([-80, 80, -90, 90])
ax_zonal.set_xticks([-40,40])
ax_zonal.set_yticks([])
ax_zonal.yaxis.grid(True,linestyle='--',which='major',color='k',alpha=0.4,zorder=0)
ax_zonal.set_title('Zonal Mean [%]',fontsize=8)

#%%% Make titles
ax_CESM2.set_title('Combined Yield Loss From all Sources',fontsize=9)  
plt.savefig('./total_yield_change_output.pdf')            

