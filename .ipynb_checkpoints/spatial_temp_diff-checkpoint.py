import os
os.environ['PROJ_LIB'] = r'C:\Users\Jared\Anaconda3\pkgs\proj4-5.2.0-ha925a31_1\Library\share'
import matplotlib
from pylab import *
from mpl_toolkits.basemap import Basemap #Basemap isn't working
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

rc('font',family='serif')
####################################################################################################
#### Custom colorbar
colors = [(222,235,247),(247,251,255),(255,245,240),(254,224,210),(252,187,161),(252,146,114),(251,106,74),(239,59,44),(203,24,29),(165,15,21),(103,0,13)]

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
    return cmap

my_cmap = make_cmap(colors,bit=True)
####################################################################################################

####################################################################################################
file               = '4jared/CESM2_spatial_temp_diff.nc'
f1                 = Dataset(file, 'r')
CESM2_timesteps    = f1.variables['year'][:]  # time
CESM2_lat          = f1.variables['latitude'][:]
CESM2_lon          = f1.variables['longitude'][:]
CESM2_sim1         = f1.variables['CESM2_ann_avg_temp1'][:]
CESM2_sim2         = f1.variables['CESM2_ann_avg_temp2'][:]
f1.close()
####################################################################################################

####################################################################################################
file               = '4jared/GFDL_spatial_temp_diff.nc'
f1                 = Dataset(file, 'r')
GFDL_timesteps     = f1.variables['year'][:]  # time
GFDL_lat           = f1.variables['latitude'][:]
GFDL_lon           = f1.variables['longitude'][:]
GFDL_sim1          = f1.variables['GFDL_ann_avg_temp1'][:]
GFDL_sim2          = f1.variables['GFDL_ann_avg_temp2'][:]
f1.close()
####################################################################################################

####################################################################################################
file               = '4jared/GISS_spatial_temp_diff.nc'
f1                 = Dataset(file, 'r')
GISS_timesteps     = f1.variables['year'][:]  # time
GISS_lat           = f1.variables['latitude'][:]
GISS_lon           = f1.variables['longitude'][:]
GISS_sim1          = f1.variables['GISS_ann_avg_temp1'][:]
GISS_sim2          = f1.variables['GISS_ann_avg_temp2'][:]
f1.close()
####################################################################################################

####################################################################################################
CESM2_mean_diff             = CESM2_sim1[10:40,:,:] - CESM2_sim2[10:40,:,:]
CESM2_mean_diff             = np.average(CESM2_mean_diff[:,:,:],axis=0)
CESM2_final_diff            = np.zeros((len(CESM2_lat),len(CESM2_lon)))

CESM2_lon                   = np.arange(-180,180,1.25)
CESM2_final_diff[:,144:]    = CESM2_mean_diff[:,0:144]
CESM2_final_diff[:,0:144]   = CESM2_mean_diff[:,144:]

GFDL_mean_diff              = GFDL_sim1[10:40,:,:] - GFDL_sim2[10:40,:,:]
GFDL_mean_diff              = np.average(GFDL_mean_diff[:,:,:],axis=0)
GFDL_final_diff             = np.zeros((len(GFDL_lat),len(GFDL_lon)))

GFDL_lon                    = np.arange(-180,180,1.25)
GFDL_final_diff[:,144:]     = GFDL_mean_diff[:,0:144]
GFDL_final_diff[:,0:144]    = GFDL_mean_diff[:,144:]

GISS_mean_diff              = GISS_sim1[10:40,:,:] - GISS_sim2[10:40,:,:]
GISS_final_diff             = np.average(GISS_mean_diff[:,:,:],axis=0)
####################################################################################################

####################################################################################################
fig = plt.figure(figsize=(10,2.2),frameon=False)

ax_CESM2     = fig.add_axes([0.01,0.15,0.28,0.80])
ax_GFDL      = fig.add_axes([0.30,0.15,0.28,0.80])
ax_GISS      = fig.add_axes([0.59,0.15,0.28,0.80])
ax_zonal     = fig.add_axes([0.88,0.225,0.11,0.646])
ax_colorbar1 = fig.add_axes([0.07,0.10,0.743,0.06])

m = Basemap(projection='robin',lon_0=0,resolution='c',ax=ax_CESM2)
m.drawmapboundary(fill_color='#BDBDBD')
m.fillcontinents(color='#BDBDBD',lake_color='#BDBDBD',zorder=0)
m.drawcoastlines(color='#000000', linewidth=0.5)
m.drawcountries(color='#000000', linewidth=0.5)

n = Basemap(projection='robin',lon_0=0,resolution='c',ax=ax_GFDL)
n.drawmapboundary(fill_color='#BDBDBD')
n.fillcontinents(color='#BDBDBD',lake_color='#BDBDBD',zorder=0)
n.drawcoastlines(color='#000000', linewidth=0.5)
n.drawcountries(color='#000000', linewidth=0.5)

o = Basemap(projection='robin',lon_0=0,resolution='c',ax=ax_GISS)
o.drawmapboundary(fill_color='#BDBDBD')
o.fillcontinents(color='#BDBDBD',lake_color='#BDBDBD',zorder=0)
o.drawcoastlines(color='#000000', linewidth=0.5)
o.drawcountries(color='#000000', linewidth=0.5)

colorbar_levels = [-0.20,-0.10,0.0,0.10,0.20,0.30,0.40,0.50,0.60]
norm = mpl.colors.BoundaryNorm(colorbar_levels,my_cmap.N)

CESM2_lon   = CESM2_lon[:] + 0.625
CESM2_lat   = CESM2_lat[:] + 0.4712
x,y         = np.meshgrid(CESM2_lon,CESM2_lat)
kx,ky       = n(x,y)

map1 = m.pcolormesh(kx,ky,CESM2_final_diff[:,:],cmap=my_cmap,vmin=-0.20,vmax=0.60,ax=ax_CESM2,norm=norm)

GFDL_lon    = GFDL_lon[:] + 0.625
GFDL_lat    = GFDL_lat[:] + 1.0
x,y         = np.meshgrid(GFDL_lon,GFDL_lat)
lx,ly       = m(x,y)

map1 = n.pcolormesh(lx,ly,GFDL_final_diff[:,:],cmap=my_cmap,vmin=-0.20,vmax=0.60,ax=ax_GFDL,norm=norm)

GISS_lon    = GISS_lon[:] + 1.25
GISS_lat    = GISS_lat[:] + 1.0
x,y         = np.meshgrid(GISS_lon,GISS_lat)
lx,ly       = m(x,y)

map1 = o.pcolormesh(lx,ly,GISS_final_diff[:,:],cmap=my_cmap,vmin=-0.20,vmax=0.60,ax=ax_GISS,norm=norm)

cbar1 = plt.colorbar(map1,cax=ax_colorbar1,orientation='horizontal',ticks=colorbar_levels)
cbar1.ax.tick_params(labelsize=10)

ax_zonal.plot([-1000,1000],[-60,-60],'#808080',alpha=0.9,linestyle=':',linewidth=1)
ax_zonal.plot([-1000,1000],[-30,-30],'#808080',alpha=0.9,linestyle=':',linewidth=1)
ax_zonal.plot([-1000,1000],[0,0],'#808080',alpha=0.9,linestyle=':',linewidth=1)
ax_zonal.plot([-1000,1000],[30,30],'#808080',alpha=0.9,linestyle=':',linewidth=1)
ax_zonal.plot([-1000,1000],[60,60],'#808080',alpha=0.9,linestyle=':',linewidth=1)
ax_zonal.plot([0,0],[-1000,1000],'k',alpha=0.9,linestyle='-',linewidth=1)


ax_zonal.plot(np.average(CESM2_final_diff[:,:],axis=1),CESM2_lat[:],'b',alpha=0.9,label='CESM2',linestyle='-',linewidth=2)
ax_zonal.plot(np.average(GFDL_final_diff[:,:],axis=1),GFDL_lat[:],'g',alpha=0.9,label='GFDL',linestyle='-',linewidth=2)
ax_zonal.plot(np.average(GISS_final_diff[:,:],axis=1),GISS_lat[:],'r',alpha=0.9,label='GISS',linestyle='-',linewidth=2)


ax_zonal.legend(loc=2,numpoints=1,fontsize=4,framealpha=1.0)
ax_zonal.tick_params(labelsize=8)
ax_zonal.axis([-0.8, 0.8, -90, 90])
ax_zonal.set_xticks([-0.4,0.4])
ax_zonal.set_yticks([])
ax_zonal.yaxis.grid(True,linestyle='--',which='major',color='k',alpha=0.4,zorder=0)
ax_zonal.set_title('Zonal Mean [K]',fontsize=8)

ax_CESM2.set_title('CESM2 $\Delta$ Annual Avg. Temp [K]',fontsize=9)
ax_GFDL.set_title('GFDL $\Delta$ Annual Avg. Temp [K]',fontsize=9)
ax_GISS.set_title('GISS $\Delta$ Annual Avg. Temp [K]',fontsize=9)

plt.savefig('./spatial_temp_diff.pdf')
####################################################################################################