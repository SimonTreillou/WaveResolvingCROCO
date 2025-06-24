import numpy as np
from opendrift.readers import reader_ROMS_native, reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift
from netCDF4 import Dataset
from pyproj import Proj, Transformer
import pyproj
pyproj.datadir.set_data_dir('/opt/anaconda3/envs/opendrift2/share/proj') # necessary
import matplotlib.pyplot as plt
from opendrift_croco_tools import *

# === Path to your grid file ===
fname = "/Users/simon/Code/CONFIGS/IB09_SZ/IB09_3D_TRC_Dcrit_WMPER_test/rip_avg_original.nc"
add_lonlat_croco(fname)
add_mask_croco(fname)

#%%
o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information

croco = reader_ROMS_native.Reader(fname)
o.add_reader(croco)
#o.set_config('drift:horizontal_diffusivity', 0.0001)     # Difference from first run

lon_seed,lat_seed=xy_to_lonlat([280,280,280,280,280],[100,200,300,400,500])
o.seed_elements(lon=lon_seed, lat=lat_seed, radius=0, number=5,
                z=[0,0,0,0,0], time=croco.start_time)

o.run(steps=160, time_step=50,outfile='res.nc')

#print(o)
#o.plot(linecolor='z', fast=True)

#%%
nc = Dataset(fname)
lon = nc.variables['lon_rho'][:]
lat = nc.variables['lat_rho'][:]
x,y=lonlat_to_xy(lon,lat)
u = nc.variables['u'][-1,0,:,:]
u = u2rho(u)

part = Dataset('res.nc')
nb_particles = part.variables['lon'].shape[0]

plt.pcolor(x,y,u,cmap='RdBu')
for i in range(nb_particles):
    plt.plot(part.variables['lon'][i,0],part.variables['lat'][i,0],'*',markersize=20,color='white')
    plt.plot(part.variables['lon'][i,:],part.variables['lat'][i,:],'-o',markersize=2,color='green')
plt.colorbar()
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.clim(-0.3,0.3)
plt.title(r'Cross-shore velocity $u$',loc='left')
plt.show()