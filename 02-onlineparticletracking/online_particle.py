#%%
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('/Users/simon/Code/Projects/Tools')
from Diagnostics.grid import u2rho,v2rho # pyright: ignore[reportMissingImports]
from Diagnostics.spectrum import welch_spectrum_CI,plot_CIbar_loglog # pyright: ignore[reportMissingImports]

#%% Load
ds = Dataset('./rip_avg.nc')
xr = ds.variables['x_rho'][0,:]
dx = xr[1]-xr[0]
yr = ds.variables['y_rho'][:,0]
u = ds.variables['u'][-1,-1,:,:]
u = u2rho(u)
ds.close()

ds = Dataset('./floats.nc')
Xgrid=ds.variables['Xgrid'][:]
Ygrid=ds.variables['Ygrid'][:]
depth=ds.variables['depth'][:]
ds.close()

#%% Plot
plt.figure(dpi=300)
pa = plt.pcolor(xr,yr,u,cmap='RdBu',vmin=-0.3,vmax=0.3)
cba = plt.colorbar(pa)
pb = plt.scatter(Xgrid*dx,Ygrid*dx,c=depth, cmap="plasma",linewidth=0.05, edgecolor='None')
cbb = plt.colorbar(pb)
cba.set_label('velocity')
cbb.set_label('depth')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.savefig('./online_particle.png')
plt.show()