#%%
import numpy as np
import xarray as xr
import netCDF4 as nc
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/simon/Code/Projects/Tools')
from Diagnostics.grid import u2rho,v2rho
from Diagnostics.useful import wave_average

#%%
# === Load CROCO file ===
ds = nc.Dataset("./FLASH_RIP/rip_avg.nc")
ds_his = nc.Dataset("./FLASH_RIP/rip_his.nc")
ubar = ds.variables['ubar'][:,:,:]
ubar = u2rho(ubar)
vbar = ds.variables['vbar'][:,:,:]
vbar = v2rho(vbar)
ubarhis = ds_his.variables['ubar'][:,:,:]
ubarhis = u2rho(ubarhis)
vbarhis = ds_his.variables['vbar'][:,:,:]
vbarhis = v2rho(vbarhis)
h = ds.variables['h'][:,:]
zeta = ds.variables['zeta'][:,:,:]
zetahis = ds_his.variables['zeta'][:,:,:]
x = ds.variables['x_rho'][0,:]
y = ds.variables['y_rho'][:,0]
dx = x[1]-x[0]
dy = y[1]-y[0]

#%% Compute vorticity
vort = v2rho((ubarhis[:,1:,:]-ubarhis[:,:-1,:])/dy) - u2rho((vbarhis[:,:,:1]-vbarhis[:,:,:-1])/dx)
vort_avg=wave_average(vort)

plt.pcolor(vort_avg[-1,:,:], cmap="RdBu")
plt.clim([-0.3,0.3])
plt.colorbar()
plt.show()

#%% Compute vorticity rate
dt=13
rate = (vort_avg[1:,:,:]-vort_avg[:-1,:,:])/dt

plt.pcolor(rate[-1,:,:], cmap="RdBu")
plt.clim([-0.01,0.01])
plt.colorbar()
plt.show()

#%% Compute Fw
Cf=0.01
urms = (ubarhis**2 + vbarhis**2)**0.5
beta = Cf*wave_average(urms)

cff1 = beta/h * vort_avg
cff2 = beta/h**2 * ( ubar * v2rho((h[1:,:]-h[:-1,:])/dy) -   vbar * u2rho((h[:,1:]-h[:,:-1])/dx) )
Fw = -cff1 -cff2

plt.pcolor(Fw[-1,:,:], cmap="RdBu")
plt.clim([-0.02,0.02])
plt.colorbar()
plt.show()

#%% Compute fluxes
# F = rho*h*u*(g*eta + U^2/2)
g = 9.81
rho = 1024.0
U = np.sqrt(ubarhis**2 + vbarhis**2)
F = h*ubarhis*(g*zetahis+U**2/2)
G = h*vbarhis*(g*zetahis+U**2/2)

plt.pcolor(F[-10,:,:],cmap='BuGn')
#plt.clim([0,0.05])
plt.colorbar()
plt.show()

#%% Compute dissipation D
dFdx = u2rho((F[:,:,1:]-F[:,:,:-1])/dx)
dGdy = v2rho((G[:,1:,:]-G[:,:-1,:])/dy)
deltabm = - np.nanmean((dFdx +dGdy),0)
Cphi = (h*g)**0.5 # not really, should be the shock velocity
D = deltabm / ( h *Cphi)

plt.pcolor(D,cmap='BuGn')
plt.clim([0,0.05])
plt.colorbar()
plt.show()

#%%
ds.variables['tke']

#%% Compute Pw
dDdx = u2rho((D[:,1:]-D[:,:-1])/dx)
dDdy = v2rho((D[1:,:]-D[:-1,:])/dy)
costheta = F / (F**2 + G**2)**0.5
sintheta = G / (F**2 + G**2)**0.5
Pw = dDdx*sintheta - dDdy*costheta

plt.pcolor(Pw[-1,:,:],cmap='RdBu')
plt.clim([-0.02,0.02])
plt.colorbar()
plt.show()

#%% Compute Aw
cff1 = ubar * u2rho((vort_avg[:,:,1:]-vort_avg[:,:,:-1])/dx)
cff2 = vbar * v2rho((vort_avg[:,1:,:]-vort_avg[:,:-1,:])/dy)
cff3 = vort_avg * ( u2rho((ubar[:,:,1:]-ubar[:,:,:-1])/dx)  +  v2rho((vbar[:,1:,:]-vbar[:,:-1,:])/dy) )
Aw = -cff1 -cff2 -cff3

plt.pcolor(Aw[-1,:,:],cmap='RdBu')
#plt.clim([-0.2,0.2])
plt.colorbar()
plt.show()

#%% Compute Mw
Mw=rate-Aw[:-1,:,:]-Pw[:-1,:,:]-Fw[:-1,:,:]

#%%
plt.plot(np.mean(np.abs(Aw),axis=(0, 1)),label=r"$<|A_w|>_y$")
plt.plot(np.mean(np.abs(Pw),axis=(0, 1)),label=r"$<|P_w|>_y$")
plt.plot(np.mean(np.abs(Fw),axis=(0, 1)),label=r"$<|F_w|>_y$")
plt.plot(np.mean(np.abs(Mw),axis=(0, 1)),label=r"$<|M_w|>_y$")
plt.legend()
plt.xlabel(r"x (m)")
plt.ylabel(r"Vorticity eq. terms ($1/s^2$)")
plt.xlim([0,200])
plt.ylim([0,0.005])
plt.show()

#%%
plt.plot(np.mean(vort**2,axis=(0, 1)),label=r"$<\omega^2>_y$")
plt.plot(np.mean(D,axis=(0)),label=r"$<D>_y$")

plt.legend()
plt.xlabel(r"x (m)")
plt.ylabel(r"Vorticity eq. terms ($1/s^2$)")
plt.xlim([0,200])
plt.ylim([0,0.03])
plt.show()

#%%
plt.plot(np.mean(rate-Aw[:-1,:,:]-Pw[:-1,:,:]-Fw[:-1,:,:],axis=(0,1)),label=r"$<|A_w|>_y$")
plt.legend()
plt.xlabel(r"x (m)")
plt.ylabel(r"Vorticity eq. terms ($1/s^2$)")
plt.xlim([0,200])
plt.ylim([-0.005,0.005])
plt.show()