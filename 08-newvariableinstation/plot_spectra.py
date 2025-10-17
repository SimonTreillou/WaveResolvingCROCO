#%%
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from scipy import signal
import pandas as pd
fname='./08-newvariableinstation/FLASH_RIP/'
fname='./FLASH_RIP/'

#%%
ds = nc.Dataset(fname+'rip_his.nc')
time_fin = ds.variables['scrum_time'][-1]
xr = ds.variables['x_rho'][0,:]
yr = ds.variables['y_rho'][:,0]
dx = xr[1]-xr[0]
u = ds.variables['u'][-1,-1,:,:]
v = ds.variables['v'][-1,-1,:,:]
vorticity=(v[:,1:]-v[:,:-1])/dx - (u[1:,:]-u[:-1,:])/dx
ds.close()

ds = nc.Dataset(fname+'stations.nc')
time = ds.variables['scrum_time'][2:]
print('Total time:'+str(time[-1])+' s')
dt=time[1]-time[0]
X=ds.variables['Xgrid'][:]
Y=ds.variables['Ygrid'][:]
print("Position: ("+str(xr[int(X[0])])+","+str(yr[int(Y[0])])+") m")
zeta=ds.variables['zeta'][2:,:]
vort=ds.variables['vort'][2:,:,:]
ds.close()
lvl=7

#%%
plt.figure(dpi=300)
p=plt.pcolor(xr,yr,vorticity,cmap='RdBu',vmin=-0.1,vmax=0.1)
plt.scatter(xr[int(X[0])],yr[int(Y[0])],color='r',edgecolors='k')
plt.colorbar(p)
plt.xlabel(r'x [m]')
plt.ylabel(r'y [m]')
plt.title(r'$\omega_z$ at $t=$'+str(time_fin)+'s [$s^{-1}$]',loc='left')
plt.show()

#%%
lvl=7
vmax=np.max(np.abs(vort))*1.1;vmin=-vmax
plt.figure(dpi=300,figsize=(10,4))
plt.hlines(0,time.min(),time.max(),'gray',linestyle=':')
plt.plot(time,vort[:,0,lvl],'k')
plt.xlabel('Time [s]')
plt.ylabel(r'$\omega_z$ [$s^{-1}$]')
plt.title(r"$\omega_z$ at Position: ("+str(xr[int(X[0])])+","+str(yr[int(Y[0])])+") m",loc='left')
it = np.argmin(np.abs(time-time_fin))
plt.ylim(vmin,vmax)
plt.scatter(time[it],vort[it,0,lvl],color='r')
plt.show()

#%%
vmax=np.max(np.abs(zeta))*1.1;vmin=-vmax
plt.figure(dpi=300,figsize=(10,4))
plt.hlines(0,time.min(),time.max(),'gray',linestyle=':')
plt.plot(time,zeta[:,0],'k')
plt.xlabel('Time [s]')
plt.ylabel(r'$\eta$ [m]')
plt.title(r"$\eta$ at Position: ("+str(xr[int(X[0])])+","+str(yr[int(Y[0])])+") m",loc='left')
it = np.argmin(np.abs(time-time_fin))
plt.ylim(vmin,vmax)
plt.scatter(time[it],zeta[it,0],color='r')
plt.show()

#%%
N=2048
f,S=signal.csd(vort[:,0,lvl],vort[:,0,lvl],fs=1/dt, nperseg=N,noverlap=N/2)
plt.figure(dpi=300)
plt.loglog(f,S,linewidth=2,color='k')
plt.xlabel(r"Frequency [Hz]")
plt.ylabel(r"$S_{\omega_z \omega_z}$")
plt.vlines(1/13,1e-5,1e0)
plt.show()

#%%
N=512
f,S=signal.csd(zeta[:,0],zeta[:,0],fs=1/dt, nperseg=N,noverlap=N/2)
plt.figure(dpi=300)
plt.loglog(f,S,linewidth=2,color='k')
plt.xlabel(r"Frequency [Hz]")
plt.ylabel(r"$S_{\eta \eta}$")
plt.vlines(1/13,1e-5,1e0)
plt.show()