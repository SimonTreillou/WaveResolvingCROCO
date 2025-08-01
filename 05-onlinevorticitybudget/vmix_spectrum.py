from netCDF4 import Dataset
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/simon/Code/Diagnostics/Tools')
from spectrum import welch_spectrum_CI, plot_CIbar_loglog
import numpy as np

cas = "FLASH_RIP"
xpt=160
ypt=150
N=256
tini=200

ds = Dataset('./05-onlinevorticitybudget/'+cas+'/rip_his.nc')
h=ds.variables['h'][0,:]
ih0=np.argmin(np.abs(h))
X=ds.variables['x_rho'][0,xpt] - ds.variables['x_rho'][0,ih0]
Y=ds.variables['y_rho'][ypt,0]
ds.close()

ds = Dataset('./05-onlinevorticitybudget/'+cas+'/rip_diags_vrt_avg.nc')
time = ds.variables['scrum_time'][tini:]
dt=time[1]-time[0]
vmix = ds.variables['vrt_vmix'][tini:,:,:]
vmix_tseries = vmix[:,ypt,xpt]
xadv = ds.variables['vrt_xadv'][tini:,:,:]
xadv_tseries = xadv[:,ypt,xpt]
yadv = ds.variables['vrt_yadv'][tini:,:,:]
yadv_tseries = yadv[:,ypt,xpt]
prsgrd = ds.variables['vrt_prsgrd'][tini:,:,:]
prsgrd_tseries = prsgrd[:,ypt,xpt]
rate = ds.variables['vrt_rate'][tini:,:,:]
rate_tseries = rate[:,ypt,xpt]

plt.figure(dpi=300)
f,S,dof,ci_lower,ci_upper=welch_spectrum_CI(vmix_tseries,fs=1/dt,nperseg=N,alpha=0.05)
plt.loglog(f,S,label='V. mix.')
Smax=np.max(S);Smin=np.min(S)
plot_CIbar_loglog(f,S,ci_upper)

f,S,dof,ci_lower,ci_upper=welch_spectrum_CI(prsgrd_tseries,fs=1/dt,nperseg=N,alpha=0.05)
plt.loglog(f,S,label='Pres. grad.')
if np.max(S)>Smax: Smax=np.max(S)
if np.min(S)<Smin: Smin=np.min(S)

f,S,dof,ci_lower,ci_upper=welch_spectrum_CI(xadv_tseries,fs=1/dt,nperseg=N,alpha=0.05)
plt.loglog(f,S,label="x-Advection")
if np.max(S)>Smax: Smax=np.max(S)
if np.min(S)<Smin: Smin=np.min(S)

f,S,dof,ci_lower,ci_upper=welch_spectrum_CI(yadv_tseries,fs=1/dt,nperseg=N,alpha=0.05)
plt.loglog(f,S,label="y-Advection")
if np.max(S)>Smax: Smax=np.max(S)
if np.min(S)<Smin: Smin=np.min(S)

f,S,dof,ci_lower,ci_upper=welch_spectrum_CI(rate_tseries,fs=1/dt,nperseg=N,alpha=0.05)
plt.loglog(f,S,label="Rate")
if np.max(S)>Smax: Smax=np.max(S)
if np.min(S)<Smin: Smin=np.min(S)

plt.vlines(1/13,Smin,Smax,color='grey',linestyle='--')
plt.text(1/12,Smin*1.3,'Peak period',rotation=90,color='gray')
plt.ylim(Smin*0.9,Smax*1.2)
plt.legend()
plt.xlabel('Frequency [Hz]')
plt.ylabel('PSD of Vorticity terms')
plt.title('PSD of vorticity terms at (x='+str(int(X))+',y='+str(int(Y))+') m',loc='left')
plt.show()