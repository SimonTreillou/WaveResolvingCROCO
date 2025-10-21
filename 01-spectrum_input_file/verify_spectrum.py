#%%
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from scipy import signal
import pandas as pd
import sys
sys.path.append('/Users/simon/Code/Projects/Tools')
from Diagnostics.spectrum import welch_spectrum_CI,plot_CIbar_loglog # pyright: ignore[reportMissingImports]

#%% Load
ds = nc.Dataset('./FLASH_RIP/stations.nc')
time = ds.variables['scrum_time'][:]
print('Total time:'+str(time[-1])+' s')
dt=time[1]-time[0]
zeta=ds.variables['zeta'][:,0]
ds.close()

#%% Compute
N=1024
f,S,dof,ci_lower,ci_upper=welch_spectrum_CI(zeta,fs=1/dt, nperseg=N,noverlap=N/2)
df=pd.read_csv('./FLASH_RIP/input_spectrum.txt',header=None,sep=' ')
dfrq=f[1]-f[0]
Hs_mod = 4*np.sqrt(np.sum(S*dfrq))

#%% Plot in log-log
plt.figure(dpi=300)
plt.loglog(df[0],df[1],linewidth=5,alpha=0.5,label=r"Input",color='blue')
plt.loglog(f,S,linewidth=2,label=r"Model",color='red')
plot_CIbar_loglog(f,S,ci_upper,color='red',y_pos=1e-2,lw=2)
plt.legend()
plt.xlabel(r"Frequency [Hz]")
plt.ylabel(r"$S_{\eta\eta}$")
#plt.xlim([0,0.5])
plt.ylim([1e-4,1e1])
plt.show()

#%% Plot in linear
plt.figure(dpi=300)
plt.plot(df[0],df[1],linewidth=5,alpha=0.5,label=r"Input")
plt.plot(f,S,linewidth=2,label=r"Model")
plt.legend()
plt.xlabel(r"Frequency [Hz]")
plt.ylabel(r"$S_{\eta\eta}$")
plt.xlim([0,0.5])
plt.show()
# %%
