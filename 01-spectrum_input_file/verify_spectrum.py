import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from scipy import signal
import pandas as pd

ds = nc.Dataset('./01-spectrum_input_file/FLASH_RIP/stations.nc')
time = ds.variables['scrum_time'][:]
print('Total time:'+str(time[-1])+' s')
dt=time[1]-time[0]
zeta=ds.variables['zeta'][:,0]

f,S=signal.csd(zeta,zeta,fs=1/dt, nperseg=512,noverlap=256)
df=pd.read_csv('./01-spectrum_input_file/FLASH_RIP/input_spectrum.txt',header=None,sep=' ')
dfrq=f[1]-f[0]
Hs_mod = 4*np.sqrt(np.sum(S*dfrq))

plt.figure(dpi=300)
plt.loglog(df[0],df[1],linewidth=5,alpha=0.5,label=r"Input")
plt.loglog(f,S,linewidth=2,label=r"Model")
plt.legend()
plt.xlabel(r"Frequency [Hz]")
plt.ylabel(r"$S_{\eta\eta}$")
#plt.xlim([0,0.5])
plt.ylim([1e-4,1e1])
plt.show()

plt.figure(dpi=300)
plt.plot(df[0],df[1],linewidth=5,alpha=0.5,label=r"Input")
plt.plot(f,S,linewidth=2,label=r"Model")
plt.legend()
plt.xlabel(r"Frequency [Hz]")
plt.ylabel(r"$S_{\eta\eta}$")
plt.xlim([0,0.5])
plt.show()