from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

cas = "FLASH_RIP"
#cas = "FLASH_RIP_nodrag"
#cas = "FLASH_RIP_noTKE"
#cas = "MINI_RIP"

#lab={'Xadv':""}

# d/dt * curl of (u,v)
ds = Dataset('./07-modifiedonlinevorticitybudget/'+cas+'/rip_avg.nc')
time = ds.variables['scrum_time'][:]
dt = time[1] - time[0]
xr = ds.variables['x_rho'][0,:]
h = ds.variables['h'][0,:]
hu = (h[1:]+h[:-1])/2.0
xu = (xr[1:]+xr[:-1])/2.0
dx=xu[1]-xu[0]
u = ds.variables['u'][:,:,:,:]
v = ds.variables['v'][:,:,:,:]
dvdx=(v[:,:,:,1:]-v[:,:,:,:-1])/dx
dudy=(u[:,:,1:,:]-u[:,:,:-1,:])/dx
vort = dvdx-dudy
vort = np.trapz(vort*hu/u.shape[1],axis=1)
rate_vort=(vort[1:,:,:]-vort[:-1,:,:])/dt
ds.close()

ds = Dataset('./07-modifiedonlinevorticitybudget/'+cas+'/rip_his.nc')
time = ds.variables['scrum_time'][:]
dt = time[1] - time[0]
xr = ds.variables['x_rho'][0,:]
h = ds.variables['h'][0,:]
hu = (h[1:]+h[:-1])/2.0
xu = (xr[1:]+xr[:-1])/2.0
dx=xu[1]-xu[0]
u = ds.variables['u'][:,:,:,:]
v = ds.variables['v'][:,:,:,:]
dzu = hu/u.shape[1]
dzv = h/v.shape[1]
rateu= 0.5*dzu*(u[1:,:,:,:]-u[:-1,:,:,:])/(dt)
rateu = np.sum(rateu,axis=1)
ratev= 0.5*dzv*(v[1:,:,:,:]-v[:-1,:,:,:])/(dt)
ratev = np.sum(ratev,axis=1)
ratevort = (ratev[:,:,1:]-ratev[:,:,:-1])/dx - (rateu[:,1:,:]-rateu[:,:-1,:])/dx
ds.close()

ds = Dataset('./07-modifiedonlinevorticitybudget/'+cas+'/rip_his.nc')
time = ds.variables['scrum_time'][:]
dt = time[1] - time[0]
xr = ds.variables['x_rho'][0,:]
h = ds.variables['h'][0,:]
hu = (h[1:]+h[:-1])/2.0
xu = (xr[1:]+xr[:-1])/2.0
dx=xu[1]-xu[0]
u = ds.variables['ubar'][:,:,:]
v = ds.variables['vbar'][:,:,:]
rateu= (u[1:,:,:]-u[:-1,:,:])/dt
ratev= (v[1:,:,:]-v[:-1,:,:])/dt
ratevort = (ratev[:,:,1:]-ratev[:,:,:-1])/dx - (rateu[:,1:,:]-rateu[:,:-1,:])/dx
ds.close()


ds = Dataset('./05-onlinevorticitybudget/'+cas+'/rip_diags_vrt.nc')
time = ds.variables['scrum_time'][:]
vrate = ds.variables['vrt_rate'][:,:,:]
ds.close()

tini=5
print("Spin-up: "+str(time[tini])+" s")






plt.figure(dpi=300,figsize=(6,6))
plt.plot(xu,np.mean(np.abs(rate_vort[200,:,:]),axis=(0)),label="Rate offline")
plt.plot(xu,np.mean(np.abs(ratevort[200,:,:]),axis=(0)),label="Rate offline")
plt.plot(xu,np.mean(np.abs(vrate[200,:,:]),axis=(0)),label="Rate online")
plt.legend()
plt.xlabel('x [m]')
plt.ylabel(r'$\partial \omega / \partial t$ [m]')
plt.show()