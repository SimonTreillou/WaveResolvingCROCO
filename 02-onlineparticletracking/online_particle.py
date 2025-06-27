from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

#%% Load
ds = Dataset('./02-onlineparticletracking/rip_avg.nc')
xr = ds.variables['x_rho'][0,:]
yr = ds.variables['y_rho'][:,0]
u = ds.variables['u'][-1,-1,:,:]
u = u2rho(u)
ds.close()

ds = Dataset('./02-onlineparticletracking/floats.nc')
Xgrid=ds.variables['Xgrid'][:]
Ygrid=ds.variables['Ygrid'][:]
depth=ds.variables['depth'][:]
ds.close()


#%% Plot
plt.figure(dpi=300)
pa = plt.pcolor(xr,yr,u,cmap='RdBu',vmin=-0.3,vmax=0.3)
cba = plt.colorbar(pa)

dx=1.5
pb = plt.scatter(Xgrid*dx,Ygrid*dx,c=depth, cmap="plasma",linewidth=0.05, edgecolor='w')
cbb = plt.colorbar(pb)

cba.set_label('velocity')
cbb.set_label('depth')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.savefig('./02-onlineparticletracking/online_particle.png')
plt.show()


#%% Tools
def u2rho(var_u):
    if len(var_u.shape)==4:
        [T,N,Mp,L]=var_u.shape
        Lp=L+1
        Lm=L-1
        var_rho=np.zeros((T,N,Mp,Lp))
        var_rho[:,:,:,1:L-1]=0.5*(var_u[:,:,:,0:Lm-1]+var_u[:,:,:,1:L-1])
        var_rho[:,:,:,0]=var_rho[:,:,:,1]
        var_rho[:,:,:,L]=var_rho[:,:,:,L-1]
    elif len(var_u.shape)==3:
        [N,Mp,L]=var_u.shape
        Lp=L+1
        Lm=L-1
        var_rho=np.zeros((N,Mp,Lp))
        var_rho[:,:,1:L-1]=0.5*(var_u[:,:,0:Lm-1]+var_u[:,:,1:L-1])
        var_rho[:,:,0]=var_rho[:,:,1]
        var_rho[:,:,L]=var_rho[:,:,L-1]
    elif len(var_u.shape)==2:
        [Mp,L]=var_u.shape
        Lp=L+1
        Lm=L-1
        var_rho=np.zeros((Mp,Lp))
        var_rho[:,1:L-1]=0.5*(var_u[:,0:Lm-1]+var_u[:,1:L-1])
        var_rho[:,0]=var_rho[:,1]
        var_rho[:,L]=var_rho[:,L-1]
    return var_rho