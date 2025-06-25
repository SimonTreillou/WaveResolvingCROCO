from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

ds = Dataset('./05-onlinevorticitybudget/FLASH_RIP/rip_avg.nc')
xr = ds.variables['x_rho'][0,:]
xu = (xr[1:]+xr[:-1])/2.0
ds.close()

ds = Dataset('./05-onlinevorticitybudget/FLASH_RIP/rip_diags_vrt_avg.nc')
time = ds.variables['scrum_time'][:]
tini=7
print("Spin-up: "+str(time[tini])+" s")
xadv = np.mean(np.abs(ds.variables['vrt_xadv'][tini:,:,:]),(0,1))
yadv = np.mean(np.abs(ds.variables['vrt_yadv'][tini:,:,:]),(0,1))
#hdiff = np.mean(ds.variables['vrt_hdiff'][tini:,:,:],(0,1)) # 0 in this simulation
#cor = np.mean(ds.variables['vrt_cor'][tini:,:,:],(0,1)) # 0 in this simulation
prsgrd = np.mean(np.abs(ds.variables['vrt_prsgrd'][tini:,:,:]),(0,1))
#hmix = np.mean(ds.variables['vrt_hmix'][tini:,:,:],(0,1)) # only in the sponge layer offshore
vmix = np.mean(np.abs(ds.variables['vrt_vmix'][tini:,:,:]),(0,1))
rate = np.mean(np.abs(ds.variables['vrt_rate'][tini:,:,:]),(0,1))
#nudg = np.mean(ds.variables['vrt_nudg'][tini:,:,:],(0,1)) # 0 in this simulation
#wind = np.mean(ds.variables['vrt_wind'][tini:,:,:],(0,1)) # 0 in this simulation
drag = np.mean(np.abs(ds.variables['vrt_drag'][tini:,:,:]),(0,1))
fast = np.mean(np.abs(ds.variables['vrt_fast'][tini:,:,:]),(0,1))
vars_rhs = ['vrt_xadv','vrt_yadv','vrt_hdiff','vrt_cor','vrt_prsgrd','vrt_hmix','vrt_vmix','vrt_nudg','vrt_wind','vrt_drag','vrt_fast']
rhs = 0
for var in vars_rhs:
    rhs = rhs + np.mean(ds.variables[var][tini:,:,:],(0,1))

plt.figure(1,dpi=100)
plt.plot(xu,np.mean(ds.variables['vrt_rate'][tini:,:,:],(0,1))-rhs)
plt.xlabel(r'$x$ (m)')
plt.ylabel(r'$\partial \omega / \partial t - rhs$')
plt.title('Checking budget closure',loc='left')
plt.savefig('./05-onlinevorticitybudget/checking_budget.png')
plt.show()

plt.figure(2,dpi=300)
#plt.plot(xu,rate,label=r"$\partial \omega / \partial t$")
plt.plot(xu,drag,label="Drag")
#plt.plot(xu,fast,label="Fast")
plt.plot(xu,vmix,label="Vertical mix.")
plt.plot(xu,prsgrd,label="Prs. grad.")
plt.plot(xu,xadv,label="x adv.")
#plt.plot(xu,yadv,label="y adv.")
plt.xlabel(r'$x$ (m)')
plt.ylabel(r'Vorticity equation terms ($1/s^2$)')
plt.title('Vorticity budget main term magnitudes',loc='left')
plt.legend()
plt.xlim([100,270])
plt.ylim([-0.00001,0.003])
plt.savefig('./05-onlinevorticitybudget/vorticity_budget.png')
plt.show()


