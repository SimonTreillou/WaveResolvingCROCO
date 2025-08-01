from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

cas = "FLASH_RIP"
#cas = "FLASH_RIP_nodrag"
#cas = "FLASH_RIP_noTKE"
#cas = "MINI_RIP"

#lab={'Xadv':""}

ds = Dataset('./05-onlinevorticitybudget/'+cas+'/rip_avg.nc')
xr = ds.variables['x_rho'][0,:]
h = ds.variables['h'][0,:]
xu = (xr[1:]+xr[:-1])/2.0
u = np.mean(ds.variables['u'][:,:,:,130:150],axis=(0,2,3))

plt.figure(dpi=100,figsize=(5,8))
plt.plot(u,range(8))
plt.vlines(0,0,8,color='gray',linestyle='dashed')
plt.xlabel(r'$<u>_{sz}$ (m/s)')
plt.ylabel(r'$Depth$ (normalized)')
plt.xlim([-0.3,0.3])
plt.ylim([0,7])
plt.title('Surfzone-averaged undertow',loc='left')
plt.savefig('./05-onlinevorticitybudget/undertow_'+cas+'.png')
plt.show()
ds.close()

ds = Dataset('./05-onlinevorticitybudget/'+cas+'/rip_diags_vrt_avg.nc')
time = ds.variables['scrum_time'][:]
tini=5
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
vars_rhs = ['vrt_xadv','vrt_yadv','vrt_cor','vrt_prsgrd','vrt_hmix','vrt_vmix','vrt_nudg','vrt_wind','vrt_fast']
closure = ds.variables['vrt_rate'][tini:,:,:]
for var in vars_rhs:
    closure = closure - ds.variables[var][tini:,:,:]
ds.close()

plt.figure(1,dpi=100)
plt.plot(xu,np.max(np.abs(closure),(0,1)))
plt.xlabel(r'$x$ (m)')
plt.ylabel(r'$\max (\partial \omega / \partial t - rhs)$')
plt.title('Checking budget closure',loc='left')
plt.savefig('./05-onlinevorticitybudget/checking_budget_'+cas+'.png')
plt.show()

plt.figure(2,dpi=300)
#plt.plot(xu,rate,label=r"$\partial \omega / \partial t$")
plt.plot(xu,drag,label="Drag") # already included in vmix
#plt.plot(xu,fast,label="Fast")
plt.plot(xu,vmix,label="Vertical mix.")
plt.plot(xu,prsgrd,label="Prs. grad.")
plt.plot(xu,xadv,label="x adv.")
plt.plot(xu,yadv,label="y adv.")
plt.xlabel(r'$x$ (m)')
plt.ylabel(r'Vorticity equation terms ($1/s^2$)')
plt.title(r'Vorticity budget main term magnitudes ($<|\cdot|>$)',loc='left')
plt.legend()
plt.vlines(250-80,0,0.003,linestyle='dashed',color='gray')
plt.vlines(250,0,0.003,linestyle='dashed',color='gray')
#plt.xlim([100,270])
plt.ylim([0,0.003])
plt.savefig('./05-onlinevorticitybudget/vorticity_budget_'+cas+'.png')
plt.show()

ds = Dataset('./05-onlinevorticitybudget/'+cas+'/rip_diags_vrt.nc')
plt.figure(3,dpi=300)
title_var='vrt_drag'
var = ds.variables[title_var][-1,:,:]
plt.pcolor(var,cmap='RdBu')
plt.colorbar()
plt.clim([-0.0003,0.0003])
plt.xlabel(r'$x$ (cells)')
plt.ylabel(r'$y$ (cells)')
plt.vlines((250)/1.5,0,200,linestyle='dashed',color='gray')
plt.vlines((250-80)/1.5,0,200,linestyle='dashed',color='gray')
plt.title(title_var+r' instantaneous at last time step',loc='left')
plt.show()

plt.figure(4,dpi=300)
title_var='vrt_vmix'
var = ds.variables[title_var][-1,:,:]
plt.pcolor(var,cmap='RdBu')
plt.colorbar()
plt.clim([-0.01,0.01])
plt.xlabel(r'$x$ (cells)')
plt.ylabel(r'$y$ (cells)')
plt.vlines((250)/1.5,0,200,linestyle='dashed',color='gray')
plt.vlines((250-80)/1.5,0,200,linestyle='dashed',color='gray')
plt.title(title_var+r' instantaneous at last time step',loc='left')
plt.show()

plt.figure(5,dpi=300)
title_var='vrt_xadv'
var = ds.variables[title_var][-1,:,:]
plt.pcolor(var,cmap='RdBu')
plt.colorbar()
plt.clim([-0.1,0.1])
plt.xlabel(r'$x$ (cells)')
plt.ylabel(r'$y$ (cells)')
plt.vlines((250)/1.5,0,200,linestyle='dashed',color='gray')
plt.vlines((250-80)/1.5,0,200,linestyle='dashed',color='gray')
plt.title(title_var+r' instantaneous at last time step',loc='left')
plt.show()

plt.figure(6,dpi=300)
title_var='vrt_yadv'
var = ds.variables[title_var][-1,:,:]
plt.pcolor(var,cmap='RdBu')
plt.colorbar()
plt.clim([-0.1,0.1])
plt.xlabel(r'$x$ (cells)')
plt.ylabel(r'$y$ (cells)')
plt.vlines((250)/1.5,0,200,linestyle='dashed',color='gray')
plt.vlines((250-80)/1.5,0,200,linestyle='dashed',color='gray')
plt.title(title_var+r' instantaneous at last time step',loc='left')
plt.show()