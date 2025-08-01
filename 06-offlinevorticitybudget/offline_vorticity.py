import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

# === Load CROCO output ===
ds = xr.open_dataset("./06-offlinevorticitybudget/FLASH_RIP/rip_his.nc")

# === Grid size (or use 1/pm, 1/pn if available) ===
dx = 1.5  # m
dy = 1.5  # m

# === Load variables ===
u = ds['u'][:]
v = ds['v']
w = ds['w']
zeta_surf = ds['zeta']
h = ds['h']
time = ds['zeta'].coords['time']
eta_rho = ds['zeta'].coords['eta_rho']
xi_rho = ds['zeta'].coords['xi_rho']

# Drop coords to avoid alignment errors
u = u.reset_coords(drop=True)
v = v.reset_coords(drop=True)

# Pad with edge values (or constant if you prefer NaNs)
u_pad = u.pad(xi_u=(0, 2), mode='edge')  # or mode='constant', constant_values=np.nan
v_pad = v.pad(eta_v=(0, 2), mode='edge')

# Interpolate to rho-points
u_rho = (u_pad.data[:,:,:,:-1]+u_pad.data[:,:,:,1:])*0.5
v_rho = (v_pad.data[:,:,:-1,:]+v_pad.data[:,:,1:,:])*0.5
u_rho = xr.DataArray(u_rho, coords=[time, ds['s_rho'], eta_rho, xi_rho], dims=['time', 's_rho', 'eta_rho', 'xi_rho'])
v_rho = xr.DataArray(v_rho, coords=[time, ds['s_rho'], eta_rho, xi_rho], dims=['time', 's_rho', 'eta_rho', 'xi_rho'])

# === Depth-averaged u, v ===
ubar = u_rho.mean(dim="s_rho")
vbar = v_rho.mean(dim="s_rho")

# === Horizontal gradients of bar{u}, bar{v} ===
dubar_dx = ubar.diff("xi_rho") / dx
dubar_dx = dubar_dx.pad(xi_rho=(0, 1), mode="edge")

dvbar_dy = vbar.diff("eta_rho") / dy
dvbar_dy = dvbar_dy.pad(eta_rho=(0, 1), mode="edge")

# === Horizontal derivatives for vorticity ===
dv_dx = vbar.diff("xi_rho") / dx
dv_dx = dv_dx.pad(xi_rho=(1, 0), mode="edge")

du_dy = ubar.diff("eta_rho") / dy
du_dy = du_dy.pad(eta_rho=(0, 1), mode="edge")

zeta_bar = dv_dx.data - du_dy.data  # depth-averaged vertical vorticity
zeta_bar = xr.DataArray(zeta_bar, coords=[time, eta_rho, xi_rho], dims=['time', 'eta_rho', 'xi_rho'])


# === ∂ζ̄/∂t ===
dzeta_dt = zeta_bar.diff("time") / zeta_bar["time"].diff("time")
dzeta_dt = dzeta_dt.pad(time=(0, 1), mode="edge")

# === Advection term: u ∂ζ/∂x + v ∂ζ/∂y ===
dzeta_dx = zeta_bar.diff("xi_rho") / dx
dzeta_dx = dzeta_dx.pad(xi_rho=(0, 1), mode="edge")

dzeta_dy = zeta_bar.diff("eta_rho") / dy
dzeta_dy = dzeta_dy.pad(eta_rho=(0, 1), mode="edge")

adv_term = ubar * dzeta_dx.data + vbar * dzeta_dy.data
adv_term = xr.DataArray(adv_term, coords=[time, eta_rho, xi_rho], dims=['time', 'eta_rho', 'xi_rho'])

# === Compression term: ζ̄ ∇·ū ===
div_uv = dubar_dx.data + dvbar_dy.data
div_uv = xr.DataArray(div_uv, coords=[time, eta_rho, xi_rho], dims=['time', 'eta_rho', 'xi_rho'])
comp_term = zeta_bar * div_uv

# === Vortex stretching term: (ζw|surf − ζw|bot)/H ===
# ζ3D for surface and bottom
dvdx = (v_rho.diff("xi_rho") / dx).pad(xi_rho=(0, 1), mode="edge")
dudy = (u_rho.diff("eta_rho") / dy).pad(eta_rho=(0, 1), mode="edge")
zeta3D = dvdx.data - dudy.data
zeta3D = xr.DataArray(zeta3D, coords=[time, ds['s_rho'],eta_rho, xi_rho], dims=['time', 's_rho','eta_rho', 'xi_rho'])

zeta_top = zeta3D.isel(s_rho=-1)
zeta_bot = zeta3D.isel(s_rho=0)
w_top = w.isel(s_w=-1)
w_bot = w.isel(s_w=0)

H = h + zeta_surf
H = H.where(H > 0.1)  # prevent divide by 0

vstretch = (zeta_top * w_top - zeta_bot * w_bot) / H

# === Residual (should be small) ===
# All terms aligned in time
dzeta_dt = dzeta_dt.isel(time=884)  # choose a well-defined time slice
adv_term = adv_term.isel(time=884)
comp_term = comp_term.isel(time=884)
vstretch = vstretch.isel(time=884)

residual = dzeta_dt - (-adv_term - comp_term + vstretch)

# === Plot closure ===
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.title("Vorticity Budget Residual")
plt.pcolormesh(residual, cmap='seismic', shading='auto')
plt.colorbar(label='Residual (1/s²)')

plt.subplot(1, 2, 2)
plt.title("Histogram of Residual")
plt.hist(residual.values.flatten(), bins=100)
plt.tight_layout()
plt.show()


#%%

fig, axs = plt.subplots(2, 3, figsize=(15, 10), constrained_layout=True)
terms = {
    "∂ζ/∂t": dzeta_dt,
    "Advection": adv_term,
    "Compression": comp_term,
    "Vortex stretching": vstretch,
    "Residual": residual
}
vmin = -1  # adjust based on typical magnitude
vmax = 1

for ax, (label, data) in zip(axs.flat, terms.items()):
    im = ax.pcolormesh(data, cmap='seismic', vmin=vmin, vmax=vmax, shading='auto')
    ax.set_title(label)
    plt.colorbar(im, ax=ax, shrink=0.8)

# Remove last empty subplot if 5 terms
for i in range(len(terms), 6):
    fig.delaxes(axs.flat[i])

plt.suptitle("Vorticity Budget Terms at t = {:.1f} s".format(ds.time.values[1].astype(float)), fontsize=16)
plt.show()

#%% Tools
def u_to_rho(var):
    var_pad = var.pad(xi_u=(0, 2), mode='edge')
    var_rho = (var_pad.isel(xi_u=slice(0, -1)).data + var_pad.isel(xi_u=slice(1,None)).data) * 0.5
    var_rho = xr.DataArray(var_rho, coords=[time, ds['s_rho'], eta_rho, xi_rho],
                         dims=['time', 's_rho', 'eta_rho', 'xi_rho'])
    return var_rho

def v_to_rho(var):
    var_pad = var.pad(eta_v=(0, 2), mode='edge')
    var_rho = (var_pad.isel(eta_v=slice(0, -1)).data + var_pad.isel(eta_v=slice(1,None)).data) * 0.5
    var_rho = xr.DataArray(var_rho, coords=[time, ds['s_rho'], eta_rho, xi_rho],
                         dims=['time', 's_rho', 'eta_rho', 'xi_rho'])
    return var_rho
