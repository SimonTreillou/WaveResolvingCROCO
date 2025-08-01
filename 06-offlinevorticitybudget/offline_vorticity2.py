import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# === Load CROCO file ===
ds = xr.open_dataset("./06-offlinevorticitybudget/FLASH_RIP/rip_his.nc")

def extract_last_timestep(input_file, output_file):
    """
    Extracts the last time step from a CROCO NetCDF file and writes it to a new file.

    Parameters:
        input_file (str): Path to the input NetCDF file
        output_file (str): Path for the output NetCDF file with the last time step
    """
    # Open the dataset with decoding of time
    ds = xr.open_dataset(input_file)

    # Detect time dimension (often 'ocean_time' in CROCO)
    time_dim = [dim for dim in ds.dims if 'time' in dim]
    if not time_dim:
        raise ValueError("No time dimension found in dataset.")
    time_dim = time_dim[0]

    # Select last time step
    ds_last = ds.isel({time_dim: -1})

    # Keep time dimension length 1 (not squeezed)
    ds_last = ds_last.expand_dims({time_dim: 1})

    # Save to output file
    ds_last.to_netcdf(output_file)
    print(f"Last time step saved to {output_file}")

extract_last_timestep("./06-offlinevorticitybudget/FLASH_RIP/rip_his.nc", "./06-offlinevorticitybudget/FLASH_RIP/rip_last.nc")

# === Extract variables ===
u = ds['u']             # (time, s_rho, eta_u, xi_u)
v = ds['v']             # (time, s_rho, eta_v, xi_v)
ubar = ds['ubar']       # (time, eta_u, xi_u)
vbar = ds['vbar']       # (time, eta_v, xi_v)
zeta = ds['zeta']       # Free surface elevation
h = ds['h']             # Bathymetry
dx = ds['pm']**-1       # Grid spacing in x
dy = ds['pn']**-1       # Grid spacing in y

# === Time step (assumed uniform) ===
dt = float(ds['scrum_time'][1] - ds['scrum_time'][0])  # 1s typically

# === Vorticity computation (depth-averaged) ===
def compute_depth_averaged_vorticity(ubar, vbar, dx, dy):
    dvdx = (vbar.shift(xi_rho=-1) - vbar) / dx
    dudy = (ubar.shift(eta_rho=-1) - ubar) / dy
    dvdx = dvdx.interp_like(zeta)  # onto rho grid
    dudy = dudy.interp_like(zeta)
    return dvdx - dudy

omega = compute_depth_averaged_vorticity(ubar, vbar, dx, dy)

# === Local rate of change ===
domega_dt = (omega.shift(ocean_time=-1) - omega) / dt
domega_dt = domega_dt.isel(ocean_time=slice(0, -1))

# === Advection term ===
# Interpolate ubar, vbar to rho grid
ubar_rho = 0.5 * (ubar[:, :, :-1] + ubar[:, :, 1:])
vbar_rho = 0.5 * (vbar[:, :-1, :] + vbar[:, 1:, :])
ubar_rho = ubar_rho.rename({'xi_u': 'xi_rho', 'eta_u': 'eta_rho'})
vbar_rho = vbar_rho.rename({'xi_v': 'xi_rho', 'eta_v': 'eta_rho'})

domegadx = (omega.shift(xi_rho=-1) - omega) / dx
domegady = (omega.shift(eta_rho=-1) - omega) / dy

advection = ubar_rho[:, :-1, :-1] * domegadx + vbar_rho[:, :-1, :-1] * domegady

# === Stretching term: -ω (∂u/∂x + ∂v/∂y)
dudx = (ubar.shift(xi_u=-1) - ubar) / dx
dvdy = (vbar.shift(eta_v=-1) - vbar) / dy
dudx_rho = dudx.interp_like(zeta)
dvdy_rho = dvdy.interp_like(zeta)
divergence = dudx_rho + dvdy_rho
stretching = -omega[:, :-1, :-1] * divergence

# === Plotting (mean over time)
plt.figure(figsize=(12, 6))
plt.subplot(1, 3, 1)
plt.title('∂ω/∂t')
plt.pcolormesh(domega_dt.mean(dim='ocean_time'), shading='auto')
plt.colorbar()

plt.subplot(1, 3, 2)
plt.title('Advection')
plt.pcolormesh(advection.mean(dim='ocean_time'), shading='auto')
plt.colorbar()

plt.subplot(1, 3, 3)
plt.title('Stretching')
plt.pcolormesh(stretching.mean(dim='ocean_time'), shading='auto')
plt.colorbar()
plt.tight_layout()
plt.show()
