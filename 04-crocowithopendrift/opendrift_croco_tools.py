import numpy as np
from opendrift import test_data_folder
from opendrift.readers import reader_ROMS_native, reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift
from netCDF4 import Dataset
from pyproj import Proj, Transformer
import pyproj
pyproj.datadir.set_data_dir('/opt/anaconda3/envs/opendrift2/share/proj')
import xarray as xr
import matplotlib.pyplot as plt

def add_lonlat_croco(fname):
    # === Open the NetCDF file in append mode ===
    nc = Dataset(fname, 'r+')

    # === Read Cartesian coordinates ===
    x_rho = nc.variables['x_rho'][:]
    y_rho = nc.variables['y_rho'][:]
    x_u = nc.variables['xi_u'][:]
    y_v = nc.variables['eta_v'][:]

    # === Set projection: from UTM Zone 13N to WGS84 (adjust this to match your data) ===
    transformer = Transformer.from_crs("epsg:32613", "epsg:4326", always_xy=True)

    def transform_xy(x, y):
        lon, lat = transformer.transform(x.flatten(), y.flatten())
        return np.reshape(lon, x.shape), np.reshape(lat, y.shape)

    # === Convert and write lon/lat variables ===
    for name, x, y, dim in [
        ('rho', x_rho, y_rho, ('eta_rho', 'xi_rho')),
    ]:
        lon, lat = transform_xy(x, y)

        # Create variables if they don't exist
        for varname, data in [(f'lon_{name}', lon), (f'lat_{name}', lat)]:
            if varname not in nc.variables:
                nc.createVariable(varname, 'f8', dim)
            nc.variables[varname][:] = data

    if 'angle' not in nc.variables:
        nc.createVariable('angle', datatype="f8", dimensions=('eta_rho', 'xi_rho'))
        nc['angle'].long_name = 'angle between xi axis and east'
        nc['angle'].units = 'radian'
    nc['angle'][:] = nc.variables['h'][:] * 0.0

    # === Close the NetCDF file ===
    nc.close()
    print("Longitude and latitude fields added successfully.")


def add_mask_croco(fname):
    # === Open the NetCDF file in append mode ===
    nc = Dataset(fname, 'r+')
    Dcrit = nc.variables['Dcrit'][0]
    zeta = nc.variables['zeta'][:]
    h = nc.variables['h'][:]
    # as we use averaged, just check where we have Dcrit
    mask = h>0
    masku,maskv,maskp=rho2uvp(mask)

    if 'mask_rho' not in nc.variables:
        nc.createVariable('mask_rho', 'f8', ('eta_rho', 'xi_rho'))
    nc.variables['mask_rho'][:] = mask

    if 'mask_u' not in nc.variables:
        nc.createVariable('mask_u', 'f8', ('eta_rho','xi_u'))
    nc.variables['mask_u'][:] = masku

    if 'mask_v' not in nc.variables:
        nc.createVariable('mask_v', 'f8', ('eta_v','xi_rho'))
    nc.variables['mask_v'][:] = maskv

    # === Close the NetCDF file ===
    nc.close()
    print("Masks added successfully.")


def rho2uvp(rfield):
    Mp, Lp = np.shape(rfield)  # matlab use size for "shape" in python
    M = Mp - 1
    L = Lp - 1

    vfield = 0.5 * (rfield[0:M, :] + rfield[1:Mp, :])
    ufield = 0.5 * (rfield[:, 0:L] + rfield[:, 1:Lp])
    pfield = 0.5 * (ufield[0:M, :] + ufield[1:Mp, :])

    return ufield, vfield, pfield

def lonlat_to_xy(lon,lat):
    transformer = Transformer.from_crs("epsg:4326","epsg:32613", always_xy=True)
    x,y = transformer.transform(lon,lat)
    return x,y

def xy_to_lonlat(x,y):
    transformer = Transformer.from_crs("epsg:32613", "epsg:4326", always_xy=True)
    lon,lat = transformer.transform(x, y)
    return lon,lat

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