import xarray as xr

# If you originally had a Zarr store or some other object, open it appropriately:
# ds = xr.open_zarr('collated_all/quantiles/rslpismflat_icemod_quants')
# or, if you can open the same object you used before:
ds = xr.open_dataset('collated_all/quantiles/rslpismflat_icemod_quants.nc', engine='netcdf4')

# Write a genuine NetCDF4 file
ds.to_netcdf('rslpismflat_icemod_quants_REAL.nc',
             engine='netcdf4',
             format='NETCDF4_CLASSIC')
