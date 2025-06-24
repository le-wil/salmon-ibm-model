import xarray as xr

# Define the input Zarr file and output NetCDF file
zarr_file = 'output/RRtunaFADPrey_no321_npart23_nfad2_T0.50_F1.00_P1.00_I0.01_p0.0_Pa0.1.zarr'
netcdf_file = 'output/RRtunaFADPrey_no321_npart23_nfad2_T0.50_F1.00_P1.00_I0.01_p0.0_Pa0.1.nc'

# Open the Zarr file using xarray
ds = xr.open_zarr(zarr_file)

# Save the dataset to a NetCDF file
ds.to_netcdf(netcdf_file)

print(f"Converted {zarr_file} to {netcdf_file}")
