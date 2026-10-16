import sys,importlib.util
print(sys.version)
for m in ['gsw','xarray','netCDF4','h5netcdf','scipy']:
 print(m, bool(importlib.util.find_spec(m)))
