#! /usr/bin/env python
#
# This script compares RFMIP results from RTE+RRTMGP against a benchmark
#
import os
import numpy as np
import xarray as xr

ref_dir  = "ref"
tst_dir = "."

#
# Comparing reference and test results
#
tst = xr.open_dataset(os.path.join(tst_dir, "rrtmgp-clouds.nc"))
ref = xr.open_dataset(os.path.join(ref_dir, "rrtmgp-clouds.nc"))

for v in ['lw_flux_up', 'lw_flux_dn', 'sw_flux_up', 'sw_flux_dn', 'sw_flux_dir']:
  if np.all(np.isnan(tst.variables[v].values)):
    raise Exception("All test values are missing. Were the tests run?")
  if np.any(np.isnan(tst.variables[v].values)):
    raise Exception("Some test values are missing. Now that is strange.")

  diff = abs((tst-ref).variables[v].values)
  avg  = 0.5*(tst+ref).variables[v].values
  # Division raises a runtime warning when we divide by zero even if the
  #   values in those locations will be ignored.
  with np.errstate(divide='ignore', invalid='ignore'):
    frac_diff = np.where((avg > 2.*np.finfo(float).eps), diff/avg, 0)

  if diff.max() > 0:
    print('Variable %s differs (max abs difference: %e; max percent difference: %e%%)'%(v, diff.max(), 100.0 * frac_diff.max()))
  else:
    print('Variable %s: No diffs'%(v))
