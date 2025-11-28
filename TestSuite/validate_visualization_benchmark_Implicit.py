import sys
import numpy as np
import validation_tools as vt

ref_dir, data_dir = vt.processArgv(sys.argv[1:])

results = []

physicals = ['rayleigh', 'prandtl', 'ekman', 'magnetic', 'velocity', 'temperature']

out = vt.checkHdf5("visState0000.hdf5", ref_dir, data_dir, physicals, 0)
results.append(out[0])

tids = [0, 1, 2, 3, 4, 5, 6]
tols = [50, 50, 50, 50, 50, 50, 50]
datasets = ['temperature/temperature', 'velocity/velocity_r', 'velocity/velocity_theta', 'velocity/velocity_phi', 'magnetic/magnetic_r', 'magnetic/magnetic_theta', 'magnetic/magnetic_phi']
for tid, t, ds in  zip(tids, tols, datasets):
    results.append(vt.hdf5Test(*out[1], ds, tid, tol = t, peraxis = (0,1,2), threshold = 2e-12))

# Output test summary
vt.printSummary(results, tids, reftol = tols)
