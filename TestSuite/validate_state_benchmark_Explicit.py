import sys
import numpy as np
import validation_tools as vt

ref_dir, data_dir = vt.processArgv(sys.argv[1:])

results = []

physicals = ['rayleigh', 'prandtl', 'ekman', 'magnetic', 'velocity', 'temperature']

out = vt.checkHdf5("state0000.hdf5", ref_dir, data_dir, physicals, 0)
results.append(out[0])

tids = [0, 1, 2]
tols = [1, 1, 1]
datasets = ['temperature/temperature', 'velocity/velocity_tor', 'velocity/velocity_pol', 'magnetic/magnetic_tor', 'magnetic/magnetic_pol']
for tid, t, ds in  zip(tids, tols, datasets):
    results.append(vt.hdf5Test(*out[1], ds, tid, tol = t, only_existence = True))

# Output test summary
vt.printSummary(results, tids, reftol = tols)
