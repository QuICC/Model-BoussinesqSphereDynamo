import sys
import numpy as np
import validation_tools as vt

ref_dir, data_dir = vt.processArgv(sys.argv[1:])

results = []

# Tolerance per max rows
rows = list(range(0, 101, 10))
tols = [51, 61, 81, 141, 201, 261, 301, 331, 451, 491, 501] # Without n spectra
tols = [61, 81, 111, 191, 221, 271, 311, 381, 461, 541, 571] # With n spectra

prefixes = ['temperature', 'kinetic', 'magnetic']
spectra = ['l', 'm', 'n']

# check energy and spectra
for prefix in prefixes:
    # Energy
    for r, t in zip(rows,tols):
        results.append(vt.tableTest(prefix + '_energy.dat', ref_dir, data_dir, r, tol = t, max_rows = r+1))

    # Spectra
    for mode in spectra:
        for r, t in zip(rows,tols):
            results.append(vt.tableTest(prefix +  '_' + mode + f'_spectrum{r:04}.dat', ref_dir, data_dir, r, tol = t, percol = True))

# Nusselt number
for r, t in zip(rows,tols):
    results.append(vt.tableTest("nusselt.dat", ref_dir, data_dir, r, tol = t, max_rows = r+1))

# CFL
for r, t in zip(rows,tols):
    results.append(vt.tableTest("cfl.dat", ref_dir, data_dir, r, usecols=(0,1,3,5,7,8,9,10,11), tol = t, max_rows = r+1))

# Angular momentum
#for r, t in zip(rows,tols):
#    results.append(vt.tableTest("angular_momentum.dat", ref_dir, data_dir, r, tol = t, max_rows = r+1))

# Output test summary
vt.printSummary(results, rows)