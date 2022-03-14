import sys
import numpy as np
import validation_tools as vt

ref_dir, data_dir = vt.processArgv(sys.argv[1:])

results = np.zeros(2, dtype='i8')

# Tolerance for long run
long_tol = 11*1e8

# Tolerance for short run
short_rows = 20
short_tol = 101

# Temperature
#   Energy
results += vt.tableTest("temperature_energy.dat", ref_dir, data_dir, tol = short_tol, max_rows = short_rows)
results += vt.tableTest("temperature_energy.dat", ref_dir, data_dir, tol = long_tol)
#   L spectrum
results += vt.tableTest("temperature_l_spectrum0000.dat", ref_dir, data_dir, tol = short_tol)
results += vt.tableTest("temperature_l_spectrum0100.dat", ref_dir, data_dir, tol = long_tol, percol = True)
#   M spectrum
results += vt.tableTest("temperature_m_spectrum0000.dat", ref_dir, data_dir, tol = short_tol)
results += vt.tableTest("temperature_m_spectrum0100.dat", ref_dir, data_dir, tol = long_tol, percol = True)
#   M spectrum
results += vt.tableTest("temperature_n_spectrum0000.dat", ref_dir, data_dir, tol = short_tol)
results += vt.tableTest("temperature_n_spectrum0100.dat", ref_dir, data_dir, tol = long_tol, percol = True)

# Kinetic
#   energy
results += vt.tableTest("kinetic_energy.dat", ref_dir, data_dir, tol = short_tol, max_rows = short_rows)
#results += vt.tableTest("kinetic_energy.dat", ref_dir, data_dir, tol = long_tol)
#   L spectrum
results += vt.tableTest("kinetic_l_spectrum0000.dat", ref_dir, data_dir, tol = short_tol)
results += vt.tableTest("kinetic_l_spectrum0100.dat", ref_dir, data_dir, tol = long_tol, percol = True)
#   M spectrum
results += vt.tableTest("kinetic_m_spectrum0000.dat", ref_dir, data_dir, tol = short_tol)
results += vt.tableTest("kinetic_m_spectrum0100.dat", ref_dir, data_dir, tol = long_tol, percol = True)
#   N spectrum
results += vt.tableTest("kinetic_n_spectrum0000.dat", ref_dir, data_dir, tol = short_tol)
results += vt.tableTest("kinetic_n_spectrum0100.dat", ref_dir, data_dir, tol = long_tol, percol = True)

# Magnetic
#   energy
results += vt.tableTest("magnetic_energy.dat", ref_dir, data_dir, tol = short_tol, max_rows = short_rows)
#results += vt.tableTest("magnetic_energy.dat", ref_dir, data_dir, tol = long_tol)
#   L spectrum
results += vt.tableTest("magnetic_l_spectrum0000.dat", ref_dir, data_dir, tol = short_tol)
results += vt.tableTest("magnetic_l_spectrum0100.dat", ref_dir, data_dir, tol = long_tol, percol = True)
#   M spectrum
results += vt.tableTest("magnetic_m_spectrum0000.dat", ref_dir, data_dir, tol = short_tol)
results += vt.tableTest("magnetic_m_spectrum0100.dat", ref_dir, data_dir, tol = long_tol, percol = True)
#   N spectrum
results += vt.tableTest("magnetic_n_spectrum0000.dat", ref_dir, data_dir, tol = short_tol)
results += vt.tableTest("magnetic_n_spectrum0100.dat", ref_dir, data_dir, tol = long_tol, percol = True)

# Nusselt number
results += vt.tableTest("nusselt.dat", ref_dir, data_dir, tol = short_tol, max_rows = short_rows)
#results += vt.tableTest("nusselt.dat", ref_dir, data_dir, tol = short_tol)

# Angular momentum
#results += vt.tableTest("angular_momentum.dat", ref_dir, data_dir, tol = short_tol, max_rows = short_rows)
#results += vt.tableTest("angular_momentum.dat", ref_dir, data_dir, tol = long_tol)

# CFL
results += vt.tableTest("cfl.dat", ref_dir, data_dir, usecols=(0,1,3,5,6,7,8,9), tol = short_tol, max_rows = short_rows)
#results += vt.tableTest("cfl.dat", ref_dir, data_dir, usecols=(0,1,3,5,6,7,8,9), tol = long_tol)

# Output test summary
vt.printSummary(results)
