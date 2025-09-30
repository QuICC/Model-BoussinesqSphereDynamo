import sys
import numpy as np
import validation_tools as vt

ref_dir, data_dir = vt.processArgv(sys.argv[1:])

results = []

tids = [0]
tols = [1]

for tid, t in  zip(tids, tols):
    results.append(vt.checkXml("parameters_template_WLFm.cfg", ref_dir, data_dir, tid))

# Output test summary
vt.printSummary(results, tids, reftol = tols)
