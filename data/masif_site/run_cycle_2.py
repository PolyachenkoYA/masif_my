import os
import sys
import subprocess as sp
import numpy as np

# ================= cycle params ==================
pdbs = ['1D6R', '1JTG', '1R6Q', '1RKE', '1ZHH', '2HQS', '2I25', '2UUY', '3SGQ']
steps = [50, 1000]
mds = ['uR', 'uL']

steps = [1000]

# ============== input parse ==============
args = sys.argv[1:]
argc = len(args)

# ================ cycle =====================
#for pdb_id in pdbs:
for arg in args:
    for md_id in mds:
        for step in steps:
            sp.run(['python', 'run_cycle.py', pdbs[int(arg)], md_id, str(step)])
            