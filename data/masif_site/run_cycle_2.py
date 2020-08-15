import os
import sys
import subprocess as sp
import numpy as np

# ================= cycle params ==================
pdbs = ['1D6R', '1JTG', '1R6Q', '1RKE', '1Z0K', '1ZHH', '2HQS', '2I25', '2O3B', '2UUY', '3SGQ']
pdbs = ['1AK4', '1CGI', '1CLV', '1D6R', '1E96', '1FLE', '1GPW', '1JTG', '1R0R', '1R6Q', '1RKE', '1ZHH', '2HQS', '2I25', '2O3B', '2OOB', '2UUY', '3F1P', '3SGQ']
step = 1000
mds = ['uR', 'uL']

# ============== input parse ==============
args = sys.argv[1:]
argc = len(args)
pdbs_to_do = [int(i) for i in args]

# ================ cycle =====================
#for pdb_id in pdbs:
for pdb_i in pdbs_to_do:
    for md_id in mds:
        sp.run(['python', 'run_cycle.py', pdbs[pdb_i], md_id, str(step)])
            