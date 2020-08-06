import pymesh
import os
import sys
from IPython.core.debugger import set_trace
from sklearn.metrics import roc_auc_score
import importlib
import numpy as np
from default_config.masif_opts import masif_opts

import data_preparation.extract_and_triangulate_lib as ext_and_trg

"""
masif_site_label_surface.py: Color a protein ply surface file by the MaSIF-site interface score.
Pablo Gainza - LPDI STI EPFL 2019
This file is part of MaSIF.
Released under an Apache License 2.0
"""

args = sys.argv[1:]
argc = len(args)
if(not argc in [3]):
    print('usage:\n' + sys.argv[0] + '   params_model   PPI_id   ROCAUC_filename')
custom_params_file = args[0]
ppi_pair_id = args[1]
ROCAUC_filename = args[2]

params = masif_opts["site"]
custom_params = importlib.import_module(custom_params_file, package=None)
custom_params = custom_params.custom_params
ROCAUC_filepath = os.path.join(ext_and_trg.pdbs_dir, ROCAUC_filename)

for key in custom_params:
    print("Setting {} to {} ".format(key, custom_params[key]))
    params[key] = custom_params[key]

# Shape precomputation dir.
parent_in_dir = params["masif_precomputation_dir"]

[pdb_id, chain_name] = ppi_pair_id.split("_")
frame_id = pdb_id.split('-')[-1]
scrores_filepath = params["out_pred_dir"] + "/pred_" + ppi_pair_id + ".npy"
ply_file = masif_opts["ply_file_template"].format(pdb_id, chain_name)
    
try:
    mymesh = pymesh.load_mesh(ply_file)
except:
    print("File does not exist: {}".format(ply_file))
    exit(1)

try:
    scores = np.load(scrores_filepath)
except:
    print("File does not exist: {}".format(scrores_filepath))
    exit(1)

ground_truth = mymesh.get_attribute('vertex_iface')
# Compute ROC AUC for this protein. 
try:
    roc_auc = roc_auc_score(ground_truth, scores[0])
    print("ROC AUC score for protein {} : {:.2f} ".format(ppi_pair_id, roc_auc))
    
    with open(ROCAUC_filepath, 'a') as rocauc_file:
        print(frame_id, roc_auc, file=rocauc_file)
        print('ROCAUC = ' + str(roc_auc))    
except: 
    print("No ROC AUC computed for protein (possibly, no ground truth defined in input)") 
