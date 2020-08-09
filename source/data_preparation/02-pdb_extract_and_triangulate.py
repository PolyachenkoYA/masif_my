### python -W ignore $masif_source/data_preparation/02-pdb_extract_and_triangulate.py $PDB_ID\_$CHAIN1 -p f
### python -W ignore 02-pdb_extract_and_triangulate.py $PDB1_ID $CHAIN1 $PDB2_ID $CHAIN2
# 1,2 - uX, 3,4 - C

#!/usr/bin/python
import sys
import os
import numpy as np
import glob
import subprocess
from default_config.masif_opts import masif_opts
import data_preparation.extract_and_triangulate_lib as ext_and_trg
from input_output.extractPDB import extractPDB
from input_output.save_ply import save_ply
import my_utils as my
    
args = sys.argv[1:]
argc = len(args)
if(not argc in [2, 3, 4]): 
    print('Usage:\n' + sys.argv[0] + '   pdbid1_chain1   pdbid2_chain2   [pdb_path()   [gr_trh_dist (2.0)   [to_recompute ((1)/0)]]')
    sys.exit(1)
    
ground_truth_cut_dist = float(args[2]) if len(args) > 2 else 2.0
to_recompute = (args[3] == '1') if (len(args) > 3) else True
u_pdb_name, u_chain_name, u_pdb_filepath, u_chain_filepath_base, u_chain_filepath = \
    ext_and_trg.parse_names(args[0])
#ply_filepath = os.path.join(masif_opts['ply_chain_dir'], u_pdb_name + '-dc' + str(ground_truth_cut_dist) + '_' + u_chain_name + '.ply')
ply_filepath = os.path.join(masif_opts['ply_chain_dir'], u_chain_filepath_base + '.ply')
#ply_filepath = os.path.join(ext_and_trg.pdbs_dir, u_chain_filepath_base + '.ply')
C_pdb_name, C_chain_name, C_pdb_filepath, C_chain_filepath_base, C_chain_filepath = \
    ext_and_trg.parse_names(args[1])

# Extract chains of interest.
if(to_recompute):
    extractPDB(u_pdb_filepath, u_chain_filepath, u_chain_name)

# construct unbound the mesh.
u_regular_mesh, u_vertex_normals, u_vertices, u_names = \
    ext_and_trg.msms_wrap(u_chain_filepath, to_recompute=to_recompute)
u_vertex_hbond, u_vertex_hphobicity, u_vertex_charges = \
    ext_and_trg.compute_features(u_chain_filepath_base, u_vertices, u_names, u_regular_mesh, to_recompute=to_recompute)

# construct Complex the mesh.
C_regular_mesh, C_vertex_normals, C_vertices, C_names = \
    ext_and_trg.msms_wrap(C_pdb_filepath, to_recompute=to_recompute)

# identify groundtruth
iface = ext_and_trg.find_iface(C_regular_mesh, u_regular_mesh, ground_truth_cut_dist)
    
# save results
save_ply(ply_filepath, u_regular_mesh.vertices,\
         u_regular_mesh.faces, normals=u_vertex_normals, charges=u_vertex_charges,\
         normalize_charges=True, hbond=u_vertex_hbond, hphob=u_vertex_hphobicity,\
         iface=iface)
ext_and_trg.copy_tmp2dst(ply_filepath, masif_opts['ply_chain_dir'])
#ext_and_trg.copy_tmp2dst(u_chain_filepath, masif_opts['pdb_chain_dir'])

# clean the /tmp dir
subprocess.run(r'rm /' + os.path.dirname(u_chain_filepath_base) + r'/*' + u_pdb_name + r'*', shell=True)
subprocess.run(r'rm /' + os.path.dirname(u_chain_filepath_base) + r'/*' + C_pdb_name + r'*', shell=True)