from copy import deepcopy
from logging import getLogger
import os
import pandas as pd

logger = getLogger(__name__)


def add_suspect_compound_info(dic_filename_vs_dic_suspect_cmpd_vs_mz,
                              mz_tol,
                              dic_cluster_total_input_idx_MOD_vs_node_info):

    dic_cluster_total_input_idx_MOD_vs_node_info_new = deepcopy(dic_cluster_total_input_idx_MOD_vs_node_info)
    for cluster_total_input_idx, node_info in dic_cluster_total_input_idx_MOD_vs_node_info_new.items():

        # suspect screening only applies to sample layer---
        if node_info['layer'].startswith('sample'):

            # !!!!!!!!!!!! TO BE FIXED
            # currently only mz of representative spec of cluster is matched to suspect list.
            # for each spectrum in the cluster....
            # for sp_uni in node_info.spec_cluster.list_spec_uni:
            # scan all suspect mz ---------------

            # for each suspect file,
            for path, dic_cmpd_vs_mz in dic_filename_vs_dic_suspect_cmpd_vs_mz.items():
                filename = os.path.splitext(os.path.basename(path))[0]
                # for each compound in suspect
                for cmpd, mz_suspect in dic_cmpd_vs_mz.items():
                    # if suspect mz match to observed mass.
                    if abs(node_info['spec_cluster']['represen_spec_uni']['precursor_mz'] - mz_suspect) < mz_tol:
                        node_info['suspect'] = filename
                        if filename not in node_info['l_source_suspect']:
                            node_info['l_source_suspect'].append(filename)

                        if not node_info['name']:
                            node_info['name'] = f'{cmpd}_in({filename})'
                        else:
                            node_info['name'] += f'; {cmpd}_in({filename})'

    return dic_cluster_total_input_idx_MOD_vs_node_info_new
