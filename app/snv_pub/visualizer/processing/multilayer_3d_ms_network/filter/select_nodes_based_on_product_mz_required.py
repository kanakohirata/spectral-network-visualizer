from logging import getLogger
import numpy as np

logger = getLogger(__name__)


def select_nodes_based_on_product_mz_required(list_product_mz_required, mz_tol, dic_cluster_total_input_idx_vs_cluster_info):
    logger.debug('starting select_nodes_based_on_product_mz_required')
    dic_cluster_total_input_idx_vs_cluster_info_new = {}
    for total_input_idx, cl_o in dic_cluster_total_input_idx_vs_cluster_info.items():

        # for ref spec data, you can just mass
        flag_add_to_new_dic = 1

        # you have to examine product ion mz only for sample data.
        if cl_o['tag'].startswith('sample'):
            flag_add_to_new_dic = 0

            # for peaks in represen spec cluster -------------------------
            fragment_mz_list = []
            for pk in cl_o['represen_spec_uni']['peak_list_mz_int_rel']:
                fragment_mz_list.append(pk[0])

            for product_mz_required in list_product_mz_required:
                is_mz_delta_less_than_mz_tol = list(map(lambda x: abs(x - product_mz_required) < mz_tol, fragment_mz_list))
                if np.any(is_mz_delta_less_than_mz_tol):
                    flag_add_to_new_dic = 1
                    break

        # if mass defect match, add to new dictionary
        if flag_add_to_new_dic == 1:
            dic_cluster_total_input_idx_vs_cluster_info_new[total_input_idx] = cl_o

    logger.debug('finishing select_nodes_based_on_product_mz_required')
    return dic_cluster_total_input_idx_vs_cluster_info_new
