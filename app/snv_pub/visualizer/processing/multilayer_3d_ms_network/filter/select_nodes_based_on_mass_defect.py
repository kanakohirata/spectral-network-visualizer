from logging import getLogger

logger = getLogger(__name__)


def select_nodes_based_on_mass_defect(mz_tol, list_mass_defect, dic_cluster_total_input_idx_vs_cluster_info):
    # mz tolerance to examine mass defect
    dic_cluster_total_input_idx_vs_cluster_info_new = {}
    count_mass_defect_match = 0
    count_sample = 0
    count_dic = 0

    logger.debug(f'list_mass_defect {list_mass_defect}')

    for total_input_idx, cl_o in dic_cluster_total_input_idx_vs_cluster_info.items():
        count_dic += 1
        # for ref spec data, you can just mass
        flag_add_to_new_dic = 1

        # for sample data, you have to examine mass defect
        if cl_o['tag'].startswith('sample'):

            flag_add_to_new_dic = 0
            count_sample += 1

            # get mass defect for current node.
            node_mass_defect = cl_o['represen_spec_uni']['precursor_mz'] - int(cl_o['represen_spec_uni']['precursor_mz'])

            # for all mass defect value specified in config.
            for mass_df_specified in list_mass_defect:
                # examine mass defect match to config val
                if abs(node_mass_defect - mass_df_specified) < mz_tol:
                    _prec_mz = cl_o['represen_spec_uni']['precursor_mz']
                    count_mass_defect_match += 1
                    flag_add_to_new_dic = 1
                    break

        # if mass defect match, add to new dictionary
        if flag_add_to_new_dic == 1:
            dic_cluster_total_input_idx_vs_cluster_info_new[total_input_idx] = cl_o
    logger.debug(f"count dic {str(count_dic)}")
    logger.debug(f"count sample {str(count_sample)}")
    logger.debug(f" matched mass defect: {str(count_mass_defect_match)}")
    logger.debug(f"len of dic after mass defect exam: {len(dic_cluster_total_input_idx_vs_cluster_info_new)}")

    return dic_cluster_total_input_idx_vs_cluster_info_new
