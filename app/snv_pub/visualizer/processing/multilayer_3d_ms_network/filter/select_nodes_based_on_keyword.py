from logging import getLogger

logger = getLogger(__name__)


def select_nodes_based_on_keyword(filter_select_category, filter_select_keyword, dic_cluster_total_input_idx_vs_cluster_info):
    # select specific node .
    logger.debug('[G2] select specific node/cluster')
    logger.debug(f'dic_config["filter_select_category"]: {filter_select_category}')

    if filter_select_category == 'list_cmpd_classification_superclass':
        dic_cluster_total_input_idx_vs_cluster_info_new = {}

        for total_input_idx, cl_info in dic_cluster_total_input_idx_vs_cluster_info.items():
            # if the current cluster is sample OR to be selected.
            if (cl_info['tag'] == 'sample'
                    or filter_select_keyword in cl_info['represen_spec_uni']['list_cmpd_classification_superclass']):
                dic_cluster_total_input_idx_vs_cluster_info_new[total_input_idx] = cl_info

        return dic_cluster_total_input_idx_vs_cluster_info_new

    else:
        return dic_cluster_total_input_idx_vs_cluster_info
