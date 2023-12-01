from logging import getLogger

logger = getLogger(__name__)


def remove_nodes_based_on_total_input_idx(list_total_input_idx_to_remove, dic_cluster_total_input_idx_vs_cluster_info):
    dic_cluster_total_input_idx_vs_cluster_info_new = {}
    for cluster_total_input_idx, cluster_info in dic_cluster_total_input_idx_vs_cluster_info.items():
        if cluster_total_input_idx not in list_total_input_idx_to_remove:
            dic_cluster_total_input_idx_vs_cluster_info_new[cluster_total_input_idx] = cluster_info

    return dic_cluster_total_input_idx_vs_cluster_info_new


