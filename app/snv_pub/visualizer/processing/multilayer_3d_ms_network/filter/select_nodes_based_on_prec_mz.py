from logging import getLogger

logger = getLogger(__name__)


def select_nodes_based_on_prec_mz(mass_lower_limit, mass_higher_limit, dic_cluster_total_input_idx_vs_cluster_info):
    # select specific node.
    logger.debug("[G2] select specific node/cluster based on prec mz")
    dic_cluster_total_input_idx_vs_cluster_info_new = {}

    for total_input_idx, cl_o in dic_cluster_total_input_idx_vs_cluster_info.items():
        if mass_lower_limit <= cl_o["represen_spec_uni"]["precursor_mz"] <= mass_higher_limit:
            dic_cluster_total_input_idx_vs_cluster_info_new[total_input_idx] = cl_o

    return dic_cluster_total_input_idx_vs_cluster_info_new
