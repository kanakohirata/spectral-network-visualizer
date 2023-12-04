from logging import getLogger


logger = getLogger(__name__)


def make_edges_and_nodes_inter_layer(list_of_edge_for_networkx: list,
                                     dic_cluster_total_input_idx_MOD_vs_node_info: dict,
                                     fo_log_key) -> tuple:
    """

    Parameters
    ----------
    list_of_edge_for_networkx : list
    dic_cluster_total_input_idx_MOD_vs_node_info : list
    fo_log_key : file object

    Returns
    -------
    list_of_edge_for_networkx_to_show_inter_sample_ref_layer : list
    list_of_edge_for_networkx_to_show_inter_sample_layer : list
    dic_cluster_total_input_idx_MOD_vs_node_info_inter_sample_ref_layer : dict
    dic_cluster_total_input_idx_MOD_vs_node_info_inter_sample_layer : dict
    list_total_input_idx_MOD_inter_sample_ref_layer : list
    list_total_input_idx_MOD_inter_sample_layer : list
    """
    # [EDGES]
    # create list_od_edge for NetworkX for particular layer. -----------------------------------
    list_of_edge_for_networkx_to_show_inter_sample_ref_layer = []
    list_of_edge_for_networkx_to_show_inter_sample_layer = []
    # this list keeps total_input_idx_MOD for "inter sample ref layer"
    list_total_input_idx_MOD_inter_sample_ref_layer = []
    list_total_input_idx_MOD_inter_sample_layer = []

    logger.debug(f" len of list_of_edge_for_networkx_to_show {len(list_of_edge_for_networkx)}")
    for edge_for_networkx in list_of_edge_for_networkx:
        # we are only looking for inner layer edge

        if edge_for_networkx[2]["edge_type"] == "inter_sample_ref_layer":
            # if the edge is ""inter sample ref edge",
            list_of_edge_for_networkx_to_show_inter_sample_ref_layer.append(edge_for_networkx)
            list_total_input_idx_MOD_inter_sample_ref_layer.append(edge_for_networkx[0])
            list_total_input_idx_MOD_inter_sample_ref_layer.append(edge_for_networkx[1])

        if edge_for_networkx[2]["edge_type"] == "inter_sample_layer":
            list_of_edge_for_networkx_to_show_inter_sample_layer.append(edge_for_networkx)
            list_total_input_idx_MOD_inter_sample_layer.append(edge_for_networkx[0])
            list_total_input_idx_MOD_inter_sample_layer.append(edge_for_networkx[1])

    fo_log_key.write("\n list_of_edge_for_networkx_to_show_inter_sample_ref_layer : "
                     f"{len(list_of_edge_for_networkx_to_show_inter_sample_ref_layer)}")
    fo_log_key.write("\n list_of_edge_for_networkx_to_show_inter_sample_layer : "
                     f"{len(list_of_edge_for_networkx_to_show_inter_sample_layer)}")
    fo_log_key.flush()
    logger.info(f'Length of list_of_edge_for_networkx_to_show_inter_sample_ref_layer: '
                f'{len(list_of_edge_for_networkx_to_show_inter_sample_ref_layer)}')
    logger.info(f'Length of ist_of_edge_for_networkx_to_show_inter_sample_layer: '
                f'{len(list_of_edge_for_networkx_to_show_inter_sample_layer)}')

    # make it non redundant.
    list_total_input_idx_MOD_inter_sample_ref_layer = list(set(list_total_input_idx_MOD_inter_sample_ref_layer))
    list_total_input_idx_MOD_inter_sample_layer = list(set(list_total_input_idx_MOD_inter_sample_layer))

    # [NODES]
    #  create dic  _cluster_total_input_idx_MOD_vs_node_info for each layer ------------------------
    dic_cluster_total_input_idx_MOD_vs_node_info_inter_sample_ref_layer = {}
    dic_cluster_total_input_idx_MOD_vs_node_info_inter_sample_layer = {}

    for cluster_total_input_idx_MOD, node_info in dic_cluster_total_input_idx_MOD_vs_node_info.items():
        if cluster_total_input_idx_MOD in list_total_input_idx_MOD_inter_sample_ref_layer:
            dic_cluster_total_input_idx_MOD_vs_node_info_inter_sample_ref_layer[cluster_total_input_idx_MOD] = node_info

        if cluster_total_input_idx_MOD in list_total_input_idx_MOD_inter_sample_layer:
            dic_cluster_total_input_idx_MOD_vs_node_info_inter_sample_layer[cluster_total_input_idx_MOD] = node_info

    return (list_of_edge_for_networkx_to_show_inter_sample_ref_layer,
            list_of_edge_for_networkx_to_show_inter_sample_layer,
            dic_cluster_total_input_idx_MOD_vs_node_info_inter_sample_ref_layer,
            dic_cluster_total_input_idx_MOD_vs_node_info_inter_sample_layer,
            list_total_input_idx_MOD_inter_sample_ref_layer,
            list_total_input_idx_MOD_inter_sample_layer)
