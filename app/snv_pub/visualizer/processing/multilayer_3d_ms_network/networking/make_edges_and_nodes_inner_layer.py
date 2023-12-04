from logging import getLogger


logger = getLogger(__name__)


def make_edges_and_nodes_inner_layer(list_attribute_for_layer: list,
                                     list_of_edge_for_networkx: list,
                                     dic_cluster_total_input_idx_MOD_vs_node_info: dict,
                                     fo_log) -> list:
    """

    Parameters
    ----------
    list_attribute_for_layer
    list_of_edge_for_networkx
    dic_cluster_total_input_idx_MOD_vs_node_info
    fo_log

    Returns
    -------
    list_dic_edges_nodes_graph_by_layer : list
        [
            {'attribute_for_layer': sample_aaa.msp, 'list_of_edge_for_networkx': [...], 'dic_cluster_total_input_idx_MOD_vs_node_info': {...}},
            {'attribute_for_layer': Benzenoids, 'list_of_edge_for_networkx': [...], 'dic_cluster_total_input_idx_MOD_vs_node_info': {...}}, ...
        ]
    """
    # making dataset for inner layer edges, node.
    list_dic_edges_nodes_graph_by_layer = []
    # for each layer. (inner edges, and related nodes...)
    for attribute_for_layer in list_attribute_for_layer:
        dic_edges_nodes_by_layer = {"attribute_for_layer": attribute_for_layer}
        fo_log.write(f"\n-------- [attribute_for_layer] : {attribute_for_layer}")

        # [EDGES]
        # create list_od_edge for NetworkX for particular layer. -----------------------------------
        list_of_edge_for_networkx_to_show_by_layer = []

        for edge_for_networkx in list_of_edge_for_networkx:

            # we are only looking for inner layer edge
            if edge_for_networkx[2]["edge_type"] in ("inter_sample_ref_layer", "inter_ref_layer", "inter_sample_layer"):
                continue

            # if both node is on the layer,
            if attribute_for_layer.startswith("sample"):
                fo_log.write(
                    f" \n attribute match 1:"
                    f"---{dic_cluster_total_input_idx_MOD_vs_node_info[edge_for_networkx[0]]['layer']}---"
                    f" VS ---{attribute_for_layer}---")
                fo_log.write(
                    f" \n attribute match 2:"
                    f"---{dic_cluster_total_input_idx_MOD_vs_node_info[edge_for_networkx[1]]['layer']}---"
                    f" VS ---{attribute_for_layer}---")

            if (dic_cluster_total_input_idx_MOD_vs_node_info[edge_for_networkx[0]]["layer"] == attribute_for_layer
                    and dic_cluster_total_input_idx_MOD_vs_node_info[edge_for_networkx[1]]["layer"] == attribute_for_layer):
                fo_log.write(" \n MATCHED")
                list_of_edge_for_networkx_to_show_by_layer.append(edge_for_networkx)
                if attribute_for_layer.startswith("sample"):
                    logger.debug("sample edge appended..........................")

        dic_edges_nodes_by_layer["list_of_edge_for_networkx"] = list_of_edge_for_networkx_to_show_by_layer

        fo_log.write("\n  for this layer:   dic_edges_nodes_by_layer[list_of_edge_for_networkx] "
                     f"{len(list_of_edge_for_networkx_to_show_by_layer)}")
        logger.info(f"For '{attribute_for_layer}' layer, length of list_of_edge_for_networkx: "
                    f"{len(list_of_edge_for_networkx_to_show_by_layer)}")

        # [NODES]
        # create dic _cluster_total_input_idx_MOD_vs_node_info for each layer ------------------------
        dic_cluster_total_input_idx_MOD_vs_node_info_by_layer = {}
        for cluster_total_input_idx_MOD, node_info in dic_cluster_total_input_idx_MOD_vs_node_info.items():
            if node_info["layer"] == attribute_for_layer:
                dic_cluster_total_input_idx_MOD_vs_node_info_by_layer[cluster_total_input_idx_MOD] = node_info

        dic_edges_nodes_by_layer["dic_cluster_total_input_idx_MOD_vs_node_info"] =\
            dic_cluster_total_input_idx_MOD_vs_node_info_by_layer
        list_dic_edges_nodes_graph_by_layer.append(dic_edges_nodes_by_layer)

    return list_dic_edges_nodes_graph_by_layer
