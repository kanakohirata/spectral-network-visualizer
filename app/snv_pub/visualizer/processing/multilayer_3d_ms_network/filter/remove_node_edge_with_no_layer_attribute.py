from logging import getLogger

logger = getLogger(__name__)


def remove_node_edge_with_no_layer_attribute(dic_cluster_total_input_idx_MOD_vs_node_info,
                                             list_of_edge_for_networkx):
    """
    Remove nodes and edges with layer attribute is "none".
    Parameters
    ----------
    dic_cluster_total_input_idx_MOD_vs_node_info
    list_of_edge_for_networkx

    Returns
    -------
    dic_cluster_total_input_idx_MOD_vs_node_info_new : dict
        Updated dic_cluster_total_input_idx_MOD_vs_node_info
    list_of_edge_for_networkx_new : dict
        Updated list_of_edge_for_networkx
    """
    #####################
    #  [R] removing node/edge with no-layer attribute.
    #   Practically it removes node/layer associated with node whose chemical/toxic
    #####################

    #######################
    # if layer attribute is none, remove from further process

    dic_cluster_total_input_idx_MOD_vs_node_info_new = {}

    for total_input_idx_MOD, node_info in dic_cluster_total_input_idx_MOD_vs_node_info.items():
        node_layer_attribute = ''
        if len(total_input_idx_MOD.split('|')) > 0:
            node_layer_attribute = total_input_idx_MOD.split('|')[1]

        if node_layer_attribute != 'none':
            dic_cluster_total_input_idx_MOD_vs_node_info_new[total_input_idx_MOD] = node_info

    ####################
    #  if layer attribute is "none", remove from further process

    list_of_edge_for_networkx_new = []

    for edge in list_of_edge_for_networkx:

        node1_layer_attribute = ""
        node2_layer_attribute = ""

        if len(edge[0].split("|")) > 0:
            node1_layer_attribute = edge[0].split("|")[1]
        if len(edge[1].split("|")) > 0:
            node2_layer_attribute = edge[1].split("|")[1]

        flag_edge_to_remove = 0
        if node1_layer_attribute == "none" or node2_layer_attribute == "none":
            flag_edge_to_remove = 1

        if flag_edge_to_remove == 0:
            list_of_edge_for_networkx_new.append(edge)

    return dic_cluster_total_input_idx_MOD_vs_node_info_new, list_of_edge_for_networkx_new
