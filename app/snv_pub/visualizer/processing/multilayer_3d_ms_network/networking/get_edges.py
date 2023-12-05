from logging import getLogger


logger = getLogger(__name__)


def update_edges(list_dic_edges_nodes_graph_by_layer,
                 list_of_edge_for_networkx_inter_sample_ref_layer,
                 list_of_edge_for_networkx_inter_sample_layer):
    """

    Parameters
    ----------
    list_dic_edges_nodes_graph_by_layer
    list_of_edge_for_networkx_inter_sample_ref_layer
    list_of_edge_for_networkx_inter_sample_layer

    Returns
    -------
    list_of_edge_for_networkx_inter_sample_ref_layer_new : list
        Updated list_of_edge_for_networkx_inter_sample_ref_layer
    list_of_edge_for_networkx_inter_sample_layer_new : list
        Updated list_of_edge_for_networkx_inter_sample_layer
    """
    l_cluster_total_input_idx_MOD_in_all_subgraph = []
    # foe each layer
    for dic_edges_nodes_graph_by_layer in list_dic_edges_nodes_graph_by_layer:
        # if the current layer is sample
        for cluster_total_input_idx_MOD, v in\
                dic_edges_nodes_graph_by_layer['dic_cluster_total_input_idx_MOD_vs_node_info'].items():
            l_cluster_total_input_idx_MOD_in_all_subgraph.append(cluster_total_input_idx_MOD)

    # get edges both nodes can be found in current subgraph:
    list_of_edge_for_networkx_inter_sample_ref_layer_new = []
    for edge_for_networkx in list_of_edge_for_networkx_inter_sample_ref_layer:

        if (edge_for_networkx[0] in l_cluster_total_input_idx_MOD_in_all_subgraph
                and edge_for_networkx[1] in l_cluster_total_input_idx_MOD_in_all_subgraph):
            list_of_edge_for_networkx_inter_sample_ref_layer_new.append(edge_for_networkx)

    list_of_edge_for_networkx_inter_sample_layer_new = []
    for edge_for_networkx in list_of_edge_for_networkx_inter_sample_layer:
        if (edge_for_networkx[0] in l_cluster_total_input_idx_MOD_in_all_subgraph
                and edge_for_networkx[1] in l_cluster_total_input_idx_MOD_in_all_subgraph):
            list_of_edge_for_networkx_inter_sample_layer_new.append(edge_for_networkx)

    return list_of_edge_for_networkx_inter_sample_ref_layer_new, list_of_edge_for_networkx_inter_sample_layer_new
