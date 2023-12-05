from logging import getLogger
import networkx as nx


logger = getLogger(__name__)


def extract_subgraph_based_on_sample_global_accession(l_global_accession_for_subgraph_sample_user_selected,
                                                      node_select_subgraph_depth,
                                                      list_dic_edges_nodes_graph_by_layer,
                                                      nx_graph,
                                                      fo_log):
    """

    Parameters
    ----------
    l_global_accession_for_subgraph_sample_user_selected
    node_select_subgraph_depth
    list_dic_edges_nodes_graph_by_layer
    nx_graph
    fo_log

    Returns
    -------
    list_dic_edges_nodes_graph_by_layer_new : list
        Updated list_dic_edges_nodes_graph_by_layer
    FG_all_samples : networkx.Graph
        A networkx Graph of extracted samples.
    l_total_input_idx_mod_sample : list
    """
    l_total_idx_mod_user_select_subgraph_all_sample = []
    l_total_input_idx_mod_sample = []

    list_total_input_idx_mod_for_subgraph = []
    fo_log.write(f'nx_graph.data(): {nx_graph.edges.data()}')
    logger.info(f'nx_graph.data(): {nx_graph.edges.data()}')

    # l_global_accession_for_subgraph_sample_user_selected  is the list of "base" (starting point of subgraph)
    # node[0] is total_input_idx_mod
    for node in nx_graph.nodes(data=True):
        fo_log.write(f'\n  FG_all_samples node:{node}')
        logger.info(f'FG_all_samples node: {str(node)}')
        if node[1]['global_accession'] in l_global_accession_for_subgraph_sample_user_selected:
            list_total_input_idx_mod_for_subgraph.append(node[0])

    fo_log.write(f'CREATING SUER DEFINED SUBGRAPH:subgraph_depth:{node_select_subgraph_depth}'
                 f'   base:{list_total_input_idx_mod_for_subgraph}\n')
    fo_log.flush()
    logger.info(f'CREATING SUER DEFINED SUBGRAPH'
                f'\nconfig_o.node_select_subgraph_depth: {node_select_subgraph_depth}'
                f'\nbase: {list_total_input_idx_mod_for_subgraph}')

    # 'foundset' is a set of node names that are connected to nodes in 'list_total_input_idx_mod_for_subgraph'
    # within a depth equal to or less than 'subgraph_depth'.
    foundset = {key for source in list_total_input_idx_mod_for_subgraph for key in
                list(nx.single_source_shortest_path(nx_graph, source, cutoff=node_select_subgraph_depth).keys())}

    # update and replace nx graph and other info
    fo_log.write(f'\n\nfound set:{foundset}')
    logger.info(f'Found set: {foundset}')
    FG_all_samples = nx_graph.subgraph(foundset)

    # l_total_idx_mod_user_select_subgraph_all_sample holds all idx mod to be extracted
    for n in FG_all_samples.nodes.data():
        l_total_idx_mod_user_select_subgraph_all_sample.append(n[0])

    # Update list_dic_edges_nodes_graph_by_layer
    list_dic_edges_nodes_graph_by_layer_new = []

    # foe each layer
    for dic_edges_nodes_graph_by_layer in list_dic_edges_nodes_graph_by_layer:
        if not dic_edges_nodes_graph_by_layer['attribute_for_layer'].startswith('sample'):
            list_dic_edges_nodes_graph_by_layer_new.append(dic_edges_nodes_graph_by_layer)

        # if the current layer is "sample"------------------------------------------------------
        else:
            # create subgraph ---------------------------------------------
            _nx_graph = dic_edges_nodes_graph_by_layer["nx_graph"]
            nx_graph_sub = _nx_graph.subgraph(l_total_idx_mod_user_select_subgraph_all_sample)

            # create list of edges again.  nx.generate_edgelist  not really working ???
            list_of_edge_for_networkx_sub = []
            for e in nx_graph_sub.edges.data():
                list_of_edge_for_networkx_sub.append([e[0], e[1], e[2]])

            # updating "dic_cluster_total_input_idx_MOD_vs_node_info" -----------------------------
            # this has total input idx mod PRESENT in subgraph
            l_total_input_idx_mod_IN_subgraph = []

            for node in nx_graph_sub.nodes(data=True):
                l_total_input_idx_mod_IN_subgraph.append(node[0])
                l_total_input_idx_mod_sample.append(node[0])

            dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE = {}
            for cluster_total_input_idx_MOD, node_info in\
                    dic_edges_nodes_graph_by_layer['dic_cluster_total_input_idx_MOD_vs_node_info'].items():
                if cluster_total_input_idx_MOD in l_total_input_idx_mod_IN_subgraph:
                    dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE[cluster_total_input_idx_MOD] = node_info

            # create new dic_edges_nodes_graph_by_layer
            dic_edges_nodes_graph_by_layer_new = {
                'attribute_for_layer': dic_edges_nodes_graph_by_layer['attribute_for_layer'],
                'nx_graph': nx_graph_sub,
                'list_of_edge_for_networkx': list_of_edge_for_networkx_sub,
                'dic_cluster_total_input_idx_MOD_vs_node_info': dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE}
            list_dic_edges_nodes_graph_by_layer_new.append(dic_edges_nodes_graph_by_layer_new)

    return list_dic_edges_nodes_graph_by_layer_new, FG_all_samples, l_total_input_idx_mod_sample


def extract_ref_subgraph_connected_to_sample(list_dic_edges_nodes_graph_by_layer: list,
                                             list_total_input_idx_MOD_inter_sample_ref_layer: list,
                                             subgraph_depth: int):
    """
    Parameters
    ----------
    list_dic_edges_nodes_graph_by_layer : list
    list_total_input_idx_MOD_inter_sample_ref_layer : list
    subgraph_depth : int

    Returns
    -------
    list_dic_edges_nodes_graph_by_layer_new : list
        Updated list_dic_edges_nodes_graph_by_layer
    """
    # Update list_dic_edges_nodes_graph_by_layer
    list_dic_edges_nodes_graph_by_layer_new = []

    # look for ref layer---------------------------------
    for dic_edges_nodes_graph_by_layer in list_dic_edges_nodes_graph_by_layer:
        if dic_edges_nodes_graph_by_layer['attribute_for_layer'].startswith('sample'):
            list_dic_edges_nodes_graph_by_layer_new.append(dic_edges_nodes_graph_by_layer)
            continue

        # choosing ref layer (layer tag is not sample.)
        # selecting node
        list_total_input_idx_mod_for_base = []
        for cluster_total_input_idx_MOD, node_info in \
                dic_edges_nodes_graph_by_layer['dic_cluster_total_input_idx_MOD_vs_node_info'].items():

            # if the ref node is in inter ref-sample node list
            if cluster_total_input_idx_MOD in list_total_input_idx_MOD_inter_sample_ref_layer:
                list_total_input_idx_mod_for_base.append(cluster_total_input_idx_MOD)

        # creating subgraph.
        nx_graph = dic_edges_nodes_graph_by_layer['nx_graph']
        # 'foundset' is a set of node names that are connected to nodes in 'list_total_input_idx_mod_for_base'
        # within a depth equal to or less than 'subgraph_depth'.
        foundset = {key for source in list_total_input_idx_mod_for_base for key in
                    list(nx.single_source_shortest_path(nx_graph, source,cutoff=subgraph_depth).keys())}
        nx_graph_sub = nx_graph.subgraph(foundset)

        # create list of edges again.  nx.generate_edgelist  not really working ???
        list_of_edge_for_networkx_sub = []
        for e in nx_graph_sub.edges.data():
            list_of_edge_for_networkx_sub.append([e[0], e[1], e[2]])

        # this has total input idx mod PRESENT in subgraph
        l_total_input_idx_mod_IN_subgraph = []

        for node in nx_graph_sub.nodes(data=True):
            l_total_input_idx_mod_IN_subgraph.append(node[0])

        dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE = {}
        for cluster_total_input_idx_MOD, node_info in \
                dic_edges_nodes_graph_by_layer['dic_cluster_total_input_idx_MOD_vs_node_info'].items():
            if cluster_total_input_idx_MOD in l_total_input_idx_mod_IN_subgraph:
                dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE[cluster_total_input_idx_MOD] = node_info

        dic_edges_nodes_graph_by_layer_new = {
            'attribute_for_layer': dic_edges_nodes_graph_by_layer['attribute_for_layer'],
            'nx_graph': nx_graph_sub,
            'list_of_edge_for_networkx': list_of_edge_for_networkx_sub,
            'dic_cluster_total_input_idx_MOD_vs_node_info': dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE}

        list_dic_edges_nodes_graph_by_layer_new.append(dic_edges_nodes_graph_by_layer_new)

    return list_dic_edges_nodes_graph_by_layer_new
