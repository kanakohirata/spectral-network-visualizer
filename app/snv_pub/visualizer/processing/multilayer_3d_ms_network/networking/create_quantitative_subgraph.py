from logging import getLogger
import networkx as nx


logger = getLogger(__name__)


def create_quantitative_subgraph(subgraph_type: str,
                                 subgraph_depth: int,
                                 subgraph_num_core_nodes: int,
                                 quant_polar: str,
                                 quant_value_type: str,
                                 quant_value: float,
                                 stat_value_threshold: float,
                                 dic_global_accession_vs_mass_feature: dict,
                                 list_dic_edges_nodes_graph_by_layer: list,
                                 fo_log_key) -> tuple:
    """

    Parameters
    ----------
    subgraph_type : str
    subgraph_depth : int
    subgraph_num_core_nodes : int
    quant_polar : str
        'quant_val_both', 'quant_val_decrease' or 'quant_val_increase'
    quant_value_type : str
        'ratio' or 'loading'
    quant_value : float
    stat_value_threshold : float
    dic_global_accession_vs_mass_feature : dict
    list_dic_edges_nodes_graph_by_layer : list
    fo_log_key : file object

    Returns
    -------
    list_dic_edges_nodes_graph_by_layer_new : list
        Updated list_dic_edges_nodes_graph_by_layer
    l_total_input_idx_mod_sample : list
        A list of extracted total_input_idx_mod of sample layers.
    """
    list_dic_edges_nodes_graph_by_layer_new = []
    l_total_input_idx_mod_sample = []

    # making dictionary where key is global accession and value is "mf['value_to_show']"
    dic_global_accession_vs_value_to_show = {}
    # this has p-value for mass feature
    dic_global_accession_vs_stat_value = {}

    for k, mf in dic_global_accession_vs_mass_feature.items():
        dic_global_accession_vs_value_to_show[k] = mf['value_to_show']
        dic_global_accession_vs_stat_value[k] = mf['stat_val']

    # foe each layer
    for dic_edges_nodes_graph_by_layer in list_dic_edges_nodes_graph_by_layer:

        fo_log_key.write('\ndic_edges_nodes_graph_by_layer[attribute_for_layer]: '
                         f'{dic_edges_nodes_graph_by_layer["attribute_for_layer"]}')

        if not dic_edges_nodes_graph_by_layer['attribute_for_layer'].startswith('sample'):
            list_dic_edges_nodes_graph_by_layer_new.append(dic_edges_nodes_graph_by_layer)
            continue

        # if the current layer is "sample"------------------------------------------------------
        # this list hold global accession for base of subgraph(quant)
        l_global_accession_for_subgraph_sample_user_selected_quant = []

        if subgraph_type in ['node_quant']:
            top_n = subgraph_num_core_nodes

            if quant_polar == 'quant_val_both':
                top_n = top_n / 2

            #######################################
            # select topN (bipolar)
            #######################################
            # small to large
            if quant_polar == 'quant_val_both' or quant_polar == 'quant_val_decrease':
                count = 0
                for k, v in sorted(list(dic_global_accession_vs_value_to_show.items()), key=lambda x: x[1]):
                    if count == top_n:
                        break

                    stat_val = dic_global_accession_vs_stat_value[k]

                    if quant_value_type == 'ratio':
                        threshold_quant_value = 1.0 / quant_value
                    else:
                        threshold_quant_value = - abs(quant_value)

                    # logger.warning(f'@@@D v: {v}, threshold_quant_value: {threshold_quant_value}, count: {count}')
                    # if statistical value is lower than threshold,  AND  the value is  lower than target RATIO. NOTE this is ratio.
                    if stat_val < stat_value_threshold and v < threshold_quant_value:
                        l_global_accession_for_subgraph_sample_user_selected_quant.append(k)
                        count = count + 1

            if quant_polar == 'quant_val_both' or quant_polar == 'quant_val_increase':

                fo_log_key.write('\n quant mode: both or val increase')
                count = 0
                # large to small
                for k, v in sorted(list(dic_global_accession_vs_value_to_show.items()),
                                   reverse=True, key=lambda x: x[1]):
                    if count == top_n:
                        break

                    stat_val = dic_global_accession_vs_stat_value[k]
                    if quant_value_type == 'ratio':
                        threshold_quant_value = quant_value
                    else:
                        threshold_quant_value = quant_value

                    # logger.warning(f'@@@I v: {v}, threshold_quant_value: {threshold_quant_value}, count: {count}')
                    # if statistical value is lower than threshold,  AND  the value is  higher than target RATIO. NOTE this is ratio.
                    if stat_val < stat_value_threshold and v > threshold_quant_value:
                        l_global_accession_for_subgraph_sample_user_selected_quant.append(k)
                        count = count + 1

            l_global_accession_for_subgraph_sample_user_selected_quant =\
                list(set(l_global_accession_for_subgraph_sample_user_selected_quant))

            fo_log_key.write('\nlen(l_global_accession_for_subgraph_sample_user_selected_quant): '
                             f'{len(l_global_accession_for_subgraph_sample_user_selected_quant)}')
            fo_log_key.flush()
            logger.warning(f"In '{dic_edges_nodes_graph_by_layer['attribute_for_layer']}' layer, "
                           f'length of l_global_accession_for_subgraph_sample_user_selected_quant: '
                           f'{len(l_global_accession_for_subgraph_sample_user_selected_quant)}')

            ####################################################################
            # get base for quant-based subgraph
            ####################################################################
            # create subgraph---------------------------------------------
            if len(l_global_accession_for_subgraph_sample_user_selected_quant) > 0:
                list_total_input_idx_mod_for_subgraph = []

                nx_graph = dic_edges_nodes_graph_by_layer['nx_graph']

                # l_global_accession_for_subgraph_sample_user_selected  is the list of "base" (starting point of subgraph)
                for node in nx_graph.nodes(data=True):
                    if node[1]['global_accession'] in l_global_accession_for_subgraph_sample_user_selected_quant:
                        list_total_input_idx_mod_for_subgraph.append(node[0])

                fo_log_key.write('\n len of base ( list_total_input_idx_mod_for_subgraph) : '
                                 f'{len(list_total_input_idx_mod_for_subgraph)}')
                fo_log_key.flush()
                logger.info(f"In '{dic_edges_nodes_graph_by_layer['attribute_for_layer']}' layer, "
                            f'length of base  (list_total_input_idx_mod_for_subgraph): '
                            f'{len(list_total_input_idx_mod_for_subgraph)}')
                base = list_total_input_idx_mod_for_subgraph

                # 'foundset' is a set of node names that are connected to nodes in 'base'
                # within a depth equal to or less than 'subgraph_depth'.
                foundset = {key for source in base for key in
                            list(nx.single_source_shortest_path(nx_graph, source,
                                                                cutoff=subgraph_depth).keys())}

                # update and replace nx graph and other info
                nx_graph_sub = nx_graph.subgraph(foundset)

                # create list of edges again.  (nx.generate_edgelist  not really working ???)----------------
                list_of_edge_for_networkx_sub = []
                for e in nx_graph_sub.edges.data():
                    list_of_edge_for_networkx_sub.append([e[0], e[1], e[2]])

                # updating "dic_cluster_total_input_idx_MOD_vs_node_info"-----------------------------
                # this has total input idx mod PRESENT in subgraph
                l_total_input_idx_mod_IN_subgraph = []
                for node in nx_graph_sub.nodes(data=True):
                    l_total_input_idx_mod_IN_subgraph.append(node[0])
                    l_total_input_idx_mod_sample.append(node[0])

                dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE = {}
                for cluster_total_input_idx_MOD, node_info in\
                        dic_edges_nodes_graph_by_layer['dic_cluster_total_input_idx_MOD_vs_node_info'].items():
                    if cluster_total_input_idx_MOD in l_total_input_idx_mod_IN_subgraph:
                        dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE[
                            cluster_total_input_idx_MOD] = node_info

                # create new dic_edges_nodes_graph_by_layer
                dic_edges_nodes_graph_by_layer_new = {
                    'attribute_for_layer': dic_edges_nodes_graph_by_layer['attribute_for_layer'],
                    'nx_graph': nx_graph_sub,
                    'list_of_edge_for_networkx': list_of_edge_for_networkx_sub,
                    'dic_cluster_total_input_idx_MOD_vs_node_info': dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE
                }
                list_dic_edges_nodes_graph_by_layer_new.append(dic_edges_nodes_graph_by_layer_new)

    return list_dic_edges_nodes_graph_by_layer_new, l_total_input_idx_mod_sample
