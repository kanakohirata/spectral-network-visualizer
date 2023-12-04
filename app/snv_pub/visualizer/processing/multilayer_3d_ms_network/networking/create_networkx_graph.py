from logging import getLogger
import networkx as nx


logger = getLogger(__name__)


def create_networkx_graph(dic_cluster_total_input_idx_MOD_vs_node_info: dict, list_of_edge_for_networkx: list):
    """

    Parameters
    ----------
    dic_cluster_total_input_idx_MOD_vs_node_info : dict
    list_of_edge_for_networkx : list

    Returns
    -------
    FG : networkx.Graph
    """

    ##################
    # [U] create graph

    FG = nx.Graph()
    FG.add_edges_from(list_of_edge_for_networkx)

    #############
    # add node
    #############

    for cluster_total_input_idx_MOD, node_info in dic_cluster_total_input_idx_MOD_vs_node_info.items():
        FG.add_node(cluster_total_input_idx_MOD)

    #############################
    # create list od node id here
    #############################

    list_total_input_idx_MOD_used = [node[0] for node in FG.nodes.data()]

    # making list of nodes

    list_attribute_for_layer = []

    ##################################################
    # iterate list of node and add node with attribute
    ##################################################

    for total_input_idx_mod, node_info in dic_cluster_total_input_idx_MOD_vs_node_info.items():

        # this "attribute_for_layer"  keeps tag used for layer separation

        # here you are duplicating node IF there is more than one layer attribute is possible,
        # for example this node compound is A and B metabolite pathway.

        # add node only if the node is present in edge list.
        if total_input_idx_mod in list_total_input_idx_MOD_used:
            FG.add_node(total_input_idx_mod, attribute_for_layer=node_info["layer"],
                        global_accession=node_info["spec_cluster"]["global_accession"])

            list_attribute_for_layer.append(node_info["layer"])
            list_attribute_for_layer = list(set(list_attribute_for_layer))

    logger.debug(f"list_attribute_for_layer {list_attribute_for_layer}")
    logger.debug("finishing  [create_networkx_graph]")

    return FG


def create_sample_networkx_graph(list_dic_edges_nodes_graph_by_layer: list,
                                 list_of_edge_for_networkx_inter_sample_layer: list,
                                 fo_log):
    """

    Parameters
    ----------
    list_dic_edges_nodes_graph_by_layer : list
    list_of_edge_for_networkx_inter_sample_layer : list
    fo_log : file object

    Returns
    -------
    FG_all_samples : networkx.Graph
    """
    # first, assemble all nodes, edges from sample tag containing dataset
    dic_cluster_total_input_idx_MOD_vs_node_info_all_samples = {}
    list_of_edge_for_networkx_all_samples = []

    for dic in list_dic_edges_nodes_graph_by_layer:
        # dealing with sample layer data.
        # note if you have more than 1, you have to merge and make layout for the merged dataset.
        if dic["attribute_for_layer"].startswith("sample"):
            fo_log.write(f'\n\niterating : {dic["attribute_for_layer"]}')
            fo_log.write("\ndic[dic_cluster_total_input_idx_MOD_vs_node_info]: "
                         f'{dic["dic_cluster_total_input_idx_MOD_vs_node_info"]}')
            fo_log.write(f'\ndic[list_of_edge_for_networkx]: {dic["list_of_edge_for_networkx"]}')
            logger.info(f'Iterating: {dic["attribute_for_layer"]}\n'
                        f'\ndic["dic_cluster_total_input_idx_MOD_vs_node_info"]: '
                        f'{dic["dic_cluster_total_input_idx_MOD_vs_node_info"]}'
                        f'\ndic[list_of_edge_for_networkx]: {dic["list_of_edge_for_networkx"]}')

            dic_cluster_total_input_idx_MOD_vs_node_info_all_samples.update(
                dic["dic_cluster_total_input_idx_MOD_vs_node_info"])
            list_of_edge_for_networkx_all_samples += dic["list_of_edge_for_networkx"]
    list_of_edge_for_networkx_all_samples += list_of_edge_for_networkx_inter_sample_layer
    fo_log.write(f"\n\nlist_of_edge_for_networkx_all_samples: {list_of_edge_for_networkx_all_samples}")
    logger.info(f'list_of_edge_for_networkx_all_samples: {list_of_edge_for_networkx_all_samples}')

    # just converting to networkx graph.  doing nothing with subgraph etc.
    FG_all_samples = create_networkx_graph(dic_cluster_total_input_idx_MOD_vs_node_info_all_samples,
                                           list_of_edge_for_networkx_all_samples)

    return FG_all_samples
