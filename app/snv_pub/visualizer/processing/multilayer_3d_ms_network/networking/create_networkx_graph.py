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
