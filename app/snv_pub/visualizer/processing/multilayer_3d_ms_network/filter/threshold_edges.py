from logging import getLogger

logger = getLogger(__name__)


def threshold_edges(score_threshold, list_of_edge_for_networkx):
    ###########################
    #  [T] note. now you used threshold to select edges.
    #  which means,  some of the nodes in node list is not used in edge list.
    #  This will cause trouble later.
    #  so now you can select nodes actually present in edge list and update.

    logger.debug("function [threshold_edges] select edges using threshold")

    list_of_edge_for_networkx_new = []
    list_node_total_input_idx_mod_in_use = []

    # dic_cluster_total_input_idx_vs_cluster_info
    for edge in list_of_edge_for_networkx:

        if edge[2]['spec_sim_score'] > score_threshold:
            list_of_edge_for_networkx_new.append(edge)

            if edge[0] not in list_node_total_input_idx_mod_in_use:
                list_node_total_input_idx_mod_in_use.append(edge[0])
            if edge[1] not in list_node_total_input_idx_mod_in_use:
                list_node_total_input_idx_mod_in_use.append(edge[1])

    logger.debug("Finished function [threshold_edges] ")
    return list_of_edge_for_networkx_new, list_node_total_input_idx_mod_in_use
