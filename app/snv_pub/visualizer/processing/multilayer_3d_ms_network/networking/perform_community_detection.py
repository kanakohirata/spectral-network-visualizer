import json
from logging import getLogger
import networkx as nx
from networkx.algorithms import community
from networkx.algorithms.community import greedy_modularity_communities
from networkx.readwrite import json_graph


logger = getLogger(__name__)


def perform_community_detection(n_level_community_detection: int,
                                list_dic_edges_nodes_graph_by_layer: list,
                                list_of_edge_for_networkx_inter_sample_ref_layer: list,
                                list_of_edge_for_networkx_inter_sample_layer: list,
                                fo_log_key):
    fo_log_key.write(f'\n performing community detection at level : {n_level_community_detection}')
    logger.info(f'Performing community detection at level: {n_level_community_detection}')

    str_log = '\n\nCommunity detection===============\n'
    #####
    # [O1] first,  combine all nodes, edges ---------------------------------
    logger.debug('  [O1] combine all nodes, edges ')
    G_all_combined = nx.Graph()

    # combine node/edges for each layer. (not inter layer edges)
    for dic_edges_nodes_graph_by_layer in list_dic_edges_nodes_graph_by_layer:
        G_all_combined = nx.compose(G_all_combined, dic_edges_nodes_graph_by_layer["nx_graph"])

    # add inter sample-ref edges
    G_all_combined.add_edges_from(list_of_edge_for_networkx_inter_sample_ref_layer)

    # add inter sample-ref edges
    G_all_combined.add_edges_from(list_of_edge_for_networkx_inter_sample_layer)

    # [O3] Perform community detection ----------------------------------------------
    logger.debug('  [O3] Perform community detection')

    l_communities = []
    community_detection_method = 'greedy_modularity_communities'

    if community_detection_method == 'girvan_newman':
        communities_generator = community.girvan_newman(G_all_combined)
        # return as tuple
        for n in range(n_level_community_detection):
            l_communities = next(communities_generator)

    if community_detection_method == 'greedy_modularity_communities':
        try:
            l_communities = list(greedy_modularity_communities(G_all_combined))
        except:
            pass

    logger.debug('  [O3] community detection finished')
    logger.debug('  [O3] update edges and graph based on the detected communities')
    list_dic_edges_nodes_graph_by_layer_new = []
    if len(l_communities) > 0:

        str_log += '\n detected communities \n'
        for c in l_communities:
            str_log += f'\n{c}'

        # [O5]  check all the edges whether it fits community definition -----------------

        # this holds edges that do not fit within community.  in other word, inter community edges
        l_inter_community_edge_as_tuple = []
        logger.debug('         iterate edges in G_all_combined')
        # iterate edges----------------------
        for u, v, d in G_all_combined.edges(data=True):
            f_in_community = 0
            # iterate communities
            for commu in l_communities:
                # !!!!!!!!!! 20220903 changed
                if u in commu and v in commu:
                    f_in_community = 1
                    break
            # if the current edge is not inside community,
            if f_in_community == 0:
                l_inter_community_edge_as_tuple.append((u, v, d))

        logger.debug(f'inter-community edges {len(l_inter_community_edge_as_tuple)}')
        str_log += f'\n inter-community edges {len(l_inter_community_edge_as_tuple)}\n'

        fo_log_key.write(str_log)
        fo_log_key.flush()
        logger.debug(str_log)

        # [O7]: update edge info-------------------------------------------------------------------
        logger.debug('  [O7]  update edge info')
        # for each layer, remove edges outside communities and regenerate networkx graph.
        for dic_edges_nodes_graph_by_layer in list_dic_edges_nodes_graph_by_layer:
            l_inter_community_edge = []
            # list of edge within community, which supposed to be preserved.
            l_inner_community_edge = []
            # iterate edges
            for edge in dic_edges_nodes_graph_by_layer['list_of_edge_for_networkx']:
                f_in_community = 0
                # check the edge is in community.
                for commu in l_communities:

                    # !!!!!!!!!! 20220903 changed
                    if edge[0] in commu and edge[1] in commu:
                        l_inner_community_edge.append((edge[0], edge[1], edge[2]))
                        f_in_community = 1

                if f_in_community == 0:
                    l_inter_community_edge.append((edge[0], edge[1], edge[2]))

            # it seems you need to unfreeze graph
            nx_graph_unfrozen = nx.Graph(dic_edges_nodes_graph_by_layer['nx_graph'])
            # then remove edge.
            nx_graph_unfrozen.remove_edges_from(l_inter_community_edge)

            dic_edges_nodes_graph_by_layer_new = {
                'attribute_for_layer': dic_edges_nodes_graph_by_layer['attribute_for_layer'],
                'list_of_edge_for_networkx': l_inner_community_edge,
                'dic_cluster_total_input_idx_MOD_vs_node_info':
                    dic_edges_nodes_graph_by_layer['dic_cluster_total_input_idx_MOD_vs_node_info'],
                'nx_graph': nx_graph_unfrozen
            }
            list_dic_edges_nodes_graph_by_layer_new.append(dic_edges_nodes_graph_by_layer_new)

        # remove inter community edges for list_of_edge_for_networkx_to_show_inter_sample_ref_layer
        list_of_edge_for_networkx_inter_sample_ref_layer_new = []
        for edge in list_of_edge_for_networkx_inter_sample_ref_layer:
            # check the edge is in community.
            for commu in l_communities:
                if all(x in commu for x in [edge[0], edge[1]]):
                    list_of_edge_for_networkx_inter_sample_ref_layer_new.append((edge[0], edge[1], edge[2]))

        # logger.debug("  [O7]  remove inter community edges  for list_of_edge_for_networkx_to_show_inter_sample_layer ")
        # remove inter community edges  for list_of_edge_for_networkx_to_show_inter_sample_layer
        list_of_edge_for_networkx_inter_sample_layer_new = []
        for edge in list_of_edge_for_networkx_inter_sample_layer:
            # check the edge is in community.
            for commu in l_communities:
                if all(x in commu for x in [edge[0], edge[1]]):
                    list_of_edge_for_networkx_inter_sample_layer_new.append((edge[0], edge[1], edge[2]))

    else:
        list_dic_edges_nodes_graph_by_layer_new = list_dic_edges_nodes_graph_by_layer
        list_of_edge_for_networkx_inter_sample_ref_layer_new = list_of_edge_for_networkx_inter_sample_ref_layer
        list_of_edge_for_networkx_inter_sample_layer_new = list_of_edge_for_networkx_inter_sample_layer
    # export
    logger.debug("    exporting to json")
    data1 = json_graph.node_link_data(G_all_combined)
    fp = open("G_all_combined.json", "w")
    json.dump(data1, fp)
    logger.debug("    finished exporting to json")
    logger.debug(" Finishing community detection")

    return (list_dic_edges_nodes_graph_by_layer_new,
            list_of_edge_for_networkx_inter_sample_ref_layer_new,
            list_of_edge_for_networkx_inter_sample_layer_new)
