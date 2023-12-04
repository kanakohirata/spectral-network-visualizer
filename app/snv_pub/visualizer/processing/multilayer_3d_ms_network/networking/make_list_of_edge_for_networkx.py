from logging import getLogger


logger = getLogger(__name__)


def make_list_of_edge_for_networkx(list_edge_info):
    """
    Parameters
    ----------
    list_edge_info

    Returns
    -------
    list_of_edge_for_networkx : list
        A list of lists which have 'spec_cluster_x_total_input_idx' at 0th, 'spec_cluster_y_total_input_idx' at 1st and
        a dictionary of 'spec_sim_score', 'delta_mz' and 'edge_type' at 2nd.

        [
            ['0|sample_aaa.msp', '1|sample_aaa.msp', {'spec_sim_score': 0.9, 'delta_mz': 0.01, 'edge_type': 'inner_sample_layer'}],
            ['0|sample_aaa.msp', '2|Benzenoids', {'spec_sim_score': 0.8, 'delta_mz': 0.02, 'edge_type': 'inter_sample_ref_layer'}],
            ['2|Benzenoids', '3|Benzenoids', {'spec_sim_score': 0.9, 'delta_mz': 0.01, 'edge_type': 'inner_ref_layer'}], ...
        ]
    """
    list_of_edge_for_networkx = []
    for edge_o in list_edge_info:
        if edge_o['edge_type'] in ['inner_ref_layer', 'inner_sample_layer',
                                   'inter_sample_ref_layer', 'inter_sample_layer']:

            # make list corresponds the info of one node.
            list_for_one_node = [
                edge_o['spec_cluster_x_total_input_idx'],
                edge_o['spec_cluster_y_total_input_idx'],
                {'spec_sim_score': edge_o['spec_sim_score'],
                 'delta_mz': edge_o['delta_mz'],
                 'edge_type': edge_o['edge_type']}
            ]

            list_of_edge_for_networkx.append(list_for_one_node)

    return list_of_edge_for_networkx
