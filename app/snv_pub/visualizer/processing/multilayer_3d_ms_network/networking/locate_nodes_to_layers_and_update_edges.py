from logging import getLogger
from .. import dedicated_dictionaries


logger = getLogger(__name__)


def locate_nodes_to_layers_and_update_edges(dic_config,
                                            dic_cluster_total_input_idx_vs_cluster_info,
                                            list_edge_info) -> tuple:
    """

    Parameters
    ----------
    dic_config
    dic_cluster_total_input_idx_vs_cluster_info
    list_edge_info

    Returns
    -------
    dic_cluster_total_input_idx_MOD_vs_node_info : dict
        A dictionary of node_info which has total_input_idx_MOD (total_input_idx + layer name) as key.

        {
            '0|sample_aaa.mps':  {'total_input_idx': 0, 'total_input_idx_mod': 0|sample_aaa.mps, 'spec_cluster': cluster_info, 'layer': sample_aaa.mps, ...},
            '1|sample_aaa.mps':  {'total_input_idx': 1, 'total_input_idx_mod': 1|sample_aaa.mps, 'spec_cluster': cluster_info, 'layer': sample_aaa.mps, ...},
            '2|Benzenoids':  {'total_input_idx': 2, 'total_input_idx_mod': 2|Benzenoids, 'spec_cluster': cluster_info, 'layer': Benzenoids, ...},
            '3|Benzenoids':  {'total_input_idx': 3, 'total_input_idx_mod': 3|Benzenoids, 'spec_cluster': cluster_info, 'layer': Benzenoids, ...}, ...
        }
    list_edge_info : list
        A list of edge_info which has total_input_idx_MOD as 'spec_cluster_x_total_input_idx' and 'spec_cluster_y_total_input_idx'

        [
            {'spec_cluster_x_total_input_idx': '0|sample_aaa.mps', 'spec_cluster_y_total_input_idx': '1|sample_aaa.mps', 'spec_sim_score': 0.9, 'edge_type': 'inner_sample_layer'},
            {'spec_cluster_x_total_input_idx': '0|sample_aaa.mps', 'spec_cluster_y_total_input_idx': '2|Benzenoids', 'spec_sim_score': 0.8, 'edge_type': 'inter_sample_ref_layer'},
            {'spec_cluster_x_total_input_idx': '2|Benzenoids', 'spec_cluster_y_total_input_idx': '3|Benzenoids', 'spec_sim_score': 0.9, 'edge_type': 'inner_ref_layer'}, ...
        ]
    dic_cluster_id_vs_l_cluster_total_input_idx_MOD : dict

        {
            0 : ['0|sample_aaa.mps'],
            1: ['1|sample_aaa.mps'],
            2: ['2|Benzenoids'],
            3: ['3|Benzenoids'], ...
        }
    """
    ######################################################
    # [L] locate nodes to corresponding layers actually create NODE objects (NOT networkx node) here.
    # edit nodes according to layer specification
    # note if you have compound/cluster which belongs to two layers (metabolic pathway A and B), you have to DUPLICATE node entity in layer.
    ######################################################
    logger.debug('Start locate_nodes_to_layers_and_update_edges()')
    fo = open("log_locate_nodes_to_layers.txt", "w")

    dic_cluster_total_input_idx_MOD_vs_node_info = {}
    dic_cluster_id_vs_l_cluster_total_input_idx_MOD = {}

    fo.write(" len dic_cluster_total_input_idx_vs_cluster_info" + str(len(dic_cluster_total_input_idx_vs_cluster_info)))
    fo.write("conf_o.type_attribute_for_layer_separation:" + dic_config["type_attribute_for_layer_separation"])
    logger.info(f'Length of dic_cluster_total_input_idx_vs_cluster_info: '
                   f'{len(dic_cluster_total_input_idx_vs_cluster_info)}')
    logger.info(f'conf_o.type_attribute_for_layer_separation: {dic_config["type_attribute_for_layer_separation"]}')
    # classify node (total_input_idx) according to layer definition.-----------------
    # for each "total_input_idx" NOT MOD !.
    for total_input_idx, cl_info in dic_cluster_total_input_idx_vs_cluster_info.items():

        list_layer_name = []
        # if the node is from sample, just enter "sample"
        if cl_info["tag"] == "sample":
            filename = cl_info["represen_spec_uni"]["source_filename"].rstrip()
            list_layer_name = [f'sample_{filename}']
        else:
            # first you read what kind of attribute this compound belongs to.   like "metabolic pathway",  ["A", "B"]
            if dic_config["type_attribute_for_layer_separation"] == "list_cmpd_classification_superclass":

                fo.write("\n list_cmpd_classification_superclass")
                for l in cl_info["represen_spec_uni"]["list_cmpd_classification_superclass"]:
                    list_layer_name.append(l)
                # if there is no description, just add to none
                if len(cl_info["represen_spec_uni"]["list_cmpd_classification_superclass"]) == 0:
                    list_layer_name.append("none")

            elif dic_config["type_attribute_for_layer_separation"] == "list_cmpd_classification_class":
                fo.write("\n list_cmpd_classification_class")
                for l in cl_info["represen_spec_uni"]["list_cmpd_classification_class"]:
                    list_layer_name.append(l)
                # if there is no description, just add to none
                if len(cl_info["represen_spec_uni"]["list_cmpd_classification_class"]) == 0:
                    list_layer_name.append("none")

            elif dic_config["type_attribute_for_layer_separation"] == "list_compound_categories":
                for cat in cl_info["list_compound_categories"]:
                    list_layer_name.append(cat)
                # if there is no description, just add to none
                if len(cl_info["list_compound_categories"]) == 0:
                    list_layer_name.append("none")

            # TODO: Check attribute_for_layer
            elif dic_config["type_attribute_for_layer_separation"] == "source_filename":
                attribute_for_layer = str(cl_info["represen_spec_uni"]["source_filename"])

        l_total_input_idx_mod_for_this_node = []
        for layer_name in list_layer_name:
            total_input_idx_mod = f'{total_input_idx}|{layer_name}'

            l_total_input_idx_mod_for_this_node.append(total_input_idx_mod)

            dic_node = dedicated_dictionaries.get_initialized_dic_node_info()
            dic_node["total_input_idx"] = total_input_idx
            dic_node["total_input_idx_mod"] = total_input_idx_mod
            dic_node["spec_cluster"] = cl_info
            dic_node["layer"] = layer_name
            dic_node["suspect"] = ""
            dic_node["l_source_suspect"] = []
            dic_node["name"] = ""

            dic_cluster_total_input_idx_MOD_vs_node_info[total_input_idx_mod] = dic_node

        dic_cluster_id_vs_l_cluster_total_input_idx_MOD[cl_info["cluster_id"]] = l_total_input_idx_mod_for_this_node

    fo.write(f'\n  dic_cluster_total_input_idx_MOD_vs_node_info {len(dic_cluster_total_input_idx_MOD_vs_node_info)}')
    fo.write(f'\n dic_cluster_id_vs_l_cluster_total_input_idx_MOD '
             f'{len(dic_cluster_id_vs_l_cluster_total_input_idx_MOD)}')

    logger.info(f'Length of dic_cluster_total_input_idx_MOD_vs_node_info: '
                f'{len(dic_cluster_total_input_idx_MOD_vs_node_info)}')
    logger.info(f'Length of dic_cluster_id_vs_l_cluster_total_input_idx_MOD: '
                f'{len(dic_cluster_id_vs_l_cluster_total_input_idx_MOD)}')

    ############################################################
    # [N]Consequently, you have to update edge
    #    duplicate node and make edges/.
    ###########################################################
    logger.debug("Recreating edges with layer info")

    # you have to create list of edge info again for network X.

    # first, iterating original edge list (not duplicated according to layers.)
    list_edge_info_mod = []

    # for each edge (now no layer attribute)
    for edge_info in list_edge_info:
        # first, convert global accession to input total idx
        # these list will contain total idx (modified) of corresponding cluster.
        # note this is necessary sice total idx have been modified and duplicated for layers,

        if edge_info['spec_cluster_x_id'] in dic_cluster_id_vs_l_cluster_total_input_idx_MOD and \
                edge_info['spec_cluster_y_id'] in dic_cluster_id_vs_l_cluster_total_input_idx_MOD:

            # mod
            l_x_total_idx_mod =\
                dic_cluster_id_vs_l_cluster_total_input_idx_MOD[edge_info['spec_cluster_x_id']]
            l_y_total_idx_mod =\
                dic_cluster_id_vs_l_cluster_total_input_idx_MOD[edge_info['spec_cluster_y_id']]

            # make combination of X (total_idx_mod) with layer VS Y (total_idx_mod) with layer
            for x_idx_mod in l_x_total_idx_mod:

                for y_idx_mod in l_y_total_idx_mod:

                    dic_edge_info_update = dedicated_dictionaries.get_initialized_dic_edge_info()
                    dic_edge_info_update['spec_cluster_x_total_input_idx'] = x_idx_mod
                    dic_edge_info_update['spec_cluster_y_total_input_idx'] = y_idx_mod

                    dic_edge_info_update['spec_cluster_x_id'] = edge_info['spec_cluster_x_id']
                    dic_edge_info_update['spec_cluster_y_id'] = edge_info['spec_cluster_y_id']

                    dic_edge_info_update['spec_cluster_x_global_accession'] = edge_info['spec_cluster_x_global_accession']
                    dic_edge_info_update['spec_cluster_y_global_accession'] = edge_info['spec_cluster_y_global_accession']

                    dic_edge_info_update['spec_sim_score'] = edge_info['spec_sim_score']
                    dic_edge_info_update['delta_mz'] = edge_info['delta_mz']

                    layer_x = dic_cluster_total_input_idx_MOD_vs_node_info[x_idx_mod]['layer']
                    layer_y = dic_cluster_total_input_idx_MOD_vs_node_info[y_idx_mod]['layer']

                    edge_type = ''

                    f_x_sample = 0
                    f_y_sample = 0
                    if layer_x.startswith('sample'):
                        f_x_sample = 1
                    if layer_y.startswith('sample'):
                        f_y_sample = 1

                    # inter ref layer.    where edge is inside one of ref layers.-----
                    # if both nodes are not in sample
                    if f_x_sample == 0 and f_y_sample == 0:
                        if layer_x == layer_y:
                            edge_type = 'inner_ref_layer'
                        if layer_x != layer_y:
                            edge_type = 'inter_ref_layer'

                    elif (f_x_sample == 1 and f_y_sample == 0) or (f_x_sample == 0 and f_y_sample == 1):
                        edge_type = 'inter_sample_ref_layer'

                    elif f_x_sample == 1 and f_y_sample == 1 and layer_x == layer_y:
                        edge_type = 'inner_sample_layer'

                    elif f_x_sample == 1 and f_y_sample == 1 and layer_x != layer_y:
                        edge_type = 'inter_sample_layer'

                    if edge_type not in ['inner_ref_layer', 'inter_ref_layer', 'inner_sample_layer',
                                         'inter_sample_layer',
                                         'inter_sample_ref_layer']:
                        logger.debug('edge type is invalid')
                        dic_config['all_data_valid'] = False
                        dic_config['l_func_invalid_data'].append('multilayer_3d_network.locate_nodes_to_layers_and_update_edges')

                    dic_edge_info_update['edge_type'] = edge_type
                    list_edge_info_mod.append(dic_edge_info_update)

    # take over
    list_edge_info = list_edge_info_mod

    fo.write(f'\n finishing function dic_cluster_total_input_idx_MOD_vs_node_info'
             f'{len(dic_cluster_total_input_idx_MOD_vs_node_info)}')
    logger.info(f'Length of dic_cluster_total_input_idx_MOD_vs_node_info: '
                f'{len(dic_cluster_total_input_idx_MOD_vs_node_info)}')
    logger.debug('finishing Recreating edges with layer info')
    return dic_cluster_total_input_idx_MOD_vs_node_info, list_edge_info, dic_cluster_id_vs_l_cluster_total_input_idx_MOD

