# TO DO
#  "select_nodes_based_on_product_mz_required"
#           currently only examine representative spec.
import copy
import glob
import os.path
import sys
import zipfile
from logging import getLogger
import networkx as nx
import pandas as pd
import plotly.graph_objs as go
from django.core.cache import cache

from . import dedicated_dictionaries
from . import feature_table_reader
from . import multilayer_3d_mesh_functsions as m3d_mesh
from . import multilayer_3d_rescale_functions as m3d_rescale
from . import read_t3db_a1
from . import suspect_compound
from .my_parser.cluster_attribute_parser import read_cluster_attribute
from .my_parser.compound_info_parser import add_external_cmpd_info
from .my_parser.edge_info_parser import read_edge_info
from .my_parser.feature_table_parser import read_feature_table
from .utils import add_color_to_t3db_compound


logger = getLogger(__name__)


def select_nodes_for_subgraph(conf_o):
    #######################
    #   [J2] make list of node for subgraph
    ##################################

    #############
    #  here you are making list of global_accession of which mass feature show significance quantitative shift.
    l_global_accession_for_subgraph = []

    if subgraph_type in ["node_quant"]:
        # making dictionary where key is global accession and value is "mf.value_to_show"
        dic_global_accession_vs_value_to_show = {}
        for k, mf in dic_global_accession_vs_mass_feature.items():
            dic_global_accession_vs_value_to_show[k] = mf['value_to_show']

        top_n = subgraph_num_core_nodes
        count = 0

        for k, v in sorted(list(dic_global_accession_vs_value_to_show.items()), key=lambda x: x[1]):
            if count < top_n:
                l_global_accession_for_subgraph.append(k)
            count = count + 1

        count = 0
        for k, v in sorted(list(dic_global_accession_vs_value_to_show.items()), reverse=True, key=lambda x: x[1]):
            if count < top_n:
                l_global_accession_for_subgraph.append(k)
            count = count + 1

        l_global_accession_for_subgraph = list(set(l_global_accession_for_subgraph))

        # print  "l_global_accession_for_subgraph" ,l_global_accession_for_subgraph


def select_nodes_based_on_keyword(filter_select_category, filter_select_keyword, dic_cluster_total_input_idx_vs_cluster_info):
    # select specific node .
    logger.debug('[G2] select specific node/cluster')
    logger.debug(f'dic_config["filter_select_category"]: {filter_select_category}')

    if filter_select_category == 'list_cmpd_classification_superclass':
        dic_cluster_total_input_idx_vs_cluster_info_new = {}

        for total_input_idx, cl_info in dic_cluster_total_input_idx_vs_cluster_info.items():
            # if the current cluster is sample OR to be selected.
            if (cl_info['tag'] == 'sample'
                    or filter_select_keyword in cl_info['represen_spec_uni']['list_cmpd_classification_superclass']):
                dic_cluster_total_input_idx_vs_cluster_info_new[total_input_idx] = cl_info

        return dic_cluster_total_input_idx_vs_cluster_info_new

    else:
        return dic_cluster_total_input_idx_vs_cluster_info


def select_nodes_based_on_prec_mz(conf_o, dic_cluster_total_input_idx_vs_cluster_info):
    # select specific node .
    logger.debug("[G2] select specific node/cluster based on prec mz")
    # if conf_o.filter_select_category == "none":
    #    print " no keyword selected for node selection"
    #    return dic_cluster_total_input_idx_vs_cluster_info

    # if conf_o.filter_select_category != "none" :
    dic_cluster_total_input_idx_vs_cluster_info_new = {}

    for total_input_idx, cl_o in dic_cluster_total_input_idx_vs_cluster_info.items():
        if cl_o["represen_spec_uni"]["precursor_mz"] >= conf_o["mass_lower_limit"] and cl_o["represen_spec_uni"][
            "precursor_mz"] <= conf_o["mass_higher_limit"]:
            dic_cluster_total_input_idx_vs_cluster_info_new[total_input_idx] = cl_o

    return dic_cluster_total_input_idx_vs_cluster_info_new


# MODIFIED20220714---------------
def select_nodes_based_on_mass_defect(conf_o, dic_cluster_total_input_idx_vs_cluster_info):
    # mz tolerance to examine mass defect
    mz_tol = conf_o.mz_tol
    dic_cluster_total_input_idx_vs_cluster_info_new = {}
    count_mass_defect_match = 0
    count_sample = 0
    count_dic = 0

    logger.debug(f"conf_o.list_mass_defect {conf_o.list_mass_defect}")

    for total_input_idx, cl_o in dic_cluster_total_input_idx_vs_cluster_info.items():
        count_dic = count_dic + 1
        # for ref spec data, you can just mass
        f_add_to_new_dic = 1

        # for sample data, you have to examine mass defect
        if cl_o.tag.startswith("sample"):

            f_add_to_new_dic = 0
            count_sample = count_sample + 1

            # get mass defect for current node.
            node_mass_defect = cl_o.represen_spec_uni.precursor_mz - float(int(cl_o.represen_spec_uni.precursor_mz))
            # print( cl_o.represen_spec_uni.precursor_mz , " md: ",node_mass_defect)

            # for all mass defect value specified in config.
            for mass_df_specified in conf_o.list_mass_defect:
                # examine mass defect match to config val

                # print(  "    abs(node_mass_defect - mass_df_specified)" , abs(node_mass_defect - mass_df_specified))
                if abs(node_mass_defect - mass_df_specified) < mz_tol:
                    count_mass_defect_match = count_mass_defect_match + 1
                    f_add_to_new_dic = 1
                    # print("---------------mass defect matched")
        # if mass defect match, add to new dictionary
        if f_add_to_new_dic == 1:
            dic_cluster_total_input_idx_vs_cluster_info_new[total_input_idx] = cl_o
    logger.debug(f"count dic {str(count_dic)}")
    logger.debug(f"count samp;;e {str(count_sample)}")
    logger.debug(f" matched mass defect: {str(count_mass_defect_match)}")
    logger.debug(f"len of dic after mass defect exam: {len(dic_cluster_total_input_idx_vs_cluster_info_new)}")

    return dic_cluster_total_input_idx_vs_cluster_info_new


# ------------MODIFIED20220714


# MODIFIED20220719---------------

def select_nodes_based_on_product_mz_required(conf_o, dic_cluster_total_input_idx_vs_cluster_info):
    logger.debug("starting select_nodes_based_on_product_mz_required")
    # mz tolerance to examine mass product ion mz
    mz_tol = 0.05
    dic_cluster_total_input_idx_vs_cluster_info_new = {}
    for total_input_idx, cl_o in dic_cluster_total_input_idx_vs_cluster_info.items():

        # for ref spec data, you can just mass
        f_add_to_new_dic = 1

        # you have to examine product ion mz only for sample data.
        if cl_o.tag.startswith("sample"):
            f_add_to_new_dic = 0
            logger.debug(f" sample, total_input_idx {total_input_idx}")

            # for peaks in represen spec cluster-------------------------
            for pk in cl_o.represen_spec_uni.peak_list_mz_int_rel:

                # for all mass defect value specified in config.
                for product_mz_specified in conf_o.list_product_mz_required:
                    # examine if mz match to input.
                    if abs(pk[0] - product_mz_specified) < mz_tol:
                        f_add_to_new_dic = 1
                        # print("---------------product mz matched:"   ,pk[0], product_mz_specified  )
        # if mass defect match, add to new dictionary
        if f_add_to_new_dic == 1:
            dic_cluster_total_input_idx_vs_cluster_info_new[total_input_idx] = cl_o

    logger.debug("finishing select_nodes_based_on_product_mz_required")
    return dic_cluster_total_input_idx_vs_cluster_info_new


def locate_nodes_to_layers_and_update_edges(dic_config, dic_cluster_total_input_idx_vs_cluster_info, list_edge_info):
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
            {'spec_cluster_x_total_input_idx': '2|Benzenoids', 'spec_cluster_y_total_input_idx': '3|Benzenoids', 'spec_sim_score': 0.9, 'edge_type': 'inner_ref_layer'}
        ]
    dic_cluster_id_vs_l_cluster_total_input_idx_MOD : dict

        {
            0 : ['0|sample_aaa.mps'],
            1: ['1|sample_aaa.mps'],
            2: ['2|Benzenoids'],
            3: ['3|Benzenoids'],
        }
    """
    ######################################################
    # [L] locate nodes to corresponding layers actually create NODE objects (NOT networkx node) here.
    # edit nodes according to layer specification
    # note if you have compound/cluster which belongs to two layers (metabolic pathway A and B), you have to DUPLICATE node entity in layer.
    ######################################################
    logger.debug('Start locate_nodes_to_layers_and_update_edges()')
    fo = open("log_locate_nodes_to_layers.txt", "w")
    
    for i, cluster_info in dic_cluster_total_input_idx_vs_cluster_info.items():
        if cluster_info['global_accession'] == '1871|MSMS-Pos-Respect_reference.msp':
            print(1)

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
    print(1)
    # take over
    list_edge_info = list_edge_info_mod

    fo.write(f'\n finishing function dic_cluster_total_input_idx_MOD_vs_node_info'
             f'{len(dic_cluster_total_input_idx_MOD_vs_node_info)}')
    logger.info(f'Length of dic_cluster_total_input_idx_MOD_vs_node_info: '
                f'{len(dic_cluster_total_input_idx_MOD_vs_node_info)}')
    logger.debug('finishing Recreating edges with layer info')
    return dic_cluster_total_input_idx_MOD_vs_node_info, list_edge_info, dic_cluster_id_vs_l_cluster_total_input_idx_MOD


#############################
#  [P] Create a network using networkX
#############################
def make_list_of_edge_for_networkx(dic_config, list_edge_info):
    # print "[G] Create a network using networkX "

    list_of_edge_for_networkx = []
    ##########
    for edge_o in list_edge_info:

        if edge_o["edge_type"] in ["inner_ref_layer", "inner_sample_layer", "inter_sample_ref_layer",
                                   "inter_sample_layer"]:
            list_for_one_node = [0] * 3
            # make list corresponds the info of one node.
            list_for_one_node[0] = edge_o["spec_cluster_x_total_input_idx"]
            list_for_one_node[1] = edge_o["spec_cluster_y_total_input_idx"]

            dict_info = {'spec_sim_score': edge_o["spec_sim_score"], 'delta_mz': edge_o["delta_mz"],
                         "edge_type": edge_o["edge_type"]}
            list_for_one_node[2] = dict_info

            list_of_edge_for_networkx.append(list_for_one_node)

    return list_of_edge_for_networkx


def remove_node_edge_with_no_layer_attribute(config_o, dic_cluster_total_input_idx_MOD_vs_node_info,
                                             list_of_edge_for_networkx):
    #####################
    #  [R] removing node/edge with no-layer attribute.
    #   Practically it removes node/layer associated with node whose chemical/toxical
    ################

    #######################
    # if layer attribute is none, remove from further process

    dic_cluster_total_input_idx_MOD_vs_node_info_new = {}

    for total_input_idx_MOD, node_info in dic_cluster_total_input_idx_MOD_vs_node_info.items():

        if len(total_input_idx_MOD.split("|")) > 0:
            node_layer_attri = total_input_idx_MOD.split("|")[1]

        if node_layer_attri != "none":
            dic_cluster_total_input_idx_MOD_vs_node_info_new[total_input_idx_MOD] = node_info

    # dic_cluster_total_input_idx_MOD_vs_node_info = dic_cluster_total_input_idx_MOD_vs_node_info_new

    ####################
    #  if layer attribute is "none", remove from further process

    list_of_edge_for_networkx_new = []

    for edge in list_of_edge_for_networkx:

        node1_layer_attri = ""
        node2_layer_attri = ""

        if len(edge[0].split("|")) > 0:
            node1_layer_attri = edge[0].split("|")[1]
        if len(edge[1].split("|")) > 0:
            node2_layer_attri = edge[1].split("|")[1]

        flag_edge_to_remove = 0
        if node1_layer_attri == "none" or node2_layer_attri == "none":
            flag_edge_to_remove = 1

        if flag_edge_to_remove == 0:
            list_of_edge_for_networkx_new.append(edge)

    return dic_cluster_total_input_idx_MOD_vs_node_info_new, list_of_edge_for_networkx_new


def threshold_edges(dic_config, list_of_edge_for_networkx):
    ###########################
    #  [T] note. now you used threshold to select edges.
    #  which means,  some of nodes in node list is not used in edge list.
    #   This will cause trouble later.
    #  so now you can select nodes actually present in edge list and update.

    # this list contains edge as tuple actually pass threshol
    l_edges_in_use = []
    ##
    #  this list contains node id actually used in netwrok passing threshold.
    list_cluster_global_accession_in_use = []
    logger.debug("function [threshold_edges] select edges using threshold")
    # fo_log.write( "\nlen of list_of_edge_for_networkx before thresholding\n" + str(len(list_of_edge_for_networkx)) + "\n" )
    # this list contain edge color infor in accordance with l_edges_in_use
    l_edges_in_use_color = []

    list_of_edge_for_networkx_new = []
    list_node_total_input_idx_mod_in_use = []

    mz_tol = 0.01

    # dic_cluster_total_input_idx_vs_cluster_info
    for edge in list_of_edge_for_networkx:

        if edge[2]['spec_sim_score'] > dic_config["score_threshold"]:
            list_of_edge_for_networkx_new.append(edge)
            list_node_total_input_idx_mod_in_use.append(edge[0])
            list_node_total_input_idx_mod_in_use.append(edge[1])
            list_node_total_input_idx_mod_in_use = list(set(list_node_total_input_idx_mod_in_use))
            l_edges_in_use.append((edge[0], edge[1]))
            # now check mz difference and if mz match, specify color.
            if abs(edge[2]['delta_mz']) < mz_tol:
                l_edges_in_use_color.append(1)
            else:
                l_edges_in_use_color.append(0)

    # list_of_edge_for_networkx = list_of_edge_for_networkx_new

    logger.debug("Finished function [threshold_edges] ")
    return list_of_edge_for_networkx_new, list_node_total_input_idx_mod_in_use


def create_networkx_graph(config_o, dic_cluster_total_input_idx_MOD_vs_node_info, list_of_edge_for_networkx):
    """
    :param conf_o:
    :param dic_cluster_total_input_idx_MOD_vs_node_info:
    :param list_edge_info:
    :return:
    ######################3
    # [P] Create a network using networkX
    #############################
    print " Create a network using networkX "

    list_of_edge_for_networkx = []
    ##########
    for edge_o in list_edge_info:

        if edge_o.edge_type in ["inner_ref_layer", "inter_sample_ref_layer","inner_sample_layer"]:

            list_for_one_node = [0] * 3
            # make list corresponds the info of one node.
            list_for_one_node[0] = edge_o.spec_cluster_x_total_input_idx
            list_for_one_node[1] = edge_o.spec_cluster_y_total_input_idx

            dict_info = {'spec_sim_score': edge_o.spec_sim_score, 'delta_mz': edge_o.delta_mz}
            list_for_one_node[2] = dict_info

            list_of_edge_for_networkx.append(list_for_one_node)

    """

    ######
    # [U] create graph

    FG = nx.Graph()
    FG.add_edges_from(list_of_edge_for_networkx)
    # fo_log.write( "len of list_of_edge_for_networkx" + str(len(list_of_edge_for_networkx)) )

    #############
    # add node
    #############

    for cluster_total_input_idx_MOD, node_info in dic_cluster_total_input_idx_MOD_vs_node_info.items():
        FG.add_node(cluster_total_input_idx_MOD)

    #####
    # create list od node id here
    ################

    list_total_input_idx_MOD_used = [node[0] for node in FG.nodes.data()]

    # making list of nodes

    list_attribute_for_layer = []

    #########################
    #  iterate list of node and add node with attribute
    #################################
    str_o = ""
    # for cl_o in list_cluster_info:

    for total_input_idx_mod, node_info in dic_cluster_total_input_idx_MOD_vs_node_info.items():

        # this "attribute_for_layer "  keeps tag used for layer separation

        # if ( node_o.spec_cluster.global_accession  in list_cluster_global_accession_to_pass_threshold   ):

        # here you are duplicating node IF thre is more than 0ne layer attribute is possible,
        # for example this node compound is A and B metabo pathway.
        # print node_o

        # add node only if the node is present in edge list.
        if total_input_idx_mod in list_total_input_idx_MOD_used:
            FG.add_node(total_input_idx_mod, attribute_for_layer=node_info["layer"],
                        global_accession=node_info["spec_cluster"]["global_accession"])

            list_attribute_for_layer.append(node_info["layer"])
            list_attribute_for_layer = list(set(list_attribute_for_layer))

    logger.debug(f"list_attribute_for_layer {list_attribute_for_layer}")
    logger.debug("finishing  [create_networkx_graph]")

    return FG


def define_layers(dic_config, dic_cluster_total_input_idx_MOD_vs_node_info):
    logger.debug("starting [define_layers]")
    list_attribute_for_layer = []

    #########################
    #  iterate list of node and add node with attribute
    #################################

    for total_input_idx_mod, node_info in dic_cluster_total_input_idx_MOD_vs_node_info.items():
        list_attribute_for_layer.append(node_info["layer"])
        list_attribute_for_layer = list(set(list_attribute_for_layer))

    logger.debug(f"list_attribute_for_layer {list_attribute_for_layer}")

    ######################################################
    # [V]Now  list_attribute_for_layer  contains all attribute for LAYER SEPARATION.
    #  define layer based on the attributes

    # but first you want to define base layer
    # Tell which attribute is used to specify base layer.
    #  if  your base layer spectra has not CMPD_CLASSIFICATION_SUPERCLASS description, use    str_key_attribute_to_base_layer = '[]'

    # you want to remove base layer key from list.

    if (dic_config["str_key_attribute_to_base_layer"] in list_attribute_for_layer):
        list_attribute_for_layer.remove(dic_config["str_key_attribute_to_base_layer"])

    # create dictionary where key is str of attribute ['organic compound'] and value is layer_id
    dic_attribute_for_layer_vs_layer_id = {}

    #####!!!!!!!!!!!!!!!!   here you are making "sample"layer....
    dic_attribute_for_layer_vs_layer_id[dic_config["str_key_attribute_to_base_layer"]] = 0

    for n in range(0, len(list_attribute_for_layer)):
        # it needs to be n+1 since 0 is for base layer you already created
        dic_attribute_for_layer_vs_layer_id[list_attribute_for_layer[n]] = n + 1

    # make switched version
    dic_layer_id_vs_attribute_for_layer = {}

    for k, v in dic_attribute_for_layer_vs_layer_id.items():
        dic_layer_id_vs_attribute_for_layer[v] = k

    logger.warning(f'dic_layer_id_vs_attribute_for_layer: {dic_layer_id_vs_attribute_for_layer}')

    # for node in FG.nodes(data= True):
    #    node[1]['layer_id'] = dic_attribute_for_layer_vs_layer_id[     node[1]['attribute_for_layer'] ]

    logger.debug(f"list_attribute_for_layer {list_attribute_for_layer}")

    logger.debug("finishing [define_layer]")
    # return FG, dic_layer_id_vs_attribute_for_layer
    return dic_layer_id_vs_attribute_for_layer


"""

def layout_locate_node_to_layers( conf, FG ,layt_2d ,  dic_layer_id_vs_attribute_for_layer,list_node_total_input_idx_mod_in_use \
                                  , dic_cluster_total_input_idx_MOD_vs_node_info,dic_global_accession_vs_mass_feature,l_edges_in_use  ):

    fo = open("layout_locate_node_to_layers.txt" , "w")

    for k, v in dic_cluster_total_input_idx_MOD_vs_node_info.items():
        print("spec_cluster.list_compound_categories" , v.spec_cluster.list_compound_categories)




    ##########################################################
    # [Y1] split nodes to layers according to attribute
    # make layer id vs node id, like   layer id 1 contain node [1,2,3,4,5]

    # but first, make list of layer id present
    l_layer_id = []
    # for node_id_vs_dic in l_tup_node_id_vs_dic:
    fo_node = open("node_info.txt" , "w")

    fo_node.write( "NODE[0]  NODE[1][attribute_for_layer]   NODE[1][layerid] NODE[1][globalaccession]\n")
    for node in FG.nodes(data=True):
        l_layer_id.append(node[1]["layer_id"])
        fo_node.write( node[0]+ "\t"+ str(node[1]["attribute_for_layer"]) +"\t"+ str(node[1]["layer_id"]) +"\t"+ str(node[1]["global_accession"]) +"\n" )

    l_layer_id = list(set(l_layer_id))

    print("layer ids:", l_layer_id)

    #############################################################
    # then make list of tuple ( layer_id , list_of_node_id)
    l_tup_layer_id_vs_l_node_id = []

    for layer_id in l_layer_id:
        l_node_id_x = []
        # for tup_node_id_vs_dic in l_tup_node_id_vs_dic:
        for node in FG.nodes(data=True):
            if node[1]["layer_id"] == layer_id:
                l_node_id_x.append(node[0])
        l_tup_layer_id_vs_l_node_id.append((layer_id, l_node_id_x))

    print("l_tup_layer_id_vs_l_node_id", l_tup_layer_id_vs_l_node_id)

    fo_layer = open("layer_info.txt" , "w")

    for tup in l_tup_layer_id_vs_l_node_id:
        so = str(tup[0]) + ":" + str (tup[1]) +"\n"
        fo_layer.write( so + "\n")

    #####################################################################
    # This list contains layout, one layout is for one layer/mesh
    dic_layer_id_vs_layt = {}
    ########################################################################
    # iterate over  list  of tuple ( layer id , list node id of this layer )
    for tup_layer_id_vs_l_node_id in l_tup_layer_id_vs_l_node_id:
        curr_layer_id = tup_layer_id_vs_l_node_id[0]
        layt_for_layer = {}

        # layt_2d is dictionary where key is node id and value is x,y coordinate as a list  1: [ 2,3]
        for node_id, node_coor in layt_2d.items():
            # if the node id is found in the list of node of the layer
            if node_id in tup_layer_id_vs_l_node_id[1]:
                layt_for_layer[node_id] = node_coor

        dic_layer_id_vs_layt[curr_layer_id] = layt_for_layer



    l_layt = []

    str_o = ""
    for k , v in dic_layer_id_vs_layt.items():
        str_o = str_o +   "layer_id:" +  str(k) + " layt: " + str(v) + "\n"
    #fo_pl.write(   " \n\n\ndic_layer_id_vs_layt\n"  + str_o    )
    #fo_pl.write(" \n\n\nl_layer_id: " + str(l_layer_id))




    #################################################################
    #  [Y3] HERE you have to make coordinates for 3d mesh
    # layer (3d mesh) to put rescaled nodes and edges
    #####
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ###########
    ##  This is causeing weird behavior and errors.......   originally x3,y2
    num_x_grid =  4
    num_y_grid =  4


    #####
    #  making number of z-depth causing some issues
    # now you have to realize number of z depth
    import math
    num_z_depth = int(math.ceil(float( len(l_layer_id)) / float(num_x_grid * num_y_grid))) + 1

    print("l_layer_id:" ,l_layer_id)
    print("len(l_layer_id))" , len(l_layer_id))




    l_my_ranges_x = m3d_mesh.fractionate_range([-10, 10], num_x_grid, ratio_gap=0.02)
    l_my_ranges_y = m3d_mesh.fractionate_range([-10, 10], num_y_grid, ratio_gap=0.02)


    l_my_z = m3d_mesh.fractionate_z_flat([1, 5], num_z_depth)
    l_my_z = sorted(l_my_z, reverse=True)
    # just chech the intermediate



    ############################
    # [Y6]make layers (3d meshes)
    ################################
    # list of graph object of 3d mesh
    list_go_mesh3d = []
    # make 3d mesh coordinates ad dic where key is layer id.

    # making new layer id list which contain only upper )(reference probably


    dic_layer_id_vs_3d_mesh_coordintates = m3d_mesh.make_dic_of_zflat_3d_mesh_coordinates(l_layer_id, l_my_ranges_x,
                                                                                    l_my_ranges_y, l_my_z)

    ###############################################
    # create 3d meshes
    ############################################
    list_go_mesh3d = []

    # define color pallet
    l_colors = ["#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#9467bd", "#8c564b"]

    # dictionary where key is layer group id and value is list of layer id.
    # if you dont want to set any layer, make this dic empty
    #dic_layer_group_id_vs_l_layer_id = {1: [2, 4, 6, 12], 2: [8, 10, 19], 3: [11, 16, 17]}



    ########################
    #  mesh annotation
    ##########################
    print("dic_layer_id_vs_3d_mesh_coordintates",dic_layer_id_vs_3d_mesh_coordintates)

    l_mesh_annotations = []

    for key_id, co in dic_layer_id_vs_3d_mesh_coordintates.items():
        txt_annotation = str(dic_layer_id_vs_attribute_for_layer[key_id])
        # remove sample text from ref layer annotation
        if dic_layer_id_vs_attribute_for_layer[key_id] == "sample":
            txt_annotation = "none"
        l_mesh_annotations.append(dict(showarrow=False, x=co[0][0], y=co[1][0], z=co[2][0],
                                       text= txt_annotation,
                                       xanchor="center",
                                       xshift=10,
                                       yanchor="middle",
                                       opacity=0.9,
                                        font = dict( color="blue", size=20 )
                                       )
                                  )

    ###########################
    #   [Y8]layer info (range) for rescaling
    #
    #   dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat will have
    #      layer_id as key
    #      [x_min, x_max, y_min, y_max, z_Stat]  as value.
    #       This will be used to rescale node/edge coodinate for each layer
    ##############################
    dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat = {}
    for k , v in dic_layer_id_vs_3d_mesh_coordintates.items():
        # v is 3d_mesh_coordintates  [[0.2, 9.8, 0.2, 9.8], [-9.8, -9.8, -0.2, -0.2], [1.0, 1.0, 1.0, 1.0]]
        dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat[ k ] = [   min( v[0]) , max( v[0])  ,   min( v[1]) , max( v[1])    , min( v[2]) ]




    f_make_layer_id_0_as_base = 1

    dic_layer_group_id_vs_l_layer_id = {}

    for layer_id, co in dic_layer_id_vs_3d_mesh_coordintates.items():

        color_code = "#1f77b4"
        # assign layer coloar according to layer group.
        for layer_group_id, l_layer_id in dic_layer_group_id_vs_l_layer_id.items():
            if layer_id in l_layer_id:
                color_code = l_colors[layer_group_id]

        # if you want to make sample layer (id:0) to base
        if f_make_layer_id_0_as_base == 1 :
            # sample layer is 0
            if layer_id != 0 :
                list_go_mesh3d.append( \
                    go.Mesh3d(x=co[0], y=co[1], z=co[2], name = "3dmesh" , opacity=0.20, color=color_code, text="layer_id_" + str(layer_id)))
        # if you do NOT want to make sample layer to base
        if f_make_layer_id_0_as_base != 1:
            list_go_mesh3d.append(
                go.Mesh3d(x=co[0], y=co[1], z=co[2],  name = "3dmesh" , opacity=0.20, color=color_code, text="layer_id_" + str(layer_id)))



    print(" len of list_go_mesh3d " , len(list_go_mesh3d))


    #################################
    ## [Y9] base layer

    if f_make_layer_id_0_as_base == 1 :
        list_go_mesh3d.append(  go.Mesh3d(x=[ -10,10,-10,10], y=[-10,-10,10,10], z=[-3,-3,-3,-3], opacity=0.10, color="green", text="base_" ))


    def merge_two_dicts(x, y):
        z = x.copy()  # start with x's keys and values
        z.update(y)  # modifies z with y's keys and values & returns None
        return z


    ##########################3
    # rescaling  coordinate
    ##############################
    l_base_layer_id = [0]
    layt_arranged = {}


    # note layt_x is layt which is dictionary (   node id is key, vale is [ x_coordinate, y_coordinate]    )
    #  rearrange_networkx_2d_layout_make_3d_zstat   layt_ori, x_min, x_max, y_min, y_max, z_stat


    for layer_id, layt_x in dic_layer_id_vs_layt.items():
        #print "layer_id" ,layer_id
        #print dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat[layer_id]

        # rescaling coordinates for base layer
        # if current layt is for lyer_id 0 which is base layer,  AND you want to make it as base layer
        if f_make_layer_id_0_as_base == 1 :
                ###
                if layer_id == 0 :
                    #print "making special base layer"

                    layt_new = m3d_rescale.rearrange_networkx_2d_layout_make_3d_zstat(layt_x,-10,10,-10,10,-3)
                    layt_arranged = merge_two_dicts(layt_arranged, layt_new)

                if layer_id != 0:
                    layt_new = m3d_rescale.rearrange_networkx_2d_layout_make_3d_zstat(layt_x, \
                            dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat[layer_id][0],
                             dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat[layer_id][1],
                             dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat[layer_id][2],
                              dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat[layer_id][3],
                              dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat[layer_id][4])

                    layt_arranged = merge_two_dicts(layt_arranged, layt_new)


        # if you dont have to make base layer
        if f_make_layer_id_0_as_base == 0:
            layt_new = m3d_rescale.rearrange_networkx_2d_layout_make_3d_zstat(layt_x, \
                                dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat[layer_id][0],
                                dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat[layer_id][1],
                                dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat[layer_id][2],
                                dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat[layer_id][3],
                                dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat[layer_id][4])

            layt_arranged = merge_two_dicts(layt_arranged, layt_new)

    print("++++++++++")
    print(layt_arranged)

    #fo_pl.write(  "\n\n  Layout arranged")
    #for id , coor in layt_arranged.iteritems():
    #    fo_pl.write("\n" + str(id)  + " :" +   str(coor))


    #fo_pl.write( "\nl_cluster_total_input_idx" + str(l_cluster_total_input_idx))
    #fo_pl.write ( "\nl_edges" + str(l_edges))

    #########
    #  MAKE 3D TRACE/COORDINATE
    ####
    # if you get error here, check all nodes area actually in edges


    ########################################
    # ADDING COLOR INFO TO NETWROK
    ########################################

    # note the element of l_nodes_in_use is cluster id.
    import math

    stat_val_threshold = 0.1

    math_log_base = 2
    l_node_color = []
    str_f = "cluster_id\tfeature_id\tval"

    fo_quant = open("quant_process.tsv", "w")

    so = "CLUSTER_ID\tTITLE(FEATURE_ID\tQUANT\n"

    count_hit_feature_list = 0
    val_lowest = None

    for total_input_idx_mod in list_node_total_input_idx_mod_in_use:
        global_accession = dic_cluster_total_input_idx_MOD_vs_node_info[total_input_idx_mod].spec_cluster.global_accession
        feature_val = 1.0

        for k, mf in dic_global_accession_vs_mass_feature.items():
            # if feature_id in dic_global_accession_vs_mass_feature:
            if mf.global_accession == global_accession:
                count_hit_feature_list += 1

                # if statistical significance val (like-ttest-pval) is lower than threshold
                if mf.stat_val < stat_val_threshold:
                    feature_val = float(mf.value_to_show)

        val_processed = math.log(feature_val + 0.001, math_log_base)
        l_node_color.append(val_processed)


    print("l_node_color", l_node_color)
    sys.exit()
    ################
    # color for upper layers for compound information
    for n in range(len(list_node_total_input_idx_mod_in_use)):
        print(n, end=' ')

        fo.write( "\n" +  str(dic_cluster_total_input_idx_MOD_vs_node_info[list_node_total_input_idx_mod_in_use[n]].spec_cluster.list_compound_categories) )
        if     len(dic_cluster_total_input_idx_MOD_vs_node_info[list_node_total_input_idx_mod_in_use[n]].spec_cluster.list_compound_categories) > 0:
            l_node_color[n]= max(l_node_color)
            fo.write(   "\n tox found" )
            print("tox")
    fo.flush()


    ###############
    # define color for edges connecting same mass nodes
    mz_tol = 0.1
    l_edges_in_use_color = []

    for edge in FG.edges.data():
        str_o = ""
        if abs(edge[2]['delta_mz']) < mz_tol:
            l_edges_in_use_color.append(1)
            str_o = "\n" + str(edge) +":" + "samemass"
        else:
            l_edges_in_use_color.append(0)
            str_o = "\n" + str(edge) + ":" + "not same"
        fo.write( str_o)
    fo.flush()

    for node in FG.nodes(data=True):
        l_layer_id.append(node[1]["layer_id"])
        fo_node.write( node[0]+ "\t"+ str(node[1]["attribute_for_layer"]) +"\t"+ str(node[1]["layer_id"]) +"\t"+ str(node[1]["global_accession"]) +"\n" )





    if count_hit_feature_list < 2:
        print(" most of fragment spec info cannot be found in the specified feature table.\n Probably wrong file name for feature table?")
        print("note file name (except file extension) of fragment spec and feature table  has to be exactly  same")
    fo_quant.write(so)
    fo_quant.close()

    str_f = str_f


    print("l_node_color:", l_node_color)

    print("m3d_rescale.get_multi_traces_3d_network_x started")
    print("\n\n")
    traces = m3d_rescale.get_multi_traces_3d_network_x(layt_arranged, list_node_total_input_idx_mod_in_use,
                                                       l_edges_in_use,
                                                       dic_total_input_idx_mod_vs_node_obj=dic_cluster_total_input_idx_MOD_vs_node_info
                                                       , l_node_color=l_node_color, l_edges_color=l_edges_in_use_color)
    print("m3d_rescale.get_multi_traces_3d_network_x finished")
    return traces , list_go_mesh3d,  l_mesh_annotations


"""


def export_current_edges(conf, list_dic_edges_nodes_graph_by_layer, \
                         list_of_edge_for_networkx_to_show_inter_sample_ref_layer,
                         list_of_edge_for_networkx_to_show_inter_sample_layer):
    fo_edge_info = open(os.path.join(os.getcwd(), 'logs', "current_edge_info.tsv"), "w", encoding="utf-8")

    ########################
    # write edge info
    ##########################

    # note CLUSTERID1 ancutally means CLUSTERIDX  "X"

    fo_edge_info.write("X_TOTAL_INPUT_IDX_MOD" + "\t" + "Y_TOTAL_INPUT_IDX_MOD" \
                       + "\t" + "SPEC_SIM_SCORE" \
                       + "\t" + "DELTA_MZ" \
                       + "\t" + "EDGE_TYPE" \
                       + "\n")

    for dicx in list_dic_edges_nodes_graph_by_layer:

        for edge in dicx["list_of_edge_for_networkx"]:
            str_write = str(edge[0]) + "\t" + str(edge[1]) \
                        + "\t" + str(edge[2]["delta_mz"]) \
                        + "\t" + str(edge[2]["spec_sim_score"]) \
                        + "\t" + str(edge[2]["edge_type"]) \
                        + "\n"
            fo_edge_info.write(str_write)

    for edge in list_of_edge_for_networkx_to_show_inter_sample_ref_layer:
        str_write = str(edge[0]) + "\t" + str(edge[1]) \
                    + "\t" + str(edge[2]["delta_mz"]) \
                    + "\t" + str(edge[2]["spec_sim_score"]) \
                    + "\t" + str(edge[2]["edge_type"]) \
                    + "\n"
        fo_edge_info.write(str_write)

    for edge in list_of_edge_for_networkx_to_show_inter_sample_layer:
        str_write = str(edge[0]) + "\t" + str(edge[1]) \
                    + "\t" + str(edge[2]["delta_mz"]) \
                    + "\t" + str(edge[2]["spec_sim_score"]) \
                    + "\t" + str(edge[2]["edge_type"]) \
                    + "\n"
        fo_edge_info.write(str_write)

    fo_edge_info.close()


def export_current_nodes(conf, list_dic_edges_nodes_graph_by_layer):
    fo_node_info = open(os.path.join(os.getcwd(), 'logs', "current_node_info.tsv"), "w", encoding="utf-8")

    str_fo = "TOTAL_INPUT_IDX_MOD" \
             + "\t" + "LAYER" \
             + "\t" + "CLISTER_ID" \
             + "\t" + "GLOBAL_ACCESSION" \
             + "\t" + "DATASET" + "\t" + "DS_SPLIT_KEY" + "\t" + "INCHI" + "\t" + "INCHI_KEY" + "\t" + "ACCESSION_NUMBER" + "\t" + "CMPD_NAME" \
             + "\t" + "FILE_NAME" \
             + "\t" + "PATHWAY_UNIQUE_ID" + "\t" + "PATHWAY_COMMON_NAME" + "\t" + "PATHWAY_COMMON_NAMEX" \
             + "\t" + "RETENTION_TIME_IN_SEC" \
             + "\t" + "PRECURSOR_MZ" \
             + "\t" + "UNIQUE_LEVEL" \
             + "\t" + "PEAK_LIST_MZ_INT_REL" \
             + "\t" + "PEAK_LIST_MZ_ANNOID" \
             + "\t" + "DIC_ANNOID_STRUCT" \
             + "\t" + "NUMBER_OF_SPECTRA" \
             + "\t" + "LIST_ACCESSION_NUMBERS" \
             + "\t" + "CMPD_CLASSIFICATION_KINGDOM" + "\t" + "CMPD_CLASSIFICATION_CLASS" \
             + "\t" + "CMPD_CLASSIFICATION_SUPERCLASS" + "\t" + "CMPD_CLASSIFICATION_ALTERNATIVEPARENT" \
             + "\n"
    fo_node_info.write(str_fo)
    # write cluster attribute

    for dicx in list_dic_edges_nodes_graph_by_layer:

        for total_input_idx_MOD, node_info in dicx["dic_cluster_total_input_idx_MOD_vs_node_info"].items():
            str_fo = \
                str(node_info.total_input_idx_mod) \
                + "\t" + str(node_info.layer) \
                + "\t" + str(node_info.spec_cluster.cluster_id) \
                + "\t" + str(node_info.spec_cluster.global_accession) \
                + "\t" + str(node_info.spec_cluster.inchi) \
                + "\t" + str(node_info.spec_cluster.inchi_key) \
                + "\t" + (node_info.spec_cluster.represen_spec_uni.accession_number) \
                + "\t" + (node_info.spec_cluster.represen_spec_uni.name) \
                + "\t" + str(node_info.spec_cluster.represen_spec_uni.source_filename) \
                + "\t" + str(node_info.spec_cluster.list_external_compound_UNIQUE_ID) \
                + "\t" + str(node_info.spec_cluster.list_pathway_UNIQUE_ID) \
                + "\t" + str(node_info.spec_cluster.list_pathway_COMMON_NAME) \
                + "\t" + str(node_info.spec_cluster.represen_spec_uni.retention_time_in_sec) \
                + "\t" + str(node_info.spec_cluster.represen_spec_uni.precursor_mz) \
                + "\t" + str(node_info.spec_cluster.unique_level) \
                + "\t" + str(node_info.spec_cluster.represen_spec_uni.peak_list_mz_int_rel) \
                + "\t" + str(node_info.spec_cluster.represen_spec_uni.peak_list_mz_annoid) \
                + "\t" + str(node_info.spec_cluster.represen_spec_uni.dic_annoid_struct) \
                + "\t" + str(len(node_info.spec_cluster.list_spec_uni)) \
                + "\t" + str(node_info.spec_cluster.represen_spec_uni.list_cmpd_classification_kingdom) \
                + "\t" + str(node_info.spec_cluster.represen_spec_uni.list_cmpd_classification_class) \
                + "\t" + str(node_info.spec_cluster.represen_spec_uni.list_cmpd_classification_superclass) \
                + "\t" + str(node_info.spec_cluster.represen_spec_uni.list_cmpd_classification_alternative_parent) \
                + "\n"
            fo_node_info.write(str_fo)
    fo_node_info.close()


# end-------- export_current_nodes


"""

def remove_nodes_from_dataset(dic_dataset,config_o):


    list_dic_edges_nodes_graph_by_layer = dic_dataset["list_dic_edges_nodes_graph_by_layer"]
    dic_global_accession_vs_mass_feature = dic_dataset["dic_global_accession_vs_mass_feature"]

    dic_cluster_total_input_idx_MOD_vs_node_info = dic_dataset["dic_cluster_total_input_idx_MOD_vs_node_info"]
    list_of_edge_for_networkx_to_show_inter_sample_ref_layer = dic_dataset["list_of_edge_for_networkx_to_show_inter_sample_ref_layer"] =
    list_of_edge_for_networkx_to_show_inter_sample_layer = dic_dataset["list_of_edge_for_networkx_to_show_inter_sample_layer"]

    # remove node for dic_dataset["list_dic_edges_nodes_graph_by_layer"]
    for dic_edges_nodes_graph_by_layer in dic_dataset["list_dic_edges_nodes_graph_by_layer"] :
        list_of_edge_for_networkx_new = []
        for edge in dic_edges_nodes_graph_by_layer["list_of_edge_for_networkx"]:

            for total_input_ddx_MOD in config_o.l_total_input_idx_MOD_to_remove:
                if edge[0] !=   total_input_ddx_MOD   and  edge[1] != total_input_ddx_MOD:
                    list_of_edge_for_networkx_new.append(edge)
        # take over
        dic_edges_nodes_graph_by_layer["list_of_edge_for_networkx"]=list_of_edge_for_networkx_new

        # remove from networkx graph---------------------

        # get nodes in the graph to remove (  get intersection of l of nodes in graph and l of input idx mod)
        l_total_input_idx_NOD_in_G_to_remove = list(set(dic_edges_nodes_graph_by_layer["nx_graph"].nodes) & set(config_o.l_total_input_idx_MOD_to_remove))
        # then remove using nx function
        dic_edges_nodes_graph_by_layer["nx_graph"].remove_nodes_from(l_total_input_idx_NOD_in_G_to_remove)
"""


def create_data_for_3d_visualization(dic_processed_data, conf):
    logger.info(f'Start {sys._getframe().f_code.co_name}()')
    # Original arguments:
    # list_dic_edges_nodes_graph_by_layer, dic_global_accession_vs_mass_feature, \
    # dic_cluster_total_input_idx_MOD_vs_node_info_all, \
    # list_of_edge_for_networkx_to_show_inter_sample_ref_layer, list_of_edge_for_networkx_to_show_inter_sample_layer

    list_dic_edges_nodes_graph_by_layer = dic_processed_data["list_dic_edges_nodes_graph_by_layer"]
    dic_global_accession_vs_mass_feature = dic_processed_data["dic_global_accession_vs_mass_feature"]
    dic_cluster_total_input_idx_MOD_vs_node_info_all = dic_processed_data[
        "dic_cluster_total_input_idx_MOD_vs_node_info"]
    list_of_edge_for_networkx_to_show_inter_sample_ref_layer = dic_processed_data[
        "list_of_edge_for_networkx_to_show_inter_sample_ref_layer"]
    list_of_edge_for_networkx_to_show_inter_sample_layer = dic_processed_data[
        "list_of_edge_for_networkx_to_show_inter_sample_layer"]

    fo = open("log_create_data_for_3d_visualization.txt", "w")

    fo_pillar = open("log_create_traces_for_network_multilayer_pillar.txt", "w")
    fo_log_layt = open("log_create_traces_for_network_layt.txt", "w")
    ####################################
    #  make identical mass feature in samples as "pillar"
    import statistics

    f_do_pillar_clustering = 1

    if f_do_pillar_clustering == 1:
        # use nx for handling pillar
        G_pillar = nx.Graph()
        #################################
        # mz RT match pillar ----------------
        mz_tol_pillar = 0.01
        rt_sec_tol_pillar = 20

        for e in list_of_edge_for_networkx_to_show_inter_sample_layer:
            # e[2]["spec_sim_score"]
            # e[2]["delta_mz"]

            rt1 = dic_cluster_total_input_idx_MOD_vs_node_info_all[e[0]].spec_cluster.represen_spec_uni.retention_time_in_sec
            rt2 = dic_cluster_total_input_idx_MOD_vs_node_info_all[e[1]].spec_cluster.represen_spec_uni.retention_time_in_sec
            ## !!!!!!!!!!!!!! you should take inte account sim score
            if e[2]["spec_sim_score"] > 0.70 and abs(e[2]["delta_mz"]) < 0.1 and abs(rt1 - rt2) < 20.0:
                G_pillar.add_edge(e[0], e[1])

        # get connected ( whose mz, rt are close,  similarity is high) nodes as sets
        l_set_pillar_cluster = list(nx.connected_components(G_pillar))

        # now iterate for all pillar cluster and adjust layout of nodes who belongs to the pillar cluster-----------
        for set_pillar_cluster in l_set_pillar_cluster:
            fo_pillar.write("\nset pillar nodes" + str(set_pillar_cluster))
            logger.info(f'Set pillar nodes: {str(set_pillar_cluster)}')
            l_layt_x = []
            l_layt_y = []
            # first, get x, y coordinate for all belonging nodes,
            for dic in list_dic_edges_nodes_graph_by_layer:
                if dic["attribute_for_layer"].startswith("sample"):
                    for idx_mod in set_pillar_cluster:
                        if idx_mod in dic["layout_2d"]:
                            # x coordinate
                            l_layt_x.append(dic["layout_2d"][idx_mod][0])
                            # y coordinate
                            l_layt_y.append(dic["layout_2d"][idx_mod][1])
            # make madian of coordinate of pillars
            median_layt_x = statistics.median(l_layt_x)
            median_layt_y = statistics.median(l_layt_y)

            # now replace with median
            for dic in list_dic_edges_nodes_graph_by_layer:
                if dic["attribute_for_layer"].startswith("sample"):
                    for idx_mod in set_pillar_cluster:
                        if idx_mod in dic["layout_2d"]:
                            # x coordinate
                            dic["layout_2d"][idx_mod][0] = median_layt_x
                            # y coordinate
                            dic["layout_2d"][idx_mod][1] = median_layt_y
                            fo_pillar.write(
                                "\n_____________:" + idx_mod + str(dic["layout_2d"][idx_mod][0]) + "," + str(dic["layout_2d"][idx_mod][1]))
                            logger.info(f'layout_2d of {idx_mod} x: {dic["layout_2d"][idx_mod][0]}, '
                                        f'y: {dic["layout_2d"][idx_mod][1]}')
            fo_pillar.write("\npillar median x y:" + str(median_layt_x) + "," + str(median_layt_y))
            logger.info(f'Pillar median x: {median_layt_x}, y: {median_layt_y}')

    fo_log_layt_after_pillar = open("log_create_traces_for_network_multilayer_layt_after_pillar.txt", "w")
    # now replace with median
    log_message = 'layout_2d after pillar'
    for dic in list_dic_edges_nodes_graph_by_layer:
        if dic["attribute_for_layer"].startswith("sample"):
            fo_log_layt_after_pillar.write("\n\n" + str(dic["layout_2d"]))
            log_message += f'\n{str(dic["layout_2d"])}'
    fo_log_layt_after_pillar.flush()
    logger.info(log_message)

    #####################
    # layer id vs attribute
    ######################
    l_attribute_for_layer_sample = []
    l_attribute_for_layer_ref = []

    for dic in list_dic_edges_nodes_graph_by_layer:
        if dic["attribute_for_layer"].startswith("sample"):
            l_attribute_for_layer_sample.append(dic["attribute_for_layer"])
        else:
            l_attribute_for_layer_ref.append(dic["attribute_for_layer"])

    dic_attribute_for_layer_vs_layer_id = {}

    # layer id for sample layer
    l_layer_id_sample = []
    for n in range(len(l_attribute_for_layer_sample)):
        my_id = "s" + str(n)
        dic_attribute_for_layer_vs_layer_id[l_attribute_for_layer_sample[n]] = my_id
        l_layer_id_sample.append(my_id)

    # layer id for ref layer
    l_layer_id_ref = []
    for n in range(len(l_attribute_for_layer_ref)):
        my_id = "r" + str(n)
        dic_attribute_for_layer_vs_layer_id[l_attribute_for_layer_ref[n]] = my_id
        l_layer_id_ref.append(my_id)

    logger.warning(f'dic_attribute_for_layer_vs_layer_id: {dic_attribute_for_layer_vs_layer_id}')
    logger.warning(f'l_layer_id_ref: {l_layer_id_ref}')

    # make switched version
    dic_layer_id_vs_attribute_for_layer = {}

    for k, v in dic_attribute_for_layer_vs_layer_id.items():
        dic_layer_id_vs_attribute_for_layer[v] = k

    ##########################################################
    # [Y1] split nodes to layers according to attribute
    # make layer id vs node id, like   layer id 1 contain node [1,2,3,4,5]

    # but first, make list of layer id present
    # l_layer_id = []
    # for node_id_vs_dic in l_tup_node_id_vs_dic:
    fo_node = open("node_info.txt", "w")

    fo_node.write("NODE[0]  NODE[1][attribute_for_layer]   NODE[1][layerid] NODE[1][globalaccession]\n")

    #################################################################
    #  [Y3] HERE you have to make GRID coordinates for 3d mesh
    # layer (3d mesh) to put rescaled nodes and edges
    #####

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ###########
    ##  This is causeing weird behavior and errors.......   originally x3,y2
    num_x_grid = 4
    num_y_grid = 4

    #####
    # now you have to realize number of z depth for REF
    import math
    num_z_depth = int(math.ceil(float(len(l_layer_id_ref)) / float(num_x_grid * num_y_grid))) + 1

    logger.debug(f"len(l_layer_id_ref )) {len(l_layer_id_ref)}")

    l_my_ranges_x = m3d_mesh.fractionate_range([-10, 10], num_x_grid, ratio_gap=0.02)
    l_my_ranges_y = m3d_mesh.fractionate_range([-10, 10], num_y_grid, ratio_gap=0.02)

    logger.debug(f"num_z_depth: {num_z_depth}")
    logger.debug(
        "!!!!! If you get error here, its possible that node is to small and having problem with creating layer")
    l_my_z = m3d_mesh.fractionate_z_flat([1, 5], num_z_depth)
    l_my_z = sorted(l_my_z, reverse=True)
    # just chech the intermediate
    logger.debug(f"l_my_ranges_x {l_my_ranges_x}")
    logger.debug(f"l_my_ranges_y {l_my_ranges_y}")
    logger.debug(f"l_my_z {l_my_z}")

    ############################
    # [Y6]make layers (3d meshes)
    ################################
    # list of graph object of 3d mesh
    list_go_mesh3d = []
    # make 3d mesh coordinates ad dic where key is layer id.

    # making new layer id list which contain ONLY UPPER (reference probably)
    #   key : layer_id, value : 3d mesh coordinate ([mesh_coord_x, mesh_coord_y, mesh_coord_z])
    dic_layer_id_vs_3d_mesh_coordintates = m3d_mesh.make_dic_of_zflat_3d_mesh_coordinates(l_layer_id_ref, l_my_ranges_x,
                                                                                          l_my_ranges_y, l_my_z)

    logger.warning(f'dic_layer_id_vs_3d_mesh_coordintates: {dic_layer_id_vs_3d_mesh_coordintates}')

    ###############################################
    # create 3d meshes
    ############################################

    # list_go_mesh3d = []
    # define color pallet
    l_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]
    # dictionary where key is layer group id and value is list of layer id.
    # if you dont want to set any layer, make this dic empty
    # dic_layer_group_id_vs_l_layer_id = {1: [2, 4, 6, 12], 2: [8, 10, 19], 3: [11, 16, 17]}

    ########################
    #  mesh annotation
    ##########################
    logger.debug(f"dic_layer_id_vs_3d_mesh_coordintates {dic_layer_id_vs_3d_mesh_coordintates}")
    l_dic_annotations = []
    l_mesh_annotations = []

    for key_id, co in dic_layer_id_vs_3d_mesh_coordintates.items():
        txt_annotation = str(dic_layer_id_vs_attribute_for_layer[key_id])
        # remove sample text from ref layer annotation
        if dic_layer_id_vs_attribute_for_layer[key_id] == "sample":
            txt_annotation = "none"

        dic_annotation = {}
        dic_annotation["x"] = co[0][0]
        dic_annotation["y"] = co[1][0]
        dic_annotation["z"] = co[2][0]
        dic_annotation["text"] = txt_annotation
        l_dic_annotations.append(dic_annotation)

        l_mesh_annotations.append(dict(showarrow=False, x=co[0][0], y=co[1][0], z=co[2][0],
                                       text=txt_annotation,
                                       xanchor="center",
                                       xshift=10,
                                       yanchor="middle",
                                       opacity=0.9,
                                       font=dict(color="blue", size=20)
                                       )
                                  )

    ###########################
    #   [Y8]layer info (range) for rescaling
    #
    #   dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat will have
    #      layer_id as key
    #      [x_min, x_max, y_min, y_max, z_Stat]  as value.
    #       This will be used to rescale node/edge coodinate for each layer
    ##############################
    # this list holds all data necessary to create 3d mesh ()
    l_dic_mesh3d_data = []

    dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat = {}
    for k, v in dic_layer_id_vs_3d_mesh_coordintates.items():
        # v is 3d_mesh_coordintates  [[0.2, 9.8, 0.2, 9.8], [-9.8, -9.8, -0.2, -0.2], [1.0, 1.0, 1.0, 1.0]]
        dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat[k] = [min(v[0]), max(v[0]), min(v[1]), max(v[1]), min(v[2])]

    # flag  sample layer ( layer id starts with "s") as base.
    f_make_layer_id_sx_as_base = 1

    dic_layer_group_id_vs_l_layer_id = {}

    for layer_id, co in dic_layer_id_vs_3d_mesh_coordintates.items():
        color_code = "#1f77b4"

        # assign layer coloar according to layer group.
        # for layer_group_id, l_layer_id in dic_layer_group_id_vs_l_layer_id.iteritems():
        #    if layer_id in l_layer_id:
        #        color_code = l_colors[layer_group_id]

        # if you want to make sample layer (id:0) to base
        if f_make_layer_id_sx_as_base == 1:
            # dealing with reference layer
            #  sample layer id is always startwith "s" (s1,s2)
            if not layer_id.startswith("s"):
                dic_mesh_data = {}
                dic_mesh_data["l_x"] = co[0]
                dic_mesh_data["l_y"] = co[1]
                dic_mesh_data["l_z"] = co[2]
                dic_mesh_data["color"] = color_code
                dic_mesh_data["text"] = "layer_id_" + str(layer_id)

                l_dic_mesh3d_data.append(dic_mesh_data)

                list_go_mesh3d.append(go.Mesh3d(x=co[0], y=co[1], z=co[2], opacity=0.20, color=color_code,
                                                text="layer_id_" + str(layer_id)))

        # if you do NOT want to make sample layer to base
        if f_make_layer_id_sx_as_base != 1:
            dic_mesh_data = {}
            dic_mesh_data["l_x"] = co[0]
            dic_mesh_data["l_y"] = co[1]
            dic_mesh_data["l_z"] = co[2]
            dic_mesh_data["color"] = color_code
            dic_mesh_data["text"] = "layer_id_" + str(layer_id)
            l_dic_mesh3d_data.append(dic_mesh_data)
            list_go_mesh3d.append(
                go.Mesh3d(x=co[0], y=co[1], z=co[2], opacity=0.20, color=color_code, text="layer_id_" + str(layer_id)))

    logger.debug(f" len of list_go_mesh3d {len(list_go_mesh3d)}")

    #################################
    ## [Y9] making base layer
    ##############################

    dic_layer_id_vs_z_depth = {}
    n_base_depth = 0
    for layer_id_sample in l_layer_id_sample:
        n_base_depth = n_base_depth + 1
        if f_make_layer_id_sx_as_base == 1:
            # change depth of base layerS
            z_val = -2 - n_base_depth * 2

            dic_mesh_data = {}
            dic_mesh_data["l_x"] = [-10, 10, -10, 10]
            dic_mesh_data["l_y"] = [-10, -10, 10, 10]
            dic_mesh_data["l_z"] = [z_val, z_val, z_val, z_val]
            dic_mesh_data["color"] = "green"
            dic_mesh_data["text"] = "layer_id_" + str(layer_id_sample)

            l_dic_mesh3d_data.append(dic_mesh_data)

            list_go_mesh3d.append(
                go.Mesh3d(x=[-10, 10, -10, 10], y=[-10, -10, 10, 10], z=[z_val, z_val, z_val, z_val], opacity=0.1,
                          color="green", text="layer_id_" + str(layer_id_sample)))
            dic_layer_id_vs_z_depth[layer_id_sample] = z_val
    ###########################
    fo.write("n_base_depth:" + str(n_base_depth))
    fo.write("\ndic_layer_id_vs_z_depth (now only sample):" + str(dic_layer_id_vs_z_depth))
    logger.info(f'n_base_depth:: {n_base_depth}')
    logger.info(f'dic_layer_id_vs_z_depth (now only sample): {str(dic_layer_id_vs_z_depth)}')

    def merge_two_dicts(x, y):
        z = x.copy()  # start with x's keys and values
        z.update(y)  # modifies z with y's keys and values & returns None
        return z

    ##########################3
    # [Y10]rescaling  coordinate
    ##############################
    l_base_layer_id = [0]
    layt_arranged = {}

    ### !!!!!!!!!!!!!!!!!!!!!
    #   list_dic_edges_nodes_graph_by_layer
    fo_xxxxxxxx = open("layout_xxxxxx.txt", "w")

    for dic in list_dic_edges_nodes_graph_by_layer:

        if f_make_layer_id_sx_as_base == 1:
            # dic_layer_id_vs_attribute_for_layer
            layer_id = dic_attribute_for_layer_vs_layer_id[dic["attribute_for_layer"]]
            logger.warning(f'layer_id: {layer_id}')

            ### for sample/base layers----------------------
            if layer_id.startswith("s"):
                logger.debug("making BASE layer")

                fo_xxxxxxxx.write("\n\n" + str(dic["layout_2d"]))
                ###############
                # !!!!!!!!!!!!!!!!!!!!!!  here you maing corrdinate shift !!!!!!!!!!!!!!!!!!
                ## for base layer, depth can be different for multiple ase layers
                layout_3d = m3d_rescale.rearrange_networkx_2d_layout_make_3d_zstat(dic["layout_2d"], -10, 10, -10, 10,
                                                                                   dic_layer_id_vs_z_depth[layer_id])
                dic["layout_3d"] = layout_3d

                fo_log_layt.write("\n\n" + str(dic["layout_3d"]))
                logger.info(f'layout_3d of {dic["attribute_for_layer"]}: {layout_3d}')

                fo_log_layt.flush()

            #### make 3d layout only for dataset that has layout 2d.   in other words, no layout 3d is made for empty network.
            if layer_id.startswith("r") and len(dic["layout_2d"]) > 0:
                layout_3d = m3d_rescale.rearrange_networkx_2d_layout_make_3d_zstat(
                    dic["layout_2d"],
                    dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat[layer_id][0],
                    dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat[layer_id][1],
                    dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat[layer_id][2],
                    dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat[layer_id][3],
                    dic_layer_id_vs_list_3d_mesh_min_max_value_xy_zstat[layer_id][4]
                )
                dic["layout_3d"] = layout_3d

            # layt_arranged = merge_two_dicts(layt_arranged, layt_new)

    logger.debug("++++++++++")
    logger.debug(layt_arranged)

    ########################################
    # [Y11]ADDING COLOR INFO TO NETWROK
    ########################################
    logger.debug("ADDING COLOR INFO TO NETWROK")
    # note the element of l_nodes_in_use is cluster id.
    import math

    stat_val_threshold = 1.0
    math_log_base = 2
    l_node_color = []
    str_f = "cluster_id\tfeature_id\tval"
    fo_quant = open("quant_process.tsv", "w")

    so = "CLUSTER_ID\tTITLE(FEATURE_ID\tQUANT\n"
    count_hit_feature_list = 0
    val_lowest = None

    logger.debug("evaluating statistcal value and thresholding")
    # list_node_total_input_idx_mod_in_use = []
    # this roughly collect all color value, and will be used to get max value later.
    l_node_color_val_all = []

    for dic in list_dic_edges_nodes_graph_by_layer:

        dic_cluster_total_input_idx_MOD_vs_node_info = dic["dic_cluster_total_input_idx_MOD_vs_node_info"]

        for cluster_total_input_idx_MOD, node_info in dic_cluster_total_input_idx_MOD_vs_node_info.items():
            spec_cluster = node_info["spec_cluster"]

            if spec_cluster['tag'] == 'sample':
                global_accession = spec_cluster["global_accession"]
                if conf['quant_value_type'] == 'ratio':
                    val_processed = 0.0000
                    for k, mf in dic_global_accession_vs_mass_feature.items():

                        # print ("mf val show", k,  mf.stat_val , mf.value_to_show)
                        # if feature_id in dic_global_accession_vs_mass_feature:
                        if mf['global_accession'] == global_accession:
                            # logger.debug(f"{k}, {mf.stat_val}, {type(mf.stat_val)}, {mf.value_to_show}, {type(mf.value_to_show)}")
                            count_hit_feature_list += 1

                            # if statistical significance val (like-ttest-pval) is lower than threshold
                            if mf['stat_val'] < stat_val_threshold and mf['value_to_show'] != 0:
                                try:
                                    feature_val = abs(float(mf['value_to_show']))  # TODO: feature_val must be > 0
                                    # logger.debug(f"threshold pass {mf.stat_val}, {feature_val}, {float(mf.value_to_show)}")
                                    logger.debug(f'k: {k}, feature_val: {feature_val}')
                                    val_processed = math.log(feature_val, math_log_base)
                                except (ValueError, TypeError):
                                    pass
                    node_info["color"] = val_processed
                    l_node_color_val_all.append(val_processed)
                else:
                    feature_val = 0.0000
                    for k, mf in dic_global_accession_vs_mass_feature.items():
                        if mf['global_accession'] == global_accession:
                            count_hit_feature_list += 1

                            if mf['stat_val'] < stat_val_threshold:
                                try:
                                    feature_val = float(mf['value_to_show'])
                                except (ValueError, TypeError):
                                    pass
                    node_info["color"] = feature_val
                    l_node_color_val_all.append(feature_val)
            else:
                if spec_cluster['extra_node_color']:
                    node_info["color"] = spec_cluster['extra_node_color']
                else:
                    node_info["color"] = 0.0000

    logger.debug(f"l_node_color_val_all {l_node_color_val_all}")

    #############
    #  external data-dependednt color coding is not working
    #  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ################
    # color for upper layers for compound information

    # for dic in list_dic_edges_nodes_graph_by_layer:
    #     for k, node_info in dic["dic_cluster_total_input_idx_MOD_vs_node_info"].items():
    #         if len(node_info["spec_cluster"]["list_compound_categories"]) > 0:
    #             if (conf['type_attribute_for_layer_separation'] == 'list_compound_categories'
    #                     and conf['color_toxic_compound']):  # TODO: Color toxic compound
    #                 node_info['color'] = '#00e08e'
    #             else:
    #                 node_info['color'] = max(l_node_color_val_all)  # TODO: node_info['color'] must be number
    #             # node_info['color'] = "white"

    ###############
    # define color for edges connecting same mass nodes
    mz_tol = 0.1
    rt_in_sec_tol = 30.0
    l_edges_in_use_color = []

    fo_log_xxx = open("log_xxx.txt", "w")

    ###=============
    for dic in list_dic_edges_nodes_graph_by_layer:
        list_of_edge_for_networkx = dic["list_of_edge_for_networkx"]
        fo_log_xxx.write("\n\n")
        fo_log_xxx.write(str(list_of_edge_for_networkx))
        logger.info(f'list_of_edge_for_networkx of {dic["attribute_for_layer"]}: {list_of_edge_for_networkx}')
        for edge in list_of_edge_for_networkx:
            if abs(edge[2]['delta_mz']) < mz_tol:
                edge[2]['color'] = "red"
            else:
                edge[2]['color'] = "black"

    # color for INTER_LAYER network
    for edge in list_of_edge_for_networkx_to_show_inter_sample_ref_layer:
        if abs(edge[2]['delta_mz']) < mz_tol:
            edge[2]['color'] = "purple"
        else:
            edge[2]['color'] = "blue"

    # for inter sample edges----------------------------------------------------
    list_of_edge_for_networkx_to_show_inter_sample_layer_new = []

    if f_do_pillar_clustering == 1:
        # look for edges both nodes are in pillar cluster
        for edge in list_of_edge_for_networkx_to_show_inter_sample_layer:
            for set_pillar_cluster in l_set_pillar_cluster:
                # if nodes in this edge can be found in the current set of pillar cluster.
                if edge[0] in set_pillar_cluster and edge[1] in set_pillar_cluster:
                    edge[2]['color'] = "orange"
                    list_of_edge_for_networkx_to_show_inter_sample_layer_new.append(edge)

        f_show_inter_layer_identical_compound_edges = 1

        if f_show_inter_layer_identical_compound_edges == 1:
            list_of_edge_for_networkx_to_show_inter_sample_layer = list_of_edge_for_networkx_to_show_inter_sample_layer_new

    if f_do_pillar_clustering == 0:
        for edge in list_of_edge_for_networkx_to_show_inter_sample_layer:

            node_info_1 = dic_cluster_total_input_idx_MOD_vs_node_info_all[edge[0]]
            node_info_2 = dic_cluster_total_input_idx_MOD_vs_node_info_all[edge[1]]

            delta_rt_in_sec = abs(
                node_info_1.spec_cluster.represen_spec_uni.retention_time_in_sec - node_info_2.spec_cluster.represen_spec_uni.retention_time_in_sec)

            if abs(edge[2]['delta_mz']) < mz_tol and delta_rt_in_sec < rt_in_sec_tol:
                edge[2]['color'] = "orange"
                list_of_edge_for_networkx_to_show_inter_sample_layer_new.append(edge)
            else:
                edge[2]['color'] = "green"

        # flag to show only inter sample layer edges showing identical compound. (same prec mz and same rt)
        f_show_inter_layer_identical_compound_edges = 0

        if f_show_inter_layer_identical_compound_edges == 1:
            list_of_edge_for_networkx_to_show_inter_sample_layer = list_of_edge_for_networkx_to_show_inter_sample_layer_new

    ## [Y15] additional info -----------------------------------

    ## define edge text --------------------------------

    for dic in list_dic_edges_nodes_graph_by_layer:
        l_edge_text = []
        list_of_edge_for_networkx = dic["list_of_edge_for_networkx"]
        for edge in list_of_edge_for_networkx:
            l_edge_text.append("spec_sim_score:" + str(edge[2]['spec_sim_score']))
        dic["l_edge_text"] = l_edge_text

    ###===============

    if count_hit_feature_list < 2:
        logger.debug(
            " most of fragment spec info cannot be found in the specified feature table.\n Probably wrong file name for feature table?")
        logger.debug(
            "note file name (except file extension) of fragment spec and feature table  has to be exactly  same")
    fo_quant.write(so)
    fo_quant.close()

    logger.debug(f"l_node_color {l_node_color}")

    """
    traces = m3d_rescale.get_multi_traces_3d_network_x(layt_arranged, list_node_total_input_idx_mod_in_use,
                                                       l_edges_in_use,
                                                       dic_total_input_idx_mod_vs_node_obj=dic_cluster_total_input_idx_MOD_vs_node_info
                                                       , l_node_color=l_node_color, l_edges_color=l_edges_in_use_color)
    """

    dic_rescaled_geometry_for_3d_network = m3d_rescale.get_rescaled_geometry_for_3d_network_x(
        list_dic_edges_nodes_graph_by_layer,
        list_of_edge_for_networkx_to_show_inter_sample_ref_layer,
        list_of_edge_for_networkx_to_show_inter_sample_layer)

    # traces = m3d_rescale.get_multi_traces_3d_network_x(  list_dic_edges_nodes_graph_by_layer , \
    #                                                     list_of_edge_for_networkx_to_show_inter_sample_ref_layer ,list_of_edge_for_networkx_to_show_inter_sample_layer)

    log_message = f'Length of each data in dic_rescaled_geometry_for_3d_network["list_dic_edges_nodes_graph_by_layer"]' \
                  f'\nattribute_for_layer\tdic_cluster_total_input_idx_MOD_vs_node_info\tlayout_2d\tlayout_3d'
    for dic in dic_rescaled_geometry_for_3d_network["list_dic_edges_nodes_graph_by_layer"]:
        fo.write(
            "\n" + dic["attribute_for_layer"] + " : " + str(len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])) \
            + " : " + str(len(dic["layout_2d"])) + " : " + str(len(dic["layout_2d"])))

        log_message += f'\n{dic["attribute_for_layer"]}\t{len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])}' \
                       f'\t{len(dic["layout_2d"])}\t{len(dic["layout_2d"])}'

    logger.info(log_message)

    dic_dataset_for_3d_network_visualization = {}

    dic_dataset_for_3d_network_visualization["list_dic_edges_nodes_graph_by_layer"] = \
    dic_rescaled_geometry_for_3d_network["list_dic_edges_nodes_graph_by_layer"]

    dic_dataset_for_3d_network_visualization["dic_edges_inter_sample_ref_layer"] = dic_rescaled_geometry_for_3d_network[
        "dic_edges_inter_sample_ref_layer"]
    dic_dataset_for_3d_network_visualization["dic_edges_inter_sample_layer"] = dic_rescaled_geometry_for_3d_network[
        "dic_edges_inter_sample_layer"]

    dic_dataset_for_3d_network_visualization["l_dic_mesh3d_data"] = l_dic_mesh3d_data
    dic_dataset_for_3d_network_visualization["l_dic_annotations"] = l_dic_annotations
    dic_dataset_for_3d_network_visualization['dic_layer_id_vs_attribute_for_layer'] = dic_layer_id_vs_attribute_for_layer  # TODO: Added temporarily (K. Hirata, 20221214)

    return dic_dataset_for_3d_network_visualization

    # return traces , list_go_mesh3d,  l_mesh_annotations


def read_data_for_multilayer_3d_network(dic_config):
    ###################
    # read config
    logger.info(f'Start {sys._getframe().f_code.co_name}()')
    logger.debug("[C]reading cluster attribute")

    fo_log = open("log_read_data_for_multilayer_3d_network.txt", "w")

    #####################
    # read input data

    ##########################################################################################
    # section D
    # [D] read spec cluster info
    ##########################################################################################
    logger.debug("* main [D]reading cluster attribute")
    dic_cluster_total_input_idx_vs_cluster_info = read_cluster_attribute(dic_config['filename_cluster_info'])

    # reading external info
    ###########################################
    # [F] Read external info for node/cluster
    ###########################################
    logger.warning("* main [D]reading external compound info")
    add_external_cmpd_info(dic_cluster_total_input_idx_vs_cluster_info, dic_config['foldername_ext_cmpd_info'])
    logger.warning("* main [D]finished reading external compound info")

    # Add color to T3DB compounds
    add_color_to_t3db_compound(dic_config['color_toxic_compound'], dic_cluster_total_input_idx_vs_cluster_info)

    # make original copy
    dic_cluster_total_input_idx_vs_cluster_info_original = copy.deepcopy(dic_cluster_total_input_idx_vs_cluster_info)

    logger.warning(f"* main[D]len dic_cluster_total_input_idx_vs_cluster_info "
                   f"{len(dic_cluster_total_input_idx_vs_cluster_info)}")

    ##########################################################################################
    # section [H]
    # [H] read EDGES info
    ##########################################################################################
    l_edges, list_edge_info = read_edge_info(dic_config['filename_edge_info'], dic_config['score_threshold'])

    fo_log.write("\n\n[F1]   list_edge_info_original")
    log_message = f'[F1]   list_edge_info_original' \
                  f'\nspec_cluster_x_global_accession <-> spec_cluster_y_global_accession'
    for e in list_edge_info:
        fo_log.write("\n" + e["spec_cluster_x_global_accession"] + " <-> " + e["spec_cluster_y_global_accession"])
        log_message += f'\n{e["spec_cluster_x_global_accession"]} <-> {e["spec_cluster_y_global_accession"]}'

    logger.info(log_message)

    #############################
    # [J] read feature table
    #############################
    logger.debug("[J] reading feature table")
    #  --------------MODIFIED20220714
    # read_feature table.
    # TODO: Check feature_table_parser
    dic_global_accession_vs_mass_feature_original = read_feature_table(dic_config)
    dic_global_accession_vs_mass_feature = read_feature_table(dic_config)

    # select nodes based on keyword
    dic_cluster_total_input_idx_vs_cluster_info_new =\
        select_nodes_based_on_keyword(dic_config['filter_select_category'],
                                      dic_config['filter_select_keyword'],
                                      dic_cluster_total_input_idx_vs_cluster_info)
    logger.warning(f"select nodes based on keyword : len dic_cluster_total_input_idx_vs_cluster_info_new "
                   f"{len(dic_cluster_total_input_idx_vs_cluster_info_new)}")

    # take over
    dic_cluster_total_input_idx_vs_cluster_info = dic_cluster_total_input_idx_vs_cluster_info_new
    logger.debug(
        f"len dic_cluster_total_input_idx_vs_cluster_info: {str(len(dic_cluster_total_input_idx_vs_cluster_info))}")

    # locate nodes to layers
    dic_cluster_total_input_idx_MOD_vs_node_info, list_edge_info, dic_cluster_id_vs_l_cluster_total_input_idx_MOD = \
        locate_nodes_to_layers_and_update_edges(dic_config, dic_cluster_total_input_idx_vs_cluster_info, list_edge_info)

    logger.debug("returned from  locate_nodes_to_layers_and_update_edges")
    # make list_of_edges_for_networkx

    logger.debug("Main  [make list_of_edges_for_networkx]")
    list_of_edge_for_networkx = make_list_of_edge_for_networkx(dic_config, list_edge_info)
    list_of_edge_for_networkx_new, list_node_total_input_idx_mod_in_use = threshold_edges(dic_config,
                                                                                          list_of_edge_for_networkx)

    dic_source_data = {}
    dic_source_data["dic_config"] = dic_config
    dic_source_data["list_edge_info"] = list_edge_info
    dic_source_data[
        "dic_cluster_total_input_idx_vs_cluster_info_original"] = dic_cluster_total_input_idx_vs_cluster_info_original
    dic_source_data["list_of_edge_for_networkx"] = list_of_edge_for_networkx
    dic_source_data["list_node_total_input_idx_mod_in_use"] = list_node_total_input_idx_mod_in_use
    dic_source_data["dic_cluster_total_input_idx_vs_cluster_info"] = dic_cluster_total_input_idx_vs_cluster_info
    dic_source_data["dic_global_accession_vs_mass_feature_original"] = dic_global_accession_vs_mass_feature_original
    dic_source_data["dic_global_accession_vs_mass_feature"] = dic_global_accession_vs_mass_feature

    fo_log.close()
    return dic_source_data


def process_3d_network_data(dic_source_data, dic_config):
    # dic_source_data_1

    dic_global_accession_vs_mass_feature = dic_source_data["dic_global_accession_vs_mass_feature"]

    logger.debug("[create_multilayer_3d_network_data]processing")

    fo_log = open("log_create_multilayer_3d_network_data.txt", "w")
    fo_log_key = open("log_KEY_multilayer_3d_network.txt", "w")

    # Read some external files
    # suspect mz file

    dic_config["dic_filename_vs_dic_suspect_cmpd_vs_mz"] = suspect_compound.read_file_suspect_mz(dic_config)

    ###   [C1]]--------------------------------------

    # deeocopy input idx and cluster info
    dic_cluster_total_input_idx_vs_cluster_info_original = copy.deepcopy(
        dic_source_data["dic_cluster_total_input_idx_vs_cluster_info"])

    fo_log.write("\n\n\n another fig update")
    fo_log.write("\n dic_cluster_total_input_idx_vs_cluster_info_original:" + str(
        len(dic_cluster_total_input_idx_vs_cluster_info_original)))

    logger.info('Another fig update')
    logger.info(f'Length of dic_cluster_total_input_idx_vs_cluster_info_original: '
                f'{len(dic_cluster_total_input_idx_vs_cluster_info_original)}')

    # Add color to T3DB compounds
    add_color_to_t3db_compound(dic_config, dic_cluster_total_input_idx_vs_cluster_info_original)

    # select nodes based on keyword
    dic_cluster_total_input_idx_vs_cluster_info_new =\
        select_nodes_based_on_keyword(dic_config['filter_select_category'],
                                      dic_config['filter_select_keyword'],
                                      dic_cluster_total_input_idx_vs_cluster_info_original)
    logger.debug(f"\n selecte nodes based on keyword : len ddic_cluster_total_input_idx_vs_cluster_info_new: {str(len(dic_cluster_total_input_idx_vs_cluster_info_new))}")

    #  [C2]
    # select nodes based on prec mz
    dic_cluster_total_input_idx_vs_cluster_info_new = select_nodes_based_on_prec_mz(dic_config,
                                                                                    dic_cluster_total_input_idx_vs_cluster_info_new)

    # take over
    dic_cluster_total_input_idx_vs_cluster_info = dic_cluster_total_input_idx_vs_cluster_info_new

    logger.debug(
        f"[C2]: len dic_cluster_total_input_idx_vs_cluster_info {str(len(dic_cluster_total_input_idx_vs_cluster_info))}")

    fo_log.write(f"\n after selecting nodes, dic_cluster_total_input_idx_vs_cluster_info: "
                 f"{len(dic_cluster_total_input_idx_vs_cluster_info)}")
    fo_log.flush()
    logger.info(f'After selecting nodes, length of dic_cluster_total_input_idx_vs_cluster_info: '
                f'{len(dic_cluster_total_input_idx_vs_cluster_info)}')

    list_edge_info_original = copy.deepcopy(dic_source_data["list_edge_info"])

    fo_log.write("\n before locate_nodes_to_layers_and_update_edges  len of list_edge_info_original:" + str(
        len(list_edge_info_original)))
    fo_log.flush()
    logger.info(f'Before locate_nodes_to_layers_and_update_edges, length of list_edge_info_original: '
                f'{len(list_edge_info_original)}')

    # MODIFIED20220714---------------
    #  [C3]-------------------------------------------
    # select nodes based on mass defect
    if len(dic_config["list_mass_defect"]) > 0:
        dic_cluster_total_input_idx_vs_cluster_info_new = select_nodes_based_on_mass_defect(dic_config,
                                                                                            dic_cluster_total_input_idx_vs_cluster_info)
        # take over
        dic_cluster_total_input_idx_vs_cluster_info = dic_cluster_total_input_idx_vs_cluster_info_new
    # ------------------------------- MODIFIED20220714---------------

    # MODIFIED20220719---------------
    # [C4] select node based on product ion
    if len(dic_config["list_product_mz_required"]) > 0:
        dic_cluster_total_input_idx_vs_cluster_info_new = select_nodes_based_on_product_mz_required(
            dic_config, dic_cluster_total_input_idx_vs_cluster_info)
        # take over
        dic_cluster_total_input_idx_vs_cluster_info = dic_cluster_total_input_idx_vs_cluster_info_new
        logger.debug(f" l of dic after product filer: {len(dic_cluster_total_input_idx_vs_cluster_info)}")

    # ------------------------------MODIFIED20220719---------------

    # remove nodes based on user defined input idx----------------------------------
    logger.debug(f"dic_config[l_total_input_idx_to_remove] {str(dic_config['l_total_input_idx_to_remove'])}")

    dic_cluster_total_input_idx_vs_cluster_info_new = {}
    for cluster_total_input_idx, cluster_info in dic_cluster_total_input_idx_vs_cluster_info.items():
        # [C6]   remove node based on user-defined global accession
        if not cluster_total_input_idx in dic_config["l_total_input_idx_to_remove"]:
            dic_cluster_total_input_idx_vs_cluster_info_new[cluster_total_input_idx] = cluster_info

    # take over
    dic_cluster_total_input_idx_vs_cluster_info = dic_cluster_total_input_idx_vs_cluster_info_new

    ## also remove edges that have nodes to be removed.
    list_edge_info_original_edit = []
    for edge_info in list_edge_info_original:
        f_match_to_remove = 0
        # if one of the nodes in this edge is supposed to be removed
        if edge_info["spec_cluster_x_total_input_idx"] in dic_config["l_total_input_idx_to_remove"] or edge_info in \
                dic_config["l_total_input_idx_to_remove"]:
            f_match_to_remove = 1
        # if there is no match, append
        if f_match_to_remove == 0:
            list_edge_info_original_edit.append(edge_info)

    ###############################33
    #  [F1]
    # update node and edges.  now node id is "MOD"  and make node info
    #     now edge has attribute "inner_sample_layer" etc.

    logger.debug("[F1] :  update node and edges.")

    fo_log.write("\n\n[F1]   list_edge_info_original")
    log_message = f'[F1] list_edge_info_original' \
                  f'\nspec_cluster_x_global_accession <-> spec_cluster_y_global_accession'
    for e in list_edge_info_original:
        fo_log.write("\n" + e["spec_cluster_x_global_accession"] + " <-> " + e["spec_cluster_y_global_accession"])
        log_message += f'\n{e["spec_cluster_x_global_accession"]} <-> {e["spec_cluster_y_global_accession"]}'
    logger.info(log_message)

    fo_log.write("\n\n[F1]   list_edge_info_original_edit")
    log_message = f'[F1] list_edge_info_original_edit' \
                  f'\nspec_cluster_x_global_accession <-> spec_cluster_y_global_accession'
    for e in list_edge_info_original_edit:
        fo_log.write("\n" + e["spec_cluster_x_global_accession"] + " <-> " + e["spec_cluster_y_global_accession"])
        log_message += f'\n{e["spec_cluster_x_global_accession"]} <-> {e["spec_cluster_y_global_accession"]}'
    logger.info(log_message)

    dic_cluster_total_input_idx_MOD_vs_node_info, list_edge_info, dic_cluster_id_vs_l_cluster_total_input_idx_MOD = \
        locate_nodes_to_layers_and_update_edges(dic_config, dic_cluster_total_input_idx_vs_cluster_info,
                                                list_edge_info_original_edit)

    fo_log.write("\n\n" + "dic_cluster_total_input_idx_vs_cluster_info")
    for idx, cl in dic_cluster_total_input_idx_vs_cluster_info.items():
        fo_log.write("\n" + str(idx))
    fo_log.flush()
    logger.info(f'Keys of dic_cluster_total_input_idx_vs_cluster_info: '
                f'{dic_cluster_total_input_idx_vs_cluster_info.keys()}')

    logger.debug("[F1] finished  locate_nodes_to_layers_and_update_edges.")
    fo_log.write("\n after locate_nodes_to_layers_and_update_edges:" + str(len(list_edge_info)))
    fo_log.flush()
    logger.info(f'After locate_nodes_to_layers_and_update_edges, length of list_edge_info: {len(list_edge_info)}')

    # make list_of_edges_for_networkx
    fo_log.write("\n\n[F1] list_edge_info updated")
    log_message = f'[F1] list_edge_info updated' \
                  f'\nspec_cluster_x_global_accession <-> spec_cluster_y_global_accession'
    for e in list_edge_info:
        fo_log.write("\n" + e["spec_cluster_x_global_accession"] + " <-> " + e["spec_cluster_y_global_accession"])
        log_message += f'\n{e["spec_cluster_x_global_accession"]} <-> {e["spec_cluster_y_global_accession"]}'
    logger.info(log_message)

    # what you get is [ node_id_X, node_id_Y, dictionary_for_attribute]
    list_of_edge_for_networkx = make_list_of_edge_for_networkx(dic_config, list_edge_info)
    fo_log.write("\n\n [F1]Show raw edges----------")
    log_message = '[F1] Show raw edges'
    for e in list_of_edge_for_networkx:
        fo_log.write("\n" + str(e))
        log_message += f'{e}'
    logger.info(log_message)

    # MODIFIED20220714-------------------
    # [F4]   add info based on  with suspect list.--------------------------------------------------
    logger.debug("starting suspect list matching")

    fo_log_key.write("\n number of suspect files: " + str(dic_config["dic_filename_vs_dic_suspect_cmpd_vs_mz"].keys()))
    logger.info(f'Number of suspect files: {len(dic_config["dic_filename_vs_dic_suspect_cmpd_vs_mz"].keys())}')
    for k, v in dic_config["dic_filename_vs_dic_suspect_cmpd_vs_mz"].items():
        fo_log_key.write("\n   " + k + " : " + str(len(v)))
        logger.info(f'Number of suspect compounds in {k}: {len(v)}')

    mz_tol_suspect = dic_config["mz_tol"]

    # iterate all { cluuster_total_input_idx, cluster_info}
    for cluster_total_input_idx, node_info in dic_cluster_total_input_idx_MOD_vs_node_info.items():

        # suspect screening only applies to sample layer---
        if node_info["layer"].startswith("sample"):

            ## !!!!!!!!!!!! TO BE FIXED
            ## currently only mz of representeive spec of cluster is matched to suspect list.
            # for each spectrum in the cluster....
            # for sp_uni in node_info.spec_cluster.list_spec_uni:
            # scan all suspect mz ---------------

            # for each suspect file,
            for path, dic_cmpd_vs_mz in dic_config["dic_filename_vs_dic_suspect_cmpd_vs_mz"].items():
                filename = os.path.splitext(os.path.basename(path))[0]
                # for each compound in suspect
                for cmpd, mz_suspect in dic_cmpd_vs_mz.items():
                    # if suspect mz match to observed mass.
                    if abs(node_info["spec_cluster"]["represen_spec_uni"][
                               "precursor_mz"] - mz_suspect) < mz_tol_suspect:
                        node_info["suspect"] = filename
                        node_info["l_source_suspect"].append(filename)
                        node_info["l_source_suspect"] = list(set(node_info["l_source_suspect"]))

                        if not node_info["name"]:
                            node_info["name"] = f'{cmpd}_in({filename})'
                        else:
                            node_info["name"] += f'_{cmpd}_in({filename})'

    # -----------------MODIFIED20220714

    logger.debug("suspcet match finished")
    fo_log_key.flush()

    fo_log_key.write("\n\n Main F4 (suspect list) len of dic_cluster_total_input_idx_MOD_vs_node_info: " + str(
        len(dic_cluster_total_input_idx_MOD_vs_node_info)) + "\n")
    logger.info(f'Main F4 (suspect list) length of dic_cluster_total_input_idx_MOD_vs_node_info: '
                f'{len(dic_cluster_total_input_idx_MOD_vs_node_info)}')

    # TODO: T3DB matching (K. Hirata, 20221214)





    # [H1]-----------------------------------------------
    ## thresholding edges -----------------------
    #   list_of_edge_for_networkx_to_show  is list of list   [node_X, node_Y, dictionary {spec_sim_score, delta_mz, edge_type}]
    #  the edge type can be "inner_ref_layer", "inter_sample_ref_layer", "inner_sample_layer"
    logger.debug("[H1] Thresholding edges")

    list_of_edge_for_networkx_to_show, list_node_total_input_idx_mod_in_use = threshold_edges(dic_config,
                                                                                              list_of_edge_for_networkx)

    # what you get is [ node_id_X, node_id_Y, dictionary_for_attribute]

    fo_log.write("\n\n [H1]Show list_of_edge_for_networkx_to_show----------")
    log_message = "\n".join(map(str, list_of_edge_for_networkx_to_show))
    logger.info(f'[H1] Show list_of_edge_for_networkx_to_show\n{log_message}')
    for e in list_of_edge_for_networkx_to_show:
        fo_log.write("\n" + str(e))

    fo_log.write("\n  after threshold edge, len of node in use :" + str(len(list_node_total_input_idx_mod_in_use)))
    fo_log.flush()
    fo_log.write("\n  list_of_edge_for_networkx_to_show :" + str(list_of_edge_for_networkx_to_show))
    fo_log.flush()
    logger.info(f'After threshold edge, length of node in use: {len(list_node_total_input_idx_mod_in_use)}')

    #######################
    #  handling layer info

    # make list of layer attribute.
    list_attribute_for_layer = []
    logger.warning(f'len(dic_cluster_total_input_idx_MOD_vs_node_info): {len(dic_cluster_total_input_idx_MOD_vs_node_info)}')

    for total_input_idx_mod, node_o in dic_cluster_total_input_idx_MOD_vs_node_info.items():
        list_attribute_for_layer.append(node_o["layer"])
        list_attribute_for_layer = list(set(list_attribute_for_layer))

    fo_log.write("\n list_attribute_for_layer :" + str(list_attribute_for_layer))
    fo_log.flush()
    logger.warning(f'list_attribute_for_layer: {list_attribute_for_layer}')

    ###########################################
    # [J1]
    # make "Graph by layer".   also isolate inter layer network (edges + layers))
    ########################################
    logger.debug("[J1] making graph by layer")
    # attribute_for_layer, list_of_edge_for_networkx, list_node_info_for_layer

    ##  making dataset for inner layer edges, node.
    ##
    list_dic_edges_nodes_graph_by_layer = []
    # for each layer. (inner edges, and realate nodes...)
    for attribute_for_layer in list_attribute_for_layer:
        dic_edges_nodes_by_layer = {}
        dic_edges_nodes_by_layer["attribute_for_layer"] = attribute_for_layer
        fo_log.write("\n-------- [attribute_for_layer] : " + attribute_for_layer)

        # [EDGES]
        # create list_od_edge for newtorkx for particular layer. -----------------------------------
        list_of_edge_for_networkx_to_show_by_layer = []

        for edge_for_networkx in list_of_edge_for_networkx_to_show:

            # we are only looking for inner layer edge
            edge_for_networkx[2]["edge_type"] != "inter_sample_ref_layer"
            # if both node is on thi layer,

            if attribute_for_layer.startswith("sample"):
                fo_log.write(" \n attribute match 1:" + "---" +
                             dic_cluster_total_input_idx_MOD_vs_node_info[edge_for_networkx[0]][
                                 "layer"] + "---" + " VS " + "---" + attribute_for_layer + "---")
                fo_log.write(" \n attribute match 2:" + "---" +
                             dic_cluster_total_input_idx_MOD_vs_node_info[edge_for_networkx[1]][
                                 "layer"] + "---" + " VS " + "---" + attribute_for_layer + "---")

            if dic_cluster_total_input_idx_MOD_vs_node_info[edge_for_networkx[0]]["layer"] == attribute_for_layer \
                    and dic_cluster_total_input_idx_MOD_vs_node_info[edge_for_networkx[1]][
                "layer"] == attribute_for_layer:

                fo_log.write(" \n MATCHEDDDD")
                list_of_edge_for_networkx_to_show_by_layer.append(edge_for_networkx)
                if attribute_for_layer.startswith("sample"):
                    logger.debug("sample edge appendes..........................")
        dic_edges_nodes_by_layer["list_of_edge_for_networkx"] = list_of_edge_for_networkx_to_show_by_layer

        fo_log.write("\n  for this layer:   dic_edges_nodes_by_layer[list_of_edge_for_networkx] " + str(
            len(list_of_edge_for_networkx_to_show_by_layer)))
        logger.info(f"For '{attribute_for_layer}' layer, length of list_of_edge_for_networkx: "
                    f"{len(list_of_edge_for_networkx_to_show_by_layer)}")

        # [NODES]
        #  create dic  _cluster_total_input_idx_MOD_vs_node_info for each layer ------------------------
        dic_cluster_total_input_idx_MOD_vs_node_info_by_layer = {}
        for cluster_total_input_idx_MOD, node_info in dic_cluster_total_input_idx_MOD_vs_node_info.items():
            if node_info["layer"] == attribute_for_layer:
                dic_cluster_total_input_idx_MOD_vs_node_info_by_layer[cluster_total_input_idx_MOD] = node_info

        dic_edges_nodes_by_layer[
            "dic_cluster_total_input_idx_MOD_vs_node_info"] = dic_cluster_total_input_idx_MOD_vs_node_info_by_layer
        list_dic_edges_nodes_graph_by_layer.append(dic_edges_nodes_by_layer)

    fo_log.write(
        "\n\n\n  list_dic_edges_nodes_graph_by_layer " + str(len(list_dic_edges_nodes_graph_by_layer)) + "\n\n\n")
    logger.info(f'length of list_dic_edges_nodes_graph_by_layer: {len(list_dic_edges_nodes_graph_by_layer)}')
    logger.debug(f"len list_dic_edges_nodes_graph_by_layer: {str(len(list_dic_edges_nodes_graph_by_layer))}")
    fo_log.write("\n [Graph by layer]       len list_dic_edges_nodes_graph_by_layer" + str(
        len(list_of_edge_for_networkx_to_show_by_layer)))
    logger.info(f'[Graph by layer] length of list_dic_edges_nodes_graph_by_layer: '
                f'{len(list_of_edge_for_networkx_to_show_by_layer)}')

    fo_log_key.write("\n\n*end of main [J1] (making list_dic_edges_nodes_graph_by_layer )")
    log_message = f'End of main [J1] (making list_dic_edges_nodes_graph_by_layer )' \
                  f'\nattribute_for_layer\tlength of dic_cluster_total_input_idx_MOD_vs_node_info' \
                  f'\tlength of list_of_edge_for_networkx'
    for dic in list_dic_edges_nodes_graph_by_layer:
        fo_log_key.write(
            "\n" + dic["attribute_for_layer"] + " : " + "dic_cluster_total_input_idx_MOD_vs_node_info:" + str(
                len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])) \
            + ",    list_of_edge_for_networkx: " + str(len(dic["list_of_edge_for_networkx"])))

        log_message += f'\n{dic["attribute_for_layer"]}\t{len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])}' \
                       f'\t{len(dic["list_of_edge_for_networkx"])}'
    logger.info(log_message)

    #################################################################
    # [L1]
    ## making dataset for "INTER  SAMPLE VS REF  LAYER"  and "INTER SAMPLE" dataset.--------------------------------
    fo_log_key.write("\n\n [L1 ]making dataset for INTER  SAMPLE VS REF  LAYER  and INTER SAMPLE dataset")
    logger.debug("[L1] making dataset for inter sample vs ref layer  and inter sample dataset")
    logger.info("[L1] making dataset for inter sample vs ref layer  and inter sample dataset")
    list_dic_edges_nodes_inter_layer = []

    # [EDGES]
    # create list_od_edge for newtorkx for particular layer. -----------------------------------
    list_of_edge_for_networkx_to_show_inter_sample_ref_layer = []
    list_of_edge_for_networkx_to_show_inter_sample_layer = []
    # this list keeps total_input_idx_MOD for "inter sample ref layer"
    list_total_input_idx_MOD_inter_sample_ref_layer = []
    list_total_input_idx_MOD_inter_sample_layer = []

    logger.debug(f" len of list_of_edge_for_networkx_to_show {len(list_of_edge_for_networkx_to_show)}")
    for edge_for_networkx in list_of_edge_for_networkx_to_show:
        # we are only looking for inner layer edge

        if edge_for_networkx[2]["edge_type"] == "inter_sample_ref_layer":
            # if the edge is ""inter sample ref edge",
            list_of_edge_for_networkx_to_show_inter_sample_ref_layer.append(edge_for_networkx)
            list_total_input_idx_MOD_inter_sample_ref_layer.append(edge_for_networkx[0])
            list_total_input_idx_MOD_inter_sample_ref_layer.append(edge_for_networkx[1])

        if edge_for_networkx[2]["edge_type"] == "inter_sample_layer":
            list_of_edge_for_networkx_to_show_inter_sample_layer.append(edge_for_networkx)
            list_total_input_idx_MOD_inter_sample_layer.append(edge_for_networkx[0])
            list_total_input_idx_MOD_inter_sample_layer.append(edge_for_networkx[1])

    fo_log_key.write("\n list_of_edge_for_networkx_to_show_inter_sample_ref_layer : " + str(
        len(list_of_edge_for_networkx_to_show_inter_sample_ref_layer)))
    fo_log_key.write("\n list_of_edge_for_networkx_to_show_inter_sample_layer : " + str(
        len(list_of_edge_for_networkx_to_show_inter_sample_layer)))
    logger.info(f'Length of list_of_edge_for_networkx_to_show_inter_sample_ref_layer: '
                f'{len(list_of_edge_for_networkx_to_show_inter_sample_ref_layer)}')
    logger.info(f'Length of ist_of_edge_for_networkx_to_show_inter_sample_layer: '
                f'{len(list_of_edge_for_networkx_to_show_inter_sample_layer)}')

    # make it non redundant.
    list_total_input_idx_MOD_inter_sample_ref_layer = list(set(list_total_input_idx_MOD_inter_sample_ref_layer))
    list_total_input_idx_MOD_inter_sample_layer = list(set(list_total_input_idx_MOD_inter_sample_layer))

    # [NODES]
    #  create dic  _cluster_total_input_idx_MOD_vs_node_info for each layer ------------------------
    dic_cluster_total_input_idx_MOD_vs_node_info_inter_sample_ref_layer = {}
    for cluster_total_input_idx_MOD, node_info in dic_cluster_total_input_idx_MOD_vs_node_info.items():
        if cluster_total_input_idx_MOD in list_total_input_idx_MOD_inter_sample_ref_layer:
            dic_cluster_total_input_idx_MOD_vs_node_info_inter_sample_ref_layer[cluster_total_input_idx_MOD] = node_info

    dic_cluster_total_input_idx_MOD_vs_node_info_inter_sample_layer = {}
    for cluster_total_input_idx_MOD, node_info in dic_cluster_total_input_idx_MOD_vs_node_info.items():
        if cluster_total_input_idx_MOD in list_total_input_idx_MOD_inter_sample_layer:
            dic_cluster_total_input_idx_MOD_vs_node_info_inter_sample_layer[cluster_total_input_idx_MOD] = node_info

    #############################################################
    ########################################################
    # M1
    # [M1] crate networkx graph for each layer
    logger.debug("[multilayer_3d_network_b1/process_3d_network_data]  [M1] crate networkx graph for each layer ")
    for dic_edges_nodes_graph_by_layer in list_dic_edges_nodes_graph_by_layer:
        fo_log.write("\n creating FG, whoing edges:" + str(dic_edges_nodes_graph_by_layer["attribute_for_layer"]) + str(
            dic_edges_nodes_graph_by_layer["list_of_edge_for_networkx"]))
        logger.info(f"Creating FG for '{dic_edges_nodes_graph_by_layer['attribute_for_layer']}, "
                    f"{dic_edges_nodes_graph_by_layer['list_of_edge_for_networkx']}'")
        FG = create_networkx_graph(dic_config,
                                   dic_edges_nodes_graph_by_layer["dic_cluster_total_input_idx_MOD_vs_node_info"], \
                                   dic_edges_nodes_graph_by_layer["list_of_edge_for_networkx"])
        dic_edges_nodes_graph_by_layer["nx_graph"] = FG

    fo_log.write("\nnumber of network (layers )  len of (dic_edges_nodes_graph_by_layer):" + str(
        len(dic_edges_nodes_graph_by_layer)) + "\n")
    fo_log.flush()
    logger.info(f'Number of network layers (length of dic_edges_nodes_graph_by_layer): '
                f'{len(dic_edges_nodes_graph_by_layer)}')

    ################################
    # [O] community detection
    logger.debug("[O]staring community detection ")

    fo_log_key.write("\n\n beginning of community detection [O]")
    log_message = f'Beginning of community detection [O]' \
                  f'\nattribute_for_layer\tlength of dic_cluster_total_input_idx_MOD_vs_node_info' \
                  f'\tlength of list_of_edge_for_networkx'
    for dic in list_dic_edges_nodes_graph_by_layer:
        fo_log_key.write(
            "\n" + dic["attribute_for_layer"] + " : " + "dic_cluster_total_input_idx_MOD_vs_node_info:" + str(
                len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])) \
            + ",    list_of_edge_for_networkx: " + str(len(dic["list_of_edge_for_networkx"])))
        log_message += f'\n{dic["attribute_for_layer"]}\t{len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])}' \
                       f'\t{len(dic["list_of_edge_for_networkx"])}'
    logger.info(log_message)

    from networkx.readwrite import json_graph
    from networkx.algorithms import community
    import json

    n_level_community_detection = dic_config["n_level_community_detection"]

    fo_log_key.write("\n level of community detection (0: no community detection): " + str(n_level_community_detection))
    logger.info(f'Level of community detection (0: no community detection): {n_level_community_detection}')

    # n_level_community_detection = 1

    if n_level_community_detection > 0:
        fo_log_key.write("\n performing community detection at level : " + str(n_level_community_detection))
        logger.info(f'Performing community detection at level: {n_level_community_detection}')

        str_log = "\n\nCommunity detection===============" + "\n"
        #####
        # [O1] first,  combine all nodes, edges ---------------------------------
        logger.debug("  [O1] combine all nodes, edges ")
        G_all_combined = nx.Graph()

        # combine node/edges for each layers. (not inter layer edges)
        for dic_edges_nodes_graph_by_layer in list_dic_edges_nodes_graph_by_layer:
            G_all_combined = nx.compose(G_all_combined, dic_edges_nodes_graph_by_layer["nx_graph"])

        # add inter sample-ref edges
        G_all_combined.add_edges_from(list_of_edge_for_networkx_to_show_inter_sample_ref_layer)

        # add inter sample-ref edges
        G_all_combined.add_edges_from(list_of_edge_for_networkx_to_show_inter_sample_layer)

        # [O3] Perform community detection ----------------------------------------------
        logger.debug("  [O3] Perform community detection")

        l_communities = []
        community_detection_method = "greedy_modularity_communities"

        if community_detection_method == "girvan_newman":
            communities_generator = community.girvan_newman(G_all_combined)
            # return as tuple
            for n in range(n_level_community_detection):
                l_communities = next(communities_generator)

        if community_detection_method == "greedy_modularity_communities":
            from networkx.algorithms.community import greedy_modularity_communities
            try:
                l_communities = list(greedy_modularity_communities(G_all_combined))
            except:
                pass

        logger.debug("  [O3] community detection finished")

        logger.debug("  [O3] update edges and graph based on the detected communiteis")
        if len(l_communities) > 0:

            str_log = str_log + "\n detected communities \n"
            for c in l_communities:
                str_log = str_log + "\n" + str(c)

            # [O5]  check all the edges whether it fits community definition -----------------

            # this holds edges that do not fit within community.  in other word, inter community edges
            l_inter_community_edge_as_tuple = []
            logger.debug("         iterate edges in G_all_combined")
            # iterate edges----------------------
            for u, v, d in G_all_combined.edges(data=True):
                f_in_community = 0
                # iterate communities
                for commu in l_communities:
                    # !!!!!!!!!! 20220903 changed
                    if u in commu and v in commu:
                        # if all(x in commu for x in [u, v]):
                        f_in_community = 1
                # if the current edge is not inside community,
                if f_in_community == 0:
                    l_inter_community_edge_as_tuple.append((u, v, d))

            logger.debug(f"inter-community edges {str(l_inter_community_edge_as_tuple)}")
            str_log = str_log + "\n inter-community edges " + str(len(l_inter_community_edge_as_tuple)) + "\n"

            fo_log_key.write(str_log)
            fo_log_key.flush()
            logger.debug(str_log)

            # [O7]: update edge info-------------------------------------------------------------------
            logger.debug("  [O7]  update edge info")
            # for each layer.  (list_dic_edges_nodes_graph_by_layer)
            for dic_edges_nodes_graph_by_layer in list_dic_edges_nodes_graph_by_layer:
                l_inter_community_edge = []
                # list of edge within community, which supposed to be preserved.
                l_inner_community_edge = []
                # iterate edges
                for edge in dic_edges_nodes_graph_by_layer["list_of_edge_for_networkx"]:
                    f_in_community = 0
                    # check the edge is in community.
                    for commu in l_communities:

                        # !!!!!!!!!! 20220903 changed
                        if edge[0] in commu and edge[1] in commu:
                            # if all(x in commu for x in [edge[0], edge[1]]):
                            l_inner_community_edge.append((edge[0], edge[1], edge[2]))
                            f_in_community = 1

                    if f_in_community == 0:
                        l_inter_community_edge.append((edge[0], edge[1], edge[2]))

                dic_edges_nodes_graph_by_layer["list_of_edge_for_networkx"] = l_inner_community_edge
                # it seems you need to unfrozen graph
                nx_graph_unfrozen = nx.Graph(dic_edges_nodes_graph_by_layer["nx_graph"])
                # then remove edge.
                nx_graph_unfrozen.remove_edges_from(l_inter_community_edge)
                dic_edges_nodes_graph_by_layer["nx_graph"] = nx_graph_unfrozen

                # log
                # str_log =  str_log  + "layer attribute: " + dic_edges_nodes_graph_by_layer["attribute_for_layer"] + "\n"
                # str_log = str_log + " num edges removed (in community detection mode: do not belong to community): " +  str(len(l_inter_community_edge))+ "\n"

            #    remove inter community edges for list_of_edge_for_networkx_to_show_inter_sample_ref_layer
            l_inner_community_edge_inter_sample_ref_layer = []
            l_inter_community_edge_inter_sample_ref_layer = []
            for edge in list_of_edge_for_networkx_to_show_inter_sample_ref_layer:
                f_in_community = 0
                # check the edge is in community.
                for commu in l_communities:
                    if all(x in commu for x in [edge[0], edge[1]]):
                        l_inner_community_edge_inter_sample_ref_layer.append((edge[0], edge[1], edge[2]))
                        f_in_community = 1

                if f_in_community == 0:
                    l_inter_community_edge_inter_sample_ref_layer.append((edge[0], edge[1], edge[2]))
            # take over
            list_of_edge_for_networkx_to_show_inter_sample_ref_layer = l_inner_community_edge_inter_sample_ref_layer

            # logger.debug("  [O7]  remove inter community edges  for list_of_edge_for_networkx_to_show_inter_sample_layer ")
            #   remove inter community edges  for list_of_edge_for_networkx_to_show_inter_sample_layer
            l_inner_community_edge_inter_sample_layer = []
            l_inter_community_edge_inter_sample_layer = []
            for edge in list_of_edge_for_networkx_to_show_inter_sample_layer:
                f_in_community = 0
                # check the edge is in community.
                for commu in l_communities:
                    if all(x in commu for x in [edge[0], edge[1]]):
                        l_inner_community_edge_inter_sample_layer.append((edge[0], edge[1], edge[2]))
                        f_in_community = 1

                if f_in_community == 0:
                    l_inter_community_edge_inter_sample_layer.append((edge[0], edge[1], edge[2]))
            # take over
            list_of_edge_for_networkx_to_show_inter_sample_layer = l_inner_community_edge_inter_sample_layer

        ## export
        logger.debug("    exporting to json")
        data1 = json_graph.node_link_data(G_all_combined)
        fp = open("(G_all_combined.json", "w")
        json.dump(data1, fp)
        logger.debug("    finished exporting to json")
        logger.debug(" Finishing community detection")

    logger.debug("[O3] community detection finished")
    #  end community detection

    fo_log_key.write("\n\n\n after community detection")
    fo_log_key.write("\n list_of_edge_for_networkx_to_show_inter_sample_ref_layer : " + str(
        len(list_of_edge_for_networkx_to_show_inter_sample_ref_layer)))
    logger.info(f'After community detection, length of list_of_edge_for_networkx_to_show_inter_sample_ref_layer: '
                f'{len(list_of_edge_for_networkx_to_show_inter_sample_ref_layer)}')

    ########
    # [N1]  SUBGRAPH
    ############################################################
    # SUBGRAPH  (sample)
    #  here you are making list of global_accession of which mass feature show significance quantitative shift.

    fo_log_key.write("\n\n [N1]  SUBGRAPH")
    log_message = f'[N1] SUBGRAPH\nattribute_for_layer\tlength of dic_cluster_total_input_idx_MOD_vs_node_info' \
                  f'\tlist_of_edge_for_networkx'
    for dic in list_dic_edges_nodes_graph_by_layer:
        fo_log_key.write(
            "\n" + dic["attribute_for_layer"] + " : " + "dic_cluster_total_input_idx_MOD_vs_node_info:" + str(
                len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])) \
            + ",    list_of_edge_for_networkx: " + str(len(dic["list_of_edge_for_networkx"])))
        log_message += f'\n{dic["attribute_for_layer"]}\t{len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])}' \
                       f'\t{len(dic["list_of_edge_for_networkx"])}'
    logger.info(log_message)
    logger.debug("[multilayer_3d_network_b1/process_3d_network_data]  [N1] creating [SUBGRAPH  (sample)] ")

    l_total_input_idx_mod_sample = []

    fo_log_key.write("\n\nconfig_o.stat_value: " + str(dic_config["stat_value"]) + "   " + str(
        type(dic_config["stat_value"])) + "\n")
    fo_log_key.write("config_o.quant_value: " + str(dic_config["quant_value"]) + "\n")
    fo_log_key.write("topn: " + str(dic_config["subgraph_num_core_nodes"]) + "\n")
    fo_log_key.write("config_o.subgraph_depth: " + str(dic_config["subgraph_depth"]) + "\n")

    logger.info(f'config_o.stat_value: {dic_config["stat_value"]}')
    logger.info(f'config_o.quant_value: {dic_config["quant_value"]}')
    logger.info(f'config_o.subgraph_num_core_nodes: {dic_config["subgraph_num_core_nodes"]}')
    logger.info(f'config_o.subgraph_depth: {dic_config["subgraph_depth"]}')
    logger.info(f'config_o.quant_polar: {dic_config["quant_polar"]}')

    ##########################
    # [N1a] Quantitative subgraph
    ############################

    # foe each layer
    for dic_edges_nodes_graph_by_layer in list_dic_edges_nodes_graph_by_layer:

        fo_log_key.write("\ndic_edges_nodes_graph_by_layer[attribute_for_layer]: " + dic_edges_nodes_graph_by_layer[
            "attribute_for_layer"])
        # if the current layer is "sample"------------------------------------------------------
        if dic_edges_nodes_graph_by_layer["attribute_for_layer"].startswith("sample"):
            l_global_accession_for_subgraph_sample_user_selected = []
            # making dictionary where key is global accession and value is "mf.value_to_show"
            dic_global_accession_vs_value_to_show = {}
            # this has p-value for mass feature
            dic_global_accession_vs_stat_value = {}

            for k, mf in dic_source_data["dic_global_accession_vs_mass_feature"].items():
                dic_global_accession_vs_value_to_show[k] = mf['value_to_show']
                dic_global_accession_vs_stat_value[k] = mf['stat_val']

            # this list hold global accession for base of sugraph(quant)
            l_global_accession_for_subgraph_sample_user_selected_quant = []

            if dic_config["subgraph_type"] in ["node_quant"]:
                top_n = dic_config["subgraph_num_core_nodes"]

                if dic_config["quant_polar"] == "quant_val_both":
                    top_n = top_n / 2

                ###########
                # select topN (bipolar)
                #######################################
                # small to large
                if dic_config["quant_polar"] == "quant_val_both" or dic_config["quant_polar"] == "quant_val_decrease":
                    count = 0
                    for k, v in sorted(list(dic_global_accession_vs_value_to_show.items()), key=lambda x: x[1]):
                        if count >= top_n:
                            break

                        stat_val = dic_global_accession_vs_stat_value[k]

                        if dic_config['quant_value_type'] == 'ratio':
                            threshold_quant_value = 1.0 / dic_config["quant_value"]
                        else:
                            threshold_quant_value = - abs(dic_config["quant_value"])

                        # logger.warning(f'@@@D v: {v}, threshold_quant_value: {threshold_quant_value}, count: {count}')
                        # if statistaial value is lower than threshold,  AND  the value is  lower thatn target RATIO. NOTE this is ratio.
                        if stat_val < dic_config["stat_value"] and v < threshold_quant_value:
                            l_global_accession_for_subgraph_sample_user_selected_quant.append(k)
                            count = count + 1

                if dic_config["quant_polar"] == "quant_val_both" or dic_config["quant_polar"] == "quant_val_increase":

                    fo_log_key.write("\n quant mode: both or val increase")
                    count = 0
                    # large to small
                    for k, v in sorted(list(dic_global_accession_vs_value_to_show.items()),
                                       reverse=True, key=lambda x: x[1]):
                        if count >= top_n:
                            break

                        stat_val = dic_global_accession_vs_stat_value[k]
                        if dic_config['quant_value_type'] == 'ratio':
                            threshold_quant_value = dic_config["quant_value"]
                        else:
                            threshold_quant_value = abs(dic_config["quant_value"])

                        # logger.warning(f'@@@I v: {v}, threshold_quant_value: {threshold_quant_value}, count: {count}')
                        # if statistaial value is lower than threshold,  AND  the value is  hier than thatn target RATIO. NOTE this is ratio.
                        if stat_val < dic_config["stat_value"] and v > threshold_quant_value:
                            l_global_accession_for_subgraph_sample_user_selected_quant.append(k)
                            count = count + 1

                l_global_accession_for_subgraph_sample_user_selected_quant = list(
                    set(l_global_accession_for_subgraph_sample_user_selected_quant))

                fo_log_key.write("\nlen(l_global_accession_for_subgraph_sample_user_selected_quant): " + str(
                    len(l_global_accession_for_subgraph_sample_user_selected_quant)))
                fo_log_key.flush()
                logger.warning(f"In '{dic_edges_nodes_graph_by_layer['attribute_for_layer']}' layer, "
                            f"length of l_global_accession_for_subgraph_sample_user_selected_quant: "
                            f"{len(l_global_accession_for_subgraph_sample_user_selected_quant)}")

                ############
                # get base for qunt-based sugraph
                ####################################################################
                # create subgraph---------------------------------------------
                if len(l_global_accession_for_subgraph_sample_user_selected_quant) > 0:
                    list_total_input_idx_mod_for_subgraph = []

                    nx_graph = dic_edges_nodes_graph_by_layer["nx_graph"]

                    # l_global_accession_for_subgraph_sample_user_selected  is the list of "base" (starting point of subgraph)
                    for node in nx_graph.nodes(data=True):
                        if node[1]['global_accession'] in l_global_accession_for_subgraph_sample_user_selected_quant:
                            list_total_input_idx_mod_for_subgraph.append(node[0])

                    fo_log_key.write("\n len of base ( list_total_input_idx_mod_for_subgraph) : " + str(
                        len(list_total_input_idx_mod_for_subgraph)))
                    fo_log_key.flush()
                    logger.info(f"In '{dic_edges_nodes_graph_by_layer['attribute_for_layer']}' layer, "
                                f"lenght of base  (list_total_input_idx_mod_for_subgraph): "
                                f"{len(list_total_input_idx_mod_for_subgraph)}")
                    base = list_total_input_idx_mod_for_subgraph

                    foundset = {key for source in base for key in
                                list(nx.single_source_shortest_path(nx_graph, source,
                                                                    cutoff=dic_config["subgraph_depth"]).keys())}

                    # update and replace nx graph and other info
                    nx_graph_sub = nx_graph.subgraph(foundset)
                    dic_edges_nodes_graph_by_layer["nx_graph"] = nx_graph_sub

                    # create list of edges again.  (nx.generate_edgelist  not really working ???)----------------
                    list_of_edge_for_networkx_sub = []
                    for e in nx_graph_sub.edges.data():
                        list_of_edge_for_networkx_sub.append([e[0], e[1], e[2]])

                    dic_edges_nodes_graph_by_layer["list_of_edge_for_networkx"] = list_of_edge_for_networkx_sub

                    # updating "dic_cluster_total_input_idx_MOD_vs_node_info"-----------------------------

                    # this has total input idx mod PRESENT in subgraph
                    l_total_input_idx_mod_IN_subgraph = []
                    for node in nx_graph_sub.nodes(data=True):
                        l_total_input_idx_mod_IN_subgraph.append(node[0])
                        l_total_input_idx_mod_sample.append(node[0])

                    dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE = {}
                    for cluster_total_input_idx_MOD, node_info in dic_edges_nodes_graph_by_layer[
                        "dic_cluster_total_input_idx_MOD_vs_node_info"].items():
                        if cluster_total_input_idx_MOD in l_total_input_idx_mod_IN_subgraph:
                            dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE[
                                cluster_total_input_idx_MOD] = node_info

                    dic_edges_nodes_graph_by_layer[
                        "dic_cluster_total_input_idx_MOD_vs_node_info"] = dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE

    ##################################
    ################################
    # [N1b] user select node-based subgraph
    ###################################

    fo_log_key.write("\n\n [N1b] user select node-based subgraph")
    log_message = f'[N1b] user select node-based subgraph\nattribute_for_layer' \
                  f'\tlength of dic_cluster_total_input_idx_MOD_vs_node_info\tlength of list_of_edge_for_networkx'
    for dic in list_dic_edges_nodes_graph_by_layer:
        fo_log_key.write(
            "\n" + dic["attribute_for_layer"] + " : " + "dic_cluster_total_input_idx_MOD_vs_node_info:" + str(
                len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])) \
            + ",    list_of_edge_for_networkx: " + str(len(dic["list_of_edge_for_networkx"])))
        log_message += f'\n{dic["attribute_for_layer"]}\t{len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])}' \
                       f'\t{len(dic["list_of_edge_for_networkx"])}'
    logger.info(log_message)

    logger.debug("[multilayer_3d_network_b1/process_3d_network_data]  [N1b] user select node-based subgraph")

    ## combine all sample data
    # first, assemble all nodes, edges from sample tag containing dataset
    dic_cluster_total_input_idx_MOD_vs_node_info_all_samples = {}
    list_of_edge_for_networkx_all_samples = []

    for dic in list_dic_edges_nodes_graph_by_layer:
        # dealing with sample layer data.
        # note if you have more than 1, you have to merge and make layout for the merged datast.
        if dic["attribute_for_layer"].startswith("sample"):
            fo_log.write("\n\nliterating : " + str(dic["attribute_for_layer"]))
            fo_log.write("\ndic[dic_cluster_total_input_idx_MOD_vs_node_info]: " + str(
                dic["dic_cluster_total_input_idx_MOD_vs_node_info"]))
            fo_log.write("\ndic[list_of_edge_for_networkx]: " + str(dic["list_of_edge_for_networkx"]))
            logger.info(f'Iterating: {dic["attribute_for_layer"]}\n'
                        f'\ndic["dic_cluster_total_input_idx_MOD_vs_node_info"]: '
                        f'{dic["dic_cluster_total_input_idx_MOD_vs_node_info"]}'
                        f'\ndic[list_of_edge_for_networkx]: {dic["list_of_edge_for_networkx"]}')

            dic_cluster_total_input_idx_MOD_vs_node_info_all_samples.update(
                dic["dic_cluster_total_input_idx_MOD_vs_node_info"])
            list_of_edge_for_networkx_all_samples = list_of_edge_for_networkx_all_samples + dic[
                "list_of_edge_for_networkx"]
    list_of_edge_for_networkx_all_samples = list_of_edge_for_networkx_all_samples + list_of_edge_for_networkx_to_show_inter_sample_layer
    fo_log.write("\n\nlist_of_edge_for_networkx_all_samples: " + str(list_of_edge_for_networkx_all_samples))
    logger.info(f'list_of_edge_for_networkx_all_samples: {list_of_edge_for_networkx_all_samples}')

    # just converting to networkx graph.  doing nothing with subgraph etc.
    FG_all_samples = create_networkx_graph(dic_config, dic_cluster_total_input_idx_MOD_vs_node_info_all_samples, \
                                           list_of_edge_for_networkx_all_samples)

    ### create subgraph  for sample-combined dataset-----------------------------------------------------
    l_global_accession_for_subgraph_sample_user_selected = dic_config["l_global_accession_for_node_select_subgraph"]
    l_total_idx_mod_user_select_subgraph_all_sample = []

    if len(l_global_accession_for_subgraph_sample_user_selected) > 0:
        list_total_input_idx_mod_for_subgraph = []
        nx_graph = FG_all_samples
        fo_log.write("nx_graph.data(): " + str(nx_graph.edges.data()))
        logger.info(f'nx_graph.data(): {nx_graph.edges.data()}')

        # l_global_accession_for_subgraph_sample_user_selected  is the list of "base" (starting point of subgraph)
        # node[0] is total_input_idx_mod
        for node in nx_graph.nodes(data=True):
            fo_log.write("\n  FG_allsamples node:" + str(node))
            logger.info(f'FG_allsamples node: {str(node)}')
            if node[1]['global_accession'] in l_global_accession_for_subgraph_sample_user_selected:
                list_total_input_idx_mod_for_subgraph.append(node[0])
        base = list_total_input_idx_mod_for_subgraph

        fo_log.write("CREATING SUER DEFINED SUBGRAPH:" + "suggraphdepth:" + str(
            dic_config["node_select_subgraph_depth"]) + "   base:" + str(base) + "\n")
        fo_log.flush()
        logger.info(f'CREATING SUER DEFINED SUBGRAPH'
                    f'\nconfig_o.node_select_subgraph_depth: {dic_config["node_select_subgraph_depth"]}'
                    f'\nbase: {base}')

        foundset = {key for source in base for key in
                    list(nx.single_source_shortest_path(nx_graph, source,
                                                        cutoff=dic_config["node_select_subgraph_depth"]).keys())}
        # update and replace nx graph and other info

        fo_log.write("\n\nfound set:" + str(foundset))
        logger.info(f'Found set: {foundset}')
        nx_graph_sub = nx_graph.subgraph(foundset)
        FG_all_samples = nx_graph_sub

        # l_total_idx_mod_user_select_subgraph_all_sample   holds all idx mod to be extracted

        for n in nx_graph_sub.nodes.data():
            l_total_idx_mod_user_select_subgraph_all_sample.append(n[0])

        # foe each layer
        for dic_edges_nodes_graph_by_layer in list_dic_edges_nodes_graph_by_layer:
            # if the current layer is "sample"------------------------------------------------------
            if dic_edges_nodes_graph_by_layer["attribute_for_layer"].startswith("sample"):

                ####################################################################
                # create subgraph---------------------------------------------
                list_total_input_idx_mod_for_subgraph = []
                nx_graph = dic_edges_nodes_graph_by_layer["nx_graph"]
                nx_graph_sub = nx_graph.subgraph(l_total_idx_mod_user_select_subgraph_all_sample)
                dic_edges_nodes_graph_by_layer["nx_graph"] = nx_graph_sub

                # create list of edges again.  nx.generate_edgelist  not really working ???
                list_of_edge_for_networkx_sub = []
                for e in nx_graph_sub.edges.data():
                    list_of_edge_for_networkx_sub.append([e[0], e[1], e[2]])

                dic_edges_nodes_graph_by_layer["list_of_edge_for_networkx"] = list_of_edge_for_networkx_sub

                # updating "dic_cluster_total_input_idx_MOD_vs_node_info"-----------------------------

                # this has total input idx mod PRESENT in subgraph
                l_total_input_idx_mod_IN_subgraph = []

                for node in nx_graph_sub.nodes(data=True):
                    l_total_input_idx_mod_IN_subgraph.append(node[0])
                    l_total_input_idx_mod_sample.append(node[0])
                dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE = {}
                for cluster_total_input_idx_MOD, node_info in dic_edges_nodes_graph_by_layer[
                    "dic_cluster_total_input_idx_MOD_vs_node_info"].items():
                    if cluster_total_input_idx_MOD in l_total_input_idx_mod_IN_subgraph:
                        dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE[cluster_total_input_idx_MOD] = node_info

                dic_edges_nodes_graph_by_layer[
                    "dic_cluster_total_input_idx_MOD_vs_node_info"] = dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE

                ##########################################

    fo_log.write("\n\n len of list_of_edge_for_networkx_to_show_inter_sample_ref_layer:" + str(
        len(list_of_edge_for_networkx_to_show_inter_sample_ref_layer)))
    fo_log.write("\n\n list_of_edge_for_networkx_to_show_inter_sample_ref_layer" + str(
        list_of_edge_for_networkx_to_show_inter_sample_ref_layer))
    logger.info(f'Length of list_of_edge_for_networkx_to_show_inter_sample_ref_layer: '
                f'{len(list_of_edge_for_networkx_to_show_inter_sample_ref_layer)}')
    logger.info(f'list_of_edge_for_networkx_to_show_inter_sample_ref_layer: '
                f'{list_of_edge_for_networkx_to_show_inter_sample_ref_layer}')

    fo_log_key.write("\n\n len of list_of_edge_for_networkx_to_show_inter_sample_ref_layer: " + str(len(
        list_of_edge_for_networkx_to_show_inter_sample_ref_layer)))
    fo_log_key.flush()

    ##########################################
    # [N1r]create subgraph for ref layer
    ##########################################

    fo_log_key.write("\n\n Just before [N1r]create subgraph for ref layer")
    log_message = f'Just before [N1r]create subgraph for ref layer' \
                  f'\nattribute_for_layer\tlength of dic_cluster_total_input_idx_MOD_vs_node_info' \
                  f'\tlength of dic_cluster_total_input_idx_MOD_vs_node_info'
    for dic in list_dic_edges_nodes_graph_by_layer:
        fo_log_key.write(
            "\n" + dic["attribute_for_layer"] + " : " + "dic_cluster_total_input_idx_MOD_vs_node_info:" + str(
                len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])) \
            + ",    list_of_edge_for_networkx: " + str(len(dic["list_of_edge_for_networkx"])))
        log_message += f'\n{dic["attribute_for_layer"]}\t{len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])}' \
                       f'\t{len(dic["list_of_edge_for_networkx"])}'
    logger.info(log_message)

    logger.debug("[multilayer_3d_network_b1/process_3d_network_data]  [N1r create subgraph for ref layer")
    # !!!!!!!!!!!!!!!! NOTE !!!!!!!!!!!!!
    #  l_total_input_idx_mod_sample is alrady only contain use-selected nodes, is user use that function. !
    logger.debug("starting [create subgraph for ref layer]")
    fo_log_key.write("\n\n[create subgraph for ref layer]\n")
    logger.info('[create subgraph for ref layer]')
    # get inter layer edges.
    l_total_input_idx_mod_ref_for_subgraph_base = []

    # iterate inter sample-ref edges.
    for inter_ref_sample_edge in list_of_edge_for_networkx_to_show_inter_sample_ref_layer:
        # if node0 is present in current l_total_input_idx_mod_sample
        if inter_ref_sample_edge[0] in l_total_input_idx_mod_sample:
            # then node1 should be on ref layer, since this is sample ref inter edge.
            l_total_input_idx_mod_ref_for_subgraph_base.append(inter_ref_sample_edge[1])
        if inter_ref_sample_edge[1] in l_total_input_idx_mod_sample:
            l_total_input_idx_mod_ref_for_subgraph_base.append(inter_ref_sample_edge[0])

    fo_log.write(
        "\n\n" + "l_total_input_idx_mod_ref_for_subgraph_base\n" + str(l_total_input_idx_mod_ref_for_subgraph_base))
    fo_log_key.write("\n" + " len of list_of_edge_for_networkx_to_show_inter_sample_ref_layer:" + str(
        len(list_of_edge_for_networkx_to_show_inter_sample_ref_layer)) + "\n")
    fo_log_key.write("\n" + " len of l_total_input_idx_mod_ref_for_subgraph_base:" + str(
        len(l_total_input_idx_mod_ref_for_subgraph_base)) + "\n")
    fo_log_key.write("\n" + " len of l_total_input_idx_mod_sample:" + str(len(l_total_input_idx_mod_sample)) + "\n")

    fo_log_key.flush()

    logger.info(f'l_total_input_idx_mod_ref_for_subgraph_base: {l_total_input_idx_mod_ref_for_subgraph_base}')
    logger.info(f'Length of list_of_edge_for_networkx_to_show_inter_sample_ref_layer: '
                f'{len(list_of_edge_for_networkx_to_show_inter_sample_ref_layer)}')
    logger.info(f'Length of l_total_input_idx_mod_ref_for_subgraph_base: '
                f'{len(l_total_input_idx_mod_ref_for_subgraph_base)}')
    logger.info(f'Length of l_total_input_idx_mod_sample: {len(l_total_input_idx_mod_sample)}')

    ########################3
    # preserve nodes, edges only connected to sample layer
    depth_ref_preserve_interlayer = 10
    f_preserve_only_ref_node_inter_layer_connected = 1

    if f_preserve_only_ref_node_inter_layer_connected > 0:
        # look for ref layer---------------------------------
        for dic_edges_nodes_graph_by_layer in list_dic_edges_nodes_graph_by_layer:

            # choosing ref layer (layer tag is not sample.)
            if not dic_edges_nodes_graph_by_layer["attribute_for_layer"].startswith("sample"):

                # selecting node
                list_total_input_idx_mod_for_base = []
                for cluster_total_input_idx_MOD, node_info in dic_edges_nodes_graph_by_layer[
                    "dic_cluster_total_input_idx_MOD_vs_node_info"].items():

                    # if the ref node is in inter ref-sample node list
                    if cluster_total_input_idx_MOD in list_total_input_idx_MOD_inter_sample_ref_layer:
                        list_total_input_idx_mod_for_base.append(cluster_total_input_idx_MOD)

                # creating subgraph.
                nx_graph = dic_edges_nodes_graph_by_layer["nx_graph"]
                base = list_total_input_idx_mod_for_base
                foundset = {key for source in base for key in
                            list(nx.single_source_shortest_path(nx_graph, source,
                                                                cutoff=depth_ref_preserve_interlayer).keys())}
                # update
                nx_graph_sub = nx_graph.subgraph(foundset)
                dic_edges_nodes_graph_by_layer["nx_graph"] = nx_graph_sub

                # create list of edges again.  nx.generate_edgelist  not really working ???
                list_of_edge_for_networkx_sub = []
                for e in nx_graph_sub.edges.data():
                    list_of_edge_for_networkx_sub.append([e[0], e[1], e[2]])

                # update_edges
                dic_edges_nodes_graph_by_layer["list_of_edge_for_networkx"] = list_of_edge_for_networkx_sub

                # this has total input idx mod PRESENT in subgraph
                l_total_input_idx_mod_IN_subgraph = []

                for node in nx_graph_sub.nodes(data=True):
                    l_total_input_idx_mod_IN_subgraph.append(node[0])
                dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE = {}

                for cluster_total_input_idx_MOD, node_info in dic_edges_nodes_graph_by_layer[
                    "dic_cluster_total_input_idx_MOD_vs_node_info"].items():
                    if cluster_total_input_idx_MOD in l_total_input_idx_mod_IN_subgraph:
                        dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE[cluster_total_input_idx_MOD] = node_info

                dic_edges_nodes_graph_by_layer[
                    "dic_cluster_total_input_idx_MOD_vs_node_info"] = dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE

    #########################################################
    # This is for extracting compounds and related, BASED ON EXTERNAL CMPD INFO

    # look for ref layer---------------------------------
    for dic_edges_nodes_graph_by_layer in list_dic_edges_nodes_graph_by_layer:
        # choosing ref layer (layer tag is not sample.)
        if not dic_edges_nodes_graph_by_layer["attribute_for_layer"].startswith("sample"):
            # if "external compound info" is selected to make subgraph
            if dic_config.get("ref_filter_ext") == "ext_cmpd_info":
                # !!!!!!!!!!!!!!!! tempo
                ref_subgraph_depth = 3
                # selecting node based on external compound info.
                list_total_input_idx_mod_for_subgraph = []
                for cluster_total_input_idx_MOD, node_info in dic_edges_nodes_graph_by_layer[
                    "dic_cluster_total_input_idx_MOD_vs_node_info"].items():
                    if len(node_info.spec_cluster.list_compound_categories) > 0:
                        list_total_input_idx_mod_for_subgraph.append(cluster_total_input_idx_MOD)
                nx_graph = dic_edges_nodes_graph_by_layer["nx_graph"]

                base = list_total_input_idx_mod_for_subgraph

                foundset = {key for source in base for key in
                            list(nx.single_source_shortest_path(nx_graph, source, cutoff=ref_subgraph_depth).keys())}

                # update and replace nx graph and other info
                nx_graph_sub = nx_graph.subgraph(foundset)
                dic_edges_nodes_graph_by_layer["nx_graph"] = nx_graph_sub

                # create list of edges again.  nx.generate_edgelist  not really working ???
                list_of_edge_for_networkx_sub = []
                for e in nx_graph_sub.edges.data():
                    list_of_edge_for_networkx_sub.append([e[0], e[1], e[2]])

                # update_edges
                dic_edges_nodes_graph_by_layer["list_of_edge_for_networkx"] = list_of_edge_for_networkx_sub

                # updating "dic_cluster_total_input_idx_MOD_vs_node_info"-----------------------------

                # this has total input idx mod PRESENT in subgraph
                l_total_input_idx_mod_IN_subgraph = []

                for node in nx_graph_sub.nodes(data=True):
                    l_total_input_idx_mod_IN_subgraph.append(node[0])
                dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE = {}
                for cluster_total_input_idx_MOD, node_info in dic_edges_nodes_graph_by_layer[
                    "dic_cluster_total_input_idx_MOD_vs_node_info"].items():

                    if cluster_total_input_idx_MOD in l_total_input_idx_mod_IN_subgraph:
                        dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE[cluster_total_input_idx_MOD] = node_info

                dic_edges_nodes_graph_by_layer[
                    "dic_cluster_total_input_idx_MOD_vs_node_info"] = dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE

    # ----------------------------------------------------------------------
    # ref layer subgraph extraction 2.
    # this time, based on ref layer subgraph connected to "usr specified" node.

    fo_log.write("\n\nlen of l_total_input_idx_mod_ref_for_subgraph_base : " + str(
        len(l_total_input_idx_mod_ref_for_subgraph_base)) + "\n")

    fo_log_key.write("\n\nlen of l_global_accession_for_subgraph_sample_user_selected : " + str(
        len(l_global_accession_for_subgraph_sample_user_selected)) + "\n")

    fo_log_key.write("\n\nlen of l_total_input_idx_mod_ref_for_subgraph_base : " + str(
        len(l_total_input_idx_mod_ref_for_subgraph_base)) + "\n")

    fo_log_key.flush()

    logger.info(f'Length of l_total_input_idx_mod_ref_for_subgraph_base: '
                f'{len(l_total_input_idx_mod_ref_for_subgraph_base)}')
    logger.info(f'Length of l_global_accession_for_subgraph_sample_user_selected'
                f'{len(l_global_accession_for_subgraph_sample_user_selected)}')

    f_show_all_ref = 0
    # doing ref layer extraction only if sugraph mode is ON for sample.
    if len(l_global_accession_for_subgraph_sample_user_selected) > 0:

        fo_log_key.write("\n performing ref layer subgraph 2 creation (user specified node)\n")
        logger.info('Performing ref layer subgraph 2 creation (user specified node)')
        # if len(l_total_input_idx_mod_ref_for_subgraph_base) > 0 :

        # this empty list will keep data for each layer to be kept
        list_dic_edges_nodes_graph_by_layer_update = []

        fo_log.write(
            "\n\n l_total_input_idx_mod_ref_for_subgraph_base:" + str(l_total_input_idx_mod_ref_for_subgraph_base))
        logger.info(f'l_total_input_idx_mod_ref_for_subgraph_base: {l_total_input_idx_mod_ref_for_subgraph_base}')

        # look for ref layer--------------------------------
        for n in range(len(list_dic_edges_nodes_graph_by_layer)):
            dic_edges_nodes_graph_by_layer = list_dic_edges_nodes_graph_by_layer[n]
            f_do_extraction = 0
            # sample layer data will be just passed.
            if dic_edges_nodes_graph_by_layer["attribute_for_layer"].startswith("sample"):
                list_dic_edges_nodes_graph_by_layer_update.append(dic_edges_nodes_graph_by_layer)

            # dealing with ref layer
            if not dic_edges_nodes_graph_by_layer["attribute_for_layer"].startswith("sample"):

                ### !!!!!!!!!!! TO DO FIX
                ref_subgraph_depth = 3
                base = []
                nx_graph = dic_edges_nodes_graph_by_layer["nx_graph"]

                # iterate node and get corresponding node to extract.
                for node in nx_graph.nodes(data=True):
                    fo_log.write("\n node" + str(node) + str(type(node)))
                    logger.info(f'Node {node}')
                    if node[0] in l_total_input_idx_mod_ref_for_subgraph_base:
                        base.append(node[0])

                fo_log.write("\n\nref sugraph base" + str(base))
                logger.info(f'Ref subgraph base {base}')

                if len(base) > 0:
                    foundset = {key for source in base for key in
                                list(nx.single_source_shortest_path(nx_graph, source,
                                                                    cutoff=ref_subgraph_depth).keys())}
                    # update and replace nx graph and other info
                    nx_graph_sub = nx_graph.subgraph(foundset)
                    dic_edges_nodes_graph_by_layer["nx_graph"] = nx_graph_sub

                    # create list of edges again.  nx.generate_edgelist  not really working ???
                    list_of_edge_for_networkx_sub = []
                    for e in nx_graph_sub.edges.data():
                        list_of_edge_for_networkx_sub.append([e[0], e[1], e[2]])
                    # update_edges
                    dic_edges_nodes_graph_by_layer["list_of_edge_for_networkx"] = list_of_edge_for_networkx_sub

                    # updating "dic_cluster_total_input_idx_MOD_vs_node_info"-----------------------------

                    # this has total input idx mod PRESENT in subgraph
                    l_total_input_idx_mod_IN_subgraph = []

                    for node in nx_graph_sub.nodes(data=True):
                        l_total_input_idx_mod_IN_subgraph.append(node[0])
                    dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE = {}
                    for cluster_total_input_idx_MOD, node_info in dic_edges_nodes_graph_by_layer[
                        "dic_cluster_total_input_idx_MOD_vs_node_info"].items():

                        if cluster_total_input_idx_MOD in l_total_input_idx_mod_IN_subgraph:
                            dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE[
                                cluster_total_input_idx_MOD] = node_info

                    dic_edges_nodes_graph_by_layer[
                        "dic_cluster_total_input_idx_MOD_vs_node_info"] = dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE

                    if len(dic_cluster_total_input_idx_MOD_vs_node_info_UPDATE) > 0:
                        fo_log.write("\nappending ref sugraph " + str(dic_edges_nodes_graph_by_layer) + "\n\n\n")
                        logger.info(f'Appending ref subgraph {dic_edges_nodes_graph_by_layer}')
                        list_dic_edges_nodes_graph_by_layer_update.append(dic_edges_nodes_graph_by_layer)

        list_dic_edges_nodes_graph_by_layer = list_dic_edges_nodes_graph_by_layer_update

    fo_log.write("\n\nmade ref subgraph" + str(list_dic_edges_nodes_graph_by_layer) + "\n\n\n")
    logger.info(f'Made ref subgraph: {list_dic_edges_nodes_graph_by_layer}')

    ######################
    # [N1t] UPDATING (SUB)GRAPH info
    # now it is possible some of graph by layer has No nodes.
    # for example, after ref compound filtration, organooxygen ref layer has NO spectra.
    #  now you want to keep "graph by layer" that actually have spectra(node)

    # output log
    fo_log_key.write("\n\n [N1t] UPDATING (SUB)GRAPH info")
    log_message = f'[N1t] UPDATING (SUB)GRAPH info\n"attribute_for_layer' \
                  f'\tlength of dic_cluster_total_input_idx_MOD_vs_node_info' \
                  f'\tlength of list_of_edge_for_networkx'
    for dic in list_dic_edges_nodes_graph_by_layer:
        fo_log_key.write(
            "\n" + dic["attribute_for_layer"] + " : " + "dic_cluster_total_input_idx_MOD_vs_node_info:" + str(
                len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])) \
            + ",    list_of_edge_for_networkx: " + str(len(dic["list_of_edge_for_networkx"])))
        log_message += f'\n{dic["attribute_for_layer"]}\t{len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])}' \
                       f'\t{len(dic["list_of_edge_for_networkx"])}'
    logger.info(log_message)

    logger.debug("[multilayer_3d_network_b1/process_3d_network_data]  update (sub)graph")
    list_dic_edges_nodes_graph_by_layer_update = []

    for dic_edges_nodes_graph_by_layer in list_dic_edges_nodes_graph_by_layer:
        if len(dic_edges_nodes_graph_by_layer["dic_cluster_total_input_idx_MOD_vs_node_info"]) > 0:
            list_dic_edges_nodes_graph_by_layer_update.append(dic_edges_nodes_graph_by_layer)

    # takeover
    list_dic_edges_nodes_graph_by_layer = list_dic_edges_nodes_graph_by_layer_update

    #####################
    ######################
    #  pick up interlayer edges that fits sugraph
    ################################
    logger.debug("[multilayer_3d_network_b1/process_3d_network_data]  [N1z] pick up interlayer edges that fits sugraph")

    l_cluster_total_input_idx_MOD_in_sample_subgraph = []
    l_cluster_total_input_idx_MOD_in_all_subgraph = []
    # foe each layer
    for dic_edges_nodes_graph_by_layer in list_dic_edges_nodes_graph_by_layer:
        # if the current layer is sample
        for cluster_total_input_idx_MOD, v in dic_edges_nodes_graph_by_layer[
            "dic_cluster_total_input_idx_MOD_vs_node_info"].items():
            l_cluster_total_input_idx_MOD_in_all_subgraph.append(cluster_total_input_idx_MOD)

    # get edges both nodes can be founc in current subgraph:
    list_of_inter_sample_ref_edge_tb_used = []
    for edge_for_networkx in list_of_edge_for_networkx_to_show_inter_sample_ref_layer:

        # if edge_for_networkx[0] in l_cluster_total_input_idx_MOD_in_sample_subgraph or  edge_for_networkx[1] in l_cluster_total_input_idx_MOD_in_sample_subgraph:
        if edge_for_networkx[0] in l_cluster_total_input_idx_MOD_in_all_subgraph and edge_for_networkx[
            1] in l_cluster_total_input_idx_MOD_in_all_subgraph:
            list_of_inter_sample_ref_edge_tb_used.append(edge_for_networkx)

    # take over
    list_of_edge_for_networkx_to_show_inter_sample_ref_layer = list_of_inter_sample_ref_edge_tb_used

    fo_log.write("list_of_edge_for_networkx_to_show_inter_sample_ref_layer \n" + str(
        list_of_edge_for_networkx_to_show_inter_sample_ref_layer) + "\n")
    logger.info(f'list_of_edge_for_networkx_to_show_inter_sample_ref_layer: '
                f'{list_of_edge_for_networkx_to_show_inter_sample_ref_layer}')

    list_of_inter_sample_edge_tb_used = []
    for edge_for_networkx in list_of_edge_for_networkx_to_show_inter_sample_layer:
        if edge_for_networkx[0] in l_cluster_total_input_idx_MOD_in_all_subgraph and edge_for_networkx[
            1] in l_cluster_total_input_idx_MOD_in_all_subgraph:
            list_of_inter_sample_edge_tb_used.append(edge_for_networkx)

    # take over
    list_of_edge_for_networkx_to_show_inter_sample_layer = list_of_inter_sample_edge_tb_used

    ### for node dictionary
    dic_cluster_total_input_idx_MOD_vs_node_info_tb_used = {}
    for cluster_total_input_idx_MOD, node_info in dic_cluster_total_input_idx_MOD_vs_node_info_inter_sample_ref_layer.items():

        if cluster_total_input_idx_MOD in l_cluster_total_input_idx_MOD_in_sample_subgraph:
            dic_cluster_total_input_idx_MOD_vs_node_info_tb_used[cluster_total_input_idx_MOD] = node_info
    # take over -----
    dic_cluster_total_input_idx_MOD_vs_node_info_inter_sample_ref_layer = dic_cluster_total_input_idx_MOD_vs_node_info_tb_used

    ########################
    # combine subgraph  and upper layer
    ########################
    logger.debug("starting [combine subgraph and upper layer]")

    ###
    # forcebly show nodes in upper layer EVEN IF the nodes are not connected to sample layer
    flag_show_all_upper = 0

    if flag_show_all_upper == 1:
        FG = nx.compose_all([FG, FG_upper])

    fo_log.write("AFTER COMPOSE FG.number_of_nodes()" + str(FG.number_of_nodes()) + "\n")
    fo_log.flush()
    logger.info(f'AFTER COMPOSE FG.number_of_nodes(): {FG.number_of_nodes()}')

    # make list of edges and nodes acutally used.-----------------------
    list_node_total_input_idx_mod_in_use = []
    for node in FG.nodes.data():
        list_node_total_input_idx_mod_in_use.append(node[0])

    l_edges_in_use = []
    for edge in FG.edges.data():
        l_edges_in_use.append((edge[0], edge[1]))

    ###################################
    # get layer attributes
    #####################################

    logger.debug("define layer")
    dic_layer_id_vs_attribute_for_layer = define_layers(dic_config, dic_cluster_total_input_idx_MOD_vs_node_info)
    logger.debug("finished to define layer")
    fo_log.write("\n dic_layer_id_vs_attribute_for_layer " + str(dic_layer_id_vs_attribute_for_layer) + "\n")
    fo_log.flush()
    logger.warning(f'dic_layer_id_vs_attribute_for_layer: {dic_layer_id_vs_attribute_for_layer}')

    ####################
    # layers in upper

    ################
    # layers in lower (samples)

    ############
    # list_dictionary     upper,   FG

    """

    ################################
    # [O] community detection
    print ("[O]staring community detection ")


    fo_log_key.write("\n\n beggining of community detection [O]")
    for dic in list_dic_edges_nodes_graph_by_layer:
        fo_log_key.write(  "\n" +  dic["attribute_for_layer"] + " : " +   "dic_cluster_total_input_idx_MOD_vs_node_info:" +  str(len( dic["dic_cluster_total_input_idx_MOD_vs_node_info"] )) \
                +   ",    list_of_edge_for_networkx: "   +  str(len(dic["list_of_edge_for_networkx"]))       )




    from networkx.readwrite import json_graph
    from networkx.algorithms import community
    import json

    str_log =  "\n\nCommunity detection==============="  + "\n"
    #####
    # [O1] first,  combine all nodes, edges ---------------------------------
    print ("  [O1] combine all nodes, edges ")
    G_all_combined = nx.Graph()

    # combine node/edges for each layers. (not inter layer edges)
    for dic_edges_nodes_graph_by_layer in list_dic_edges_nodes_graph_by_layer:
        G_all_combined = nx.compose(G_all_combined, dic_edges_nodes_graph_by_layer["nx_graph"])

    # add inter sample-ref edges
    G_all_combined.add_edges_from(list_of_edge_for_networkx_to_show_inter_sample_ref_layer)

    # add inter sample-ref edges
    G_all_combined.add_edges_from(list_of_edge_for_networkx_to_show_inter_sample_layer)



    # [O3] Perform community detection ----------------------------------------------
    print ("  [O3] Perform community detection")


    print("  [O3] community generator created ")
    n_level_community_detection = 1


    if n_level_community_detection > 0 :
        l_communities = []
        community_detection_method = "greedy_modularity_communities"


        if community_detection_method == "girvan_newman":
            communities_generator = community.girvan_newman(G_all_combined)
            # return as tuple
            for n in range(n_level_community_detection):
                l_communities = next(communities_generator)

        if community_detection_method == "greedy_modularity_communities":
            from networkx.algorithms.community import greedy_modularity_communities
            try :
                l_communities = list(greedy_modularity_communities(G_all_combined))
            except:
                pass


        print("  [O3] community detection finished")


        if len( l_communities) > 0 :

            str_log = str_log + "\n detected communities \n"
            for c in l_communities:
                str_log = str_log + "\n"  +str(c)

                # [O5]  check all the edges whether it fits community definition -----------------

            # this holds edges that do not fit within community.  in other word, inter community edges
            l_inter_community_edge_as_tuple = []

            # iterate edges----------------------
            for u, v, d in G_all_combined.edges(data=True):
                f_in_community = 0
                # iterate communities
                for commu in l_communities:
                    if all(x in commu for x in [u, v]):
                        f_in_community = 1
                # if the current edge is not inside community,
                if f_in_community == 0:
                    l_inter_community_edge_as_tuple.append((u, v, d))

            print("inter-community edges", l_inter_community_edge_as_tuple)
            str_log = str_log + "\n inter-community edges "    +    str(len(l_inter_community_edge_as_tuple)) + "\n"

            fo_log_key.write(str_log)
            fo_log_key.flush()

            # [O7]: update edge info-------------------------------------------------------------------
            print("  [O7]  update edge info")
            # for each layer.  (list_dic_edges_nodes_graph_by_layer)
            for dic_edges_nodes_graph_by_layer in list_dic_edges_nodes_graph_by_layer :
                l_inter_community_edge = []
                # list of edge within community, which supposed to be preserved.
                l_inner_community_edge = []
                # iterate edges
                for edge in dic_edges_nodes_graph_by_layer["list_of_edge_for_networkx"]:
                    f_in_community = 0
                    # check the edge is in community.
                    for commu in l_communities:
                        if all(x in commu for x in [edge[0], edge[1]]):
                            l_inner_community_edge.append((edge[0], edge[1], edge[2]))
                            f_in_community = 1

                    if f_in_community == 0:
                        l_inter_community_edge.append((edge[0], edge[1], edge[2]))

                dic_edges_nodes_graph_by_layer["list_of_edge_for_networkx"]= l_inner_community_edge
                # it seems you need to unfrozen graph
                nx_graph_unfrozen =    nx.Graph(dic_edges_nodes_graph_by_layer["nx_graph"])
                # then remove edge.
                nx_graph_unfrozen.remove_edges_from(l_inter_community_edge)
                dic_edges_nodes_graph_by_layer["nx_graph"] = nx_graph_unfrozen

                #log
                str_log =  str_log  + "layer attribute: " + dic_edges_nodes_graph_by_layer["attribute_for_layer"] + "\n"
                str_log = str_log + " num edges removed (in community detection mode: do not belong to community): " +  str(len(l_inter_community_edge))+ "\n"


            #    remove inter community edges for list_of_edge_for_networkx_to_show_inter_sample_ref_layer
            l_inner_community_edge_inter_sample_ref_layer = []
            l_inter_community_edge_inter_sample_ref_layer = []
            for edge in list_of_edge_for_networkx_to_show_inter_sample_ref_layer:
                f_in_community = 0
                # check the edge is in community.
                for commu in l_communities:
                    if all(x in commu for x in [edge[0], edge[1]]):
                        l_inner_community_edge_inter_sample_ref_layer.append((edge[0], edge[1], edge[2]))
                        f_in_community = 1

                if f_in_community == 0:
                    l_inter_community_edge_inter_sample_ref_layer.append((edge[0], edge[1], edge[2]))
            # take over
            list_of_edge_for_networkx_to_show_inter_sample_ref_layer =  l_inner_community_edge_inter_sample_ref_layer



            #   remove inter community edges  for list_of_edge_for_networkx_to_show_inter_sample_layer
            l_inner_community_edge_inter_sample_layer = []
            l_inter_community_edge_inter_sample_layer = []
            for edge in list_of_edge_for_networkx_to_show_inter_sample_layer:
                f_in_community = 0
                # check the edge is in community.
                for commu in l_communities:
                    if all(x in commu for x in [edge[0], edge[1]]):
                        l_inner_community_edge_inter_sample_layer.append((edge[0], edge[1], edge[2]))
                        f_in_community = 1

                if f_in_community == 0:
                    l_inter_community_edge_inter_sample_layer.append((edge[0], edge[1], edge[2]))
            # take over
            list_of_edge_for_networkx_to_show_inter_sample_layer =  l_inner_community_edge_inter_sample_layer




        ## export

        data1 = json_graph.node_link_data(G_all_combined)
        fp = open("(G_all_combined.json", "w")
        json.dump(data1, fp)

        print(" Finishing community detection")



    #  end community detection

    """

    ####################################
    # Main [Q] get layout, 2D
    ######################################
    fo_log_key.write("\n\n beggining of main [Q]")
    log_message = f'Beginning of main [Q]\nattribute_for_layer\tlength of dic_cluster_total_input_idx_MOD_vs_node_info' \
                  f'\tlength of list_of_edge_for_networkx'
    for dic in list_dic_edges_nodes_graph_by_layer:
        fo_log_key.write(
            "\n" + dic["attribute_for_layer"] + " : " + "dic_cluster_total_input_idx_MOD_vs_node_info:" + str(
                len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])) \
            + ",    list_of_edge_for_networkx: " + str(len(dic["list_of_edge_for_networkx"])))
        log_message += f'\n{dic["attribute_for_layer"]}\t{len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])}' \
                       f'{len(dic["list_of_edge_for_networkx"])}'
    logger.info(log_message)

    logger.debug("[multilayer_3d_network_b1/process_3d_network_data]  start getting layout, 2D ")
    ##  Layout for sample\\\\======================================================

    # first, assemble all nodes, edges from "sample tag" containing dataset-----------------------
    dic_cluster_total_input_idx_MOD_vs_node_info_all_samples = {}
    list_of_edge_for_networkx_all_samples = []

    for dic in list_dic_edges_nodes_graph_by_layer:
        # dealing with sample layer data.
        # note if you have more than 1, you have to merge and make layout for the merged datast.
        if dic["attribute_for_layer"].startswith("sample_"):
            dic_cluster_total_input_idx_MOD_vs_node_info_all_samples.update(
                dic["dic_cluster_total_input_idx_MOD_vs_node_info"])
            list_of_edge_for_networkx_all_samples = list_of_edge_for_networkx_all_samples + dic[
                "list_of_edge_for_networkx"]

    fo_x = open("logx.txt", "w")
    fo_x.write(str(len(list_of_edge_for_networkx_to_show_inter_sample_layer)) + "\n")
    logger.info(f'Length of list_of_edge_for_networkx_to_show_inter_sample_layer: '
                f'{len(list_of_edge_for_networkx_to_show_inter_sample_layer)}')

    ###############
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ### Whether you show similarity edge or only identcail edge in inter sample layer.
    f_show_only_pillar_inter_sample = 0

    list_of_edge_for_pillar = []
    count = 0
    if f_show_only_pillar_inter_sample == 1:

        for e in list_of_edge_for_networkx_to_show_inter_sample_layer:
            if abs(e[2]["delta_mz"]) < 0.1:
                list_of_edge_for_pillar.append(e)
                count = count + 1
        list_of_edge_for_networkx_to_show_inter_sample_layer = list_of_edge_for_pillar
    #

    logger.debug("[multilayer_3d_network_b1/process_3d_network_data]  performing nx spring layout for all samples ")

    # this layout ( dictionary in which key is total_input_idx_mod, value is coordintate) keeps all sample laytout and CAN BE SHARED among all samples.
    layout_2d_all_samples = nx.spring_layout(FG_all_samples, k=0.1, scale=2)

    logger.debug("[multilayer_3d_network_b1/process_3d_network_data]  finished nx spring layout for all samples")

    ###
    ## Get layout

    logger.debug("(main) get layout 2D")
    for dic in list_dic_edges_nodes_graph_by_layer:
        f_base_layer = 0
        if dic["attribute_for_layer"].startswith("sample_"):
            f_base_layer = 1

        # for sample data sets (layer). can be more than one.
        ##################
        if f_base_layer == 1:
            # NOTE !!!!  the layout for all samples "layout_2d_all_samples " contain all nodes in multiple sample layers.
            #   if multiple layers keep sharing this layout, you will bump into problem when adding Z-axis.
            # therefore you should make dedicated layout foreach layer here.
            layout_extracted = {}
            # iterate all total input idx mod for this layer
            for total_input_idx_mod, node_info in dic["dic_cluster_total_input_idx_MOD_vs_node_info"].items():
                # then feeed to new layout.
                layout_extracted[total_input_idx_mod] = layout_2d_all_samples[total_input_idx_mod]
            layout_2d = layout_extracted

        # for ref layer, you have to perform layout.
        # !!!!!!!!!!!! making ref layer shrink
        if f_base_layer == 0:
            layout_2d = nx.spring_layout(dic["nx_graph"], k=0.01, scale=2)
            # layout_2d = nx.spectral_layout(dic["nx_graph"])
        dic["layout_2d"] = layout_2d

        for idxm, o in dic["dic_cluster_total_input_idx_MOD_vs_node_info"].items():
            if not idxm in layout_2d:
                logger.debug(f"!!!!!! {idxm} not present !!!!")
                dic_config["all_data_valid"] = False
                dic_config["l_func_invalid_data"].append("read_config_file")

    logger.debug("(main) finished geting layout 2D")
    logger.debug("making traces")
    ## traces is a list plotly object, inclusing both nodes, and edges.
    ## l_mesh_annotations are plotly object.

    fo_log.write("\n\n  list_dic_edges_nodes_graph_by_layer" + "\n")
    log_message = 'list_dic_edges_nodes_graph_by_layer\nattribute_for_layer\tlist_of_edge_for_networkx'
    for dic in list_dic_edges_nodes_graph_by_layer:
        str_o = "\n"
        str_o = str_o + str(dic["attribute_for_layer"]) + " : "
        str_o = str_o + str(dic["list_of_edge_for_networkx"])
        fo_log.write(str_o)
        log_message += f'\n{dic["attribute_for_layer"]}\t{dic["list_of_edge_for_networkx"]}'

    fo_log.flush()
    logger.info(log_message)

    dic_processed_network_data = {}

    dic_processed_network_data["list_dic_edges_nodes_graph_by_layer"] = list_dic_edges_nodes_graph_by_layer
    dic_processed_network_data["dic_global_accession_vs_mass_feature"] = dic_global_accession_vs_mass_feature

    dic_processed_network_data[
        "dic_cluster_total_input_idx_MOD_vs_node_info"] = dic_cluster_total_input_idx_MOD_vs_node_info
    dic_processed_network_data[
        "list_of_edge_for_networkx_to_show_inter_sample_ref_layer"] = list_of_edge_for_networkx_to_show_inter_sample_ref_layer
    dic_processed_network_data[
        "list_of_edge_for_networkx_to_show_inter_sample_layer"] = list_of_edge_for_networkx_to_show_inter_sample_layer

    dic_processed_network_data["dic_layer_id_vs_attribute_for_layer"] = dic_layer_id_vs_attribute_for_layer

    fo_log_key.write("\n\n finishing")
    log_message = 'Finishing\nLength of each data\nattribute_for_layer\tdic_cluster_total_input_idx_MOD_vs_node_info' \
                  '\tlist_of_edge_for_networkx\tlayout_2d'

    for dic in list_dic_edges_nodes_graph_by_layer:
        fo_log_key.write(
            "\n" + dic["attribute_for_layer"] + " : " + "dic_cluster_total_input_idx_MOD_vs_node_info:" + str(
                len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])) \
            + ",    list_of_edge_for_networkx: " + str(len(dic["list_of_edge_for_networkx"])) + " : " + str(
                len(dic["layout_2d"])) + " : " + str(len(dic["layout_2d"])))
        log_message += f'{dic["attribute_for_layer"]}\t{len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])}' \
                       f'\t{len(dic["list_of_edge_for_networkx"])}\t{len(dic["layout_2d"])}'
    fo_log_key.close()
    logger.info(log_message)
    return dic_processed_network_data


