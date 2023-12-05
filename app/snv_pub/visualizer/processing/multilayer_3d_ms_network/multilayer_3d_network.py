# TO DO
#  "select_nodes_based_on_product_mz_required"
#           currently only examine representative spec.
import copy
from logging import getLogger
import networkx as nx
import os.path
import plotly.graph_objs as go
import sys
from . import multilayer_3d_mesh_functsions as m3d_mesh
from . import multilayer_3d_rescale_functions as m3d_rescale
from . import suspect_compound
from .filter import (remove_nodes_based_on_total_input_idx,
                     select_nodes_based_on_keyword,
                     select_nodes_based_on_mass_defect,
                     select_nodes_based_on_prec_mz,
                     select_nodes_based_on_product_mz_required)
from .my_parser.cluster_attribute_parser import read_cluster_attribute
from .my_parser.compound_info_parser import add_external_cmpd_info
from .my_parser.edge_info_parser import read_edge_info
from .my_parser.feature_table_parser import read_feature_table
from .networking import (create_networkx_graph,
                         create_quantitative_subgraph,
                         create_sample_networkx_graph,
                         define_layers,
                         extract_ref_subgraph_based_on_total_input_idx,
                         extract_ref_subgraph_in_external_compounds,
                         extract_subgraph_based_on_sample_global_accession,
                         locate_nodes_to_layers_and_update_edges,
                         make_edges_and_nodes_inner_layer,
                         make_edges_and_nodes_inter_layer,
                         make_list_of_edge_for_networkx,
                         perform_community_detection,
                         update_edges)
from .utils import (add_color_to_t3db_compound,
                    add_suspect_compound_info)

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
    logger.info('Start read_data_for_multilayer_3d_network()')
    logger.debug('[C]reading cluster attribute')

    fo_log = open('log_read_data_for_multilayer_3d_network.txt', "w")

    #####################
    # read input data

    ##########################################################################################
    # section D
    # [D] read spec cluster info
    ##########################################################################################
    logger.debug('* main [D]reading cluster attribute')
    dic_cluster_total_input_idx_vs_cluster_info = read_cluster_attribute(dic_config['filename_cluster_info'])

    # reading external info
    ###########################################
    # [F] Read external info for node/cluster
    ###########################################
    logger.warning('* main [D]reading external compound info')
    add_external_cmpd_info(dic_cluster_total_input_idx_vs_cluster_info, dic_config['foldername_ext_cmpd_info'])
    logger.warning('* main [D]finished reading external compound info')

    # make original copy
    dic_cluster_total_input_idx_vs_cluster_info_original = copy.deepcopy(dic_cluster_total_input_idx_vs_cluster_info)

    logger.warning(f'* main[D]len dic_cluster_total_input_idx_vs_cluster_info '
                   f'{len(dic_cluster_total_input_idx_vs_cluster_info)}')

    ##########################################################################################
    # section [H]
    # [H] read EDGES info
    ##########################################################################################
    l_edges, list_edge_info = read_edge_info(dic_config['filename_edge_info'], dic_config['score_threshold'])

    fo_log.write('\n\n[F1]   list_edge_info_original')
    log_message = f'[F1]   list_edge_info_original' \
                  f'\nspec_cluster_x_global_accession <-> spec_cluster_y_global_accession'
    for e in list_edge_info:
        fo_log.write(f'\n{e["spec_cluster_x_global_accession"]} <-> {e["spec_cluster_y_global_accession"]}')
        log_message += f'\n{e["spec_cluster_x_global_accession"]} <-> {e["spec_cluster_y_global_accession"]}'

    logger.info(log_message)

    #############################
    # [J] read feature table
    #############################
    logger.debug('[J] reading feature table')
    #  --------------MODIFIED20220714
    # read_feature table.
    # TODO: Check feature_table_parser
    dic_global_accession_vs_mass_feature_original = read_feature_table(dic_config)
    dic_global_accession_vs_mass_feature = read_feature_table(dic_config)

    # locate nodes to layers
    _, list_edge_info, _ = \
        locate_nodes_to_layers_and_update_edges(dic_config, dic_cluster_total_input_idx_vs_cluster_info, list_edge_info)

    logger.debug('returned from  locate_nodes_to_layers_and_update_edges')

    # make list_of_edges_for_networkx
    logger.debug('Main  [make list_of_edges_for_networkx]')
    list_of_edge_for_networkx = make_list_of_edge_for_networkx(list_edge_info)
    _, list_node_total_input_idx_mod_in_use = threshold_edges(dic_config['score_threshold'],
                                                              list_of_edge_for_networkx)

    dic_source_data = {
        'dic_config': dic_config,
        'list_edge_info': list_edge_info,
        'dic_cluster_total_input_idx_vs_cluster_info_original': dic_cluster_total_input_idx_vs_cluster_info_original,
        'dic_cluster_total_input_idx_vs_cluster_info': dic_cluster_total_input_idx_vs_cluster_info,
        'list_of_edge_for_networkx': list_of_edge_for_networkx,
        'list_node_total_input_idx_mod_in_use': list_node_total_input_idx_mod_in_use,
        'dic_global_accession_vs_mass_feature_original': dic_global_accession_vs_mass_feature_original,
        'dic_global_accession_vs_mass_feature': dic_global_accession_vs_mass_feature
    }

    fo_log.close()
    return dic_source_data


def process_3d_network_data(dic_source_data, dic_config):
    dic_global_accession_vs_mass_feature = dic_source_data["dic_global_accession_vs_mass_feature"]

    logger.debug("[create_multilayer_3d_network_data]processing")

    fo_log = open("log_create_multilayer_3d_network_data.txt", "w")
    fo_log_key = open("log_KEY_multilayer_3d_network.txt", "w")

    # Read some external files
    # suspect mz file
    dic_config["dic_filename_vs_dic_suspect_cmpd_vs_mz"] = suspect_compound.read_file_suspect_mz(dic_config)

    # [C1] --------------------------------------
    # deepcopy input idx and cluster info
    dic_cluster_total_input_idx_vs_cluster_info_original = copy.deepcopy(
        dic_source_data["dic_cluster_total_input_idx_vs_cluster_info"])

    fo_log.write("\n\n\n another fig update")
    fo_log.write(f"\n dic_cluster_total_input_idx_vs_cluster_info_original: "
                 f"{len(dic_cluster_total_input_idx_vs_cluster_info_original)}")

    logger.info('Another fig update')
    logger.info(f'Length of dic_cluster_total_input_idx_vs_cluster_info_original: '
                f'{len(dic_cluster_total_input_idx_vs_cluster_info_original)}')

    # Add color to T3DB compounds
    add_color_to_t3db_compound(dic_config['color_toxic_compound'], dic_cluster_total_input_idx_vs_cluster_info_original)

    # select nodes based on keyword
    dic_cluster_total_input_idx_vs_cluster_info =\
        select_nodes_based_on_keyword(dic_config['filter_select_category'],
                                      dic_config['filter_select_keyword'],
                                      dic_cluster_total_input_idx_vs_cluster_info_original)
    logger.debug(f"\n select nodes based on keyword : len dic_cluster_total_input_idx_vs_cluster_info: "
                 f"{str(len(dic_cluster_total_input_idx_vs_cluster_info))}")

    #  [C2]
    # select nodes based on prec mz
    logger.debug(f"[C2]: len dic_cluster_total_input_idx_vs_cluster_info "
                 f"{len(dic_cluster_total_input_idx_vs_cluster_info)}")
    dic_cluster_total_input_idx_vs_cluster_info =\
        select_nodes_based_on_prec_mz(dic_config["mass_lower_limit"],
                                      dic_config["mass_higher_limit"],
                                      dic_cluster_total_input_idx_vs_cluster_info)

    fo_log.write(f"\n after selecting nodes, dic_cluster_total_input_idx_vs_cluster_info: "
                 f"{len(dic_cluster_total_input_idx_vs_cluster_info)}")
    fo_log.flush()
    logger.info(f'After selecting nodes based on precursor m/z, length of dic_cluster_total_input_idx_vs_cluster_info: '
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
        dic_cluster_total_input_idx_vs_cluster_info =\
            select_nodes_based_on_mass_defect(dic_config["mz_tol"],
                                              dic_config["list_mass_defect"],
                                              dic_cluster_total_input_idx_vs_cluster_info)
    # ------------------------------- MODIFIED20220714---------------

    # MODIFIED20220719---------------
    # [C4] select node based on product ion
    if len(dic_config["list_product_mz_required"]) > 0:
        dic_cluster_total_input_idx_vs_cluster_info =\
            select_nodes_based_on_product_mz_required(dic_config["list_product_mz_required"],
                                                      dic_config["mz_tolerance_for_fragment"],
                                                      dic_cluster_total_input_idx_vs_cluster_info)

        logger.debug(f" l of dic after product filer: {len(dic_cluster_total_input_idx_vs_cluster_info)}")
    # ------------------------------MODIFIED20220719---------------

    # remove nodes based on user defined input idx ----------------------------------
    logger.debug(f"dic_config[l_total_input_idx_to_remove] {str(dic_config['l_total_input_idx_to_remove'])}")
    dic_cluster_total_input_idx_vs_cluster_info =\
        remove_nodes_based_on_total_input_idx(dic_config["l_total_input_idx_to_remove"],
                                              dic_cluster_total_input_idx_vs_cluster_info)

    ###############################
    # [F1]
    # update node and edges. now node id is "MOD" and make node info
    # now edge has attribute "inner_sample_layer" etc.

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
    for e in list_edge_info_original:
        fo_log.write("\n" + e["spec_cluster_x_global_accession"] + " <-> " + e["spec_cluster_y_global_accession"])
        log_message += f'\n{e["spec_cluster_x_global_accession"]} <-> {e["spec_cluster_y_global_accession"]}'
    logger.info(log_message)

    (dic_cluster_total_input_idx_MOD_vs_node_info,
     list_edge_info,
     dic_cluster_id_vs_l_cluster_total_input_idx_MOD) = \
        locate_nodes_to_layers_and_update_edges(dic_config,
                                                dic_cluster_total_input_idx_vs_cluster_info,
                                                list_edge_info_original)

    # Logging ------------------------------------------------------------------------
    fo_log.write("\n\ndic_cluster_total_input_idx_vs_cluster_info")
    for idx, cl in dic_cluster_total_input_idx_vs_cluster_info.items():
        fo_log.write("\n" + str(idx))
    fo_log.flush()
    logger.info(f'Keys of dic_cluster_total_input_idx_vs_cluster_info: '
                f'{dic_cluster_total_input_idx_vs_cluster_info.keys()}')

    logger.debug("[F1] finished  locate_nodes_to_layers_and_update_edges.")
    fo_log.write(f"\n after locate_nodes_to_layers_and_update_edges: {len(list_edge_info)}")
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
    # --------------------------------------------------------------------------------

    # what you get is [node_id_X, node_id_Y, dictionary_for_attribute]
    list_of_edge_for_networkx = make_list_of_edge_for_networkx(list_edge_info)
    fo_log.write("\n\n [F1]Show raw edges----------")
    log_message = '[F1] Show raw edges'
    for e in list_of_edge_for_networkx:
        fo_log.write("\n" + str(e))
        log_message += f'{e}'
    logger.info(log_message)

    # MODIFIED20220714-------------------
    # [F4] add info based on with suspect list.--------------------------------------------------
    logger.debug("starting suspect list matching")

    fo_log_key.write("\n number of suspect files: " + str(dic_config["dic_filename_vs_dic_suspect_cmpd_vs_mz"].keys()))
    logger.info(f'Number of suspect files: {len(dic_config["dic_filename_vs_dic_suspect_cmpd_vs_mz"].keys())}')
    for k, v in dic_config["dic_filename_vs_dic_suspect_cmpd_vs_mz"].items():
        fo_log_key.write("\n   " + k + " : " + str(len(v)))
        logger.info(f'Number of suspect compounds in {k}: {len(v)}')

    dic_cluster_total_input_idx_MOD_vs_node_info =\
        add_suspect_compound_info(dic_config["dic_filename_vs_dic_suspect_cmpd_vs_mz"],
                                  dic_config["mz_tol"],
                                  dic_cluster_total_input_idx_MOD_vs_node_info)
    logger.debug("suspect match finished")

    fo_log_key.write(f"\n\n Main F4 (suspect list) len of dic_cluster_total_input_idx_MOD_vs_node_info: "
                     f"{len(dic_cluster_total_input_idx_MOD_vs_node_info)}\n")
    fo_log_key.flush()
    logger.info(f'Main F4 (suspect list) length of dic_cluster_total_input_idx_MOD_vs_node_info: '
                f'{len(dic_cluster_total_input_idx_MOD_vs_node_info)}')

    # TODO: T3DB matching (K. Hirata, 20221214)

    # [H1] -----------------------------------------------
    # thresholding edges -----------------------
    #  list_of_edge_for_networkx_to_show  is list of list[node_X, node_Y, dictionary{spec_sim_score, delta_mz, edge_type}]
    #  the edge type can be "inner_ref_layer", "inter_sample_ref_layer", "inner_sample_layer"
    logger.debug("[H1] Thresholding edges")

    list_of_edge_for_networkx_to_show, list_node_total_input_idx_mod_in_use =\
        threshold_edges(dic_config['score_threshold'], list_of_edge_for_networkx)

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
    # make "Graph by layer".
    # also isolate inter layer network (edges + layers))
    ########################################
    logger.debug("[J1] making graph by layer")
    list_dic_edges_nodes_graph_by_layer =\
        make_edges_and_nodes_inner_layer(list_attribute_for_layer,
                                         list_of_edge_for_networkx_to_show,
                                         dic_cluster_total_input_idx_MOD_vs_node_info,
                                         fo_log)

    fo_log.write(f"\n\n\n  list_dic_edges_nodes_graph_by_layer {len(list_dic_edges_nodes_graph_by_layer)}\n\n\n")
    logger.info(f'length of list_dic_edges_nodes_graph_by_layer: {len(list_dic_edges_nodes_graph_by_layer)}')
    logger.debug(f"len list_dic_edges_nodes_graph_by_layer: {len(list_dic_edges_nodes_graph_by_layer)}")

    fo_log_key.write("\n\n*end of main [J1] (making list_dic_edges_nodes_graph_by_layer )")
    log_message = 'End of main [J1] (making list_dic_edges_nodes_graph_by_layer )' \
                  '\nattribute_for_layer\tlength of dic_cluster_total_input_idx_MOD_vs_node_info' \
                  '\tlength of list_of_edge_for_networkx'
    for dic in list_dic_edges_nodes_graph_by_layer:
        fo_log_key.write(
            f'\n{dic["attribute_for_layer"]} : dic_cluster_total_input_idx_MOD_vs_node_info:'
            f'{len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])},    '
            f'list_of_edge_for_networkx: {len(dic["list_of_edge_for_networkx"])}')

        log_message += f'\n{dic["attribute_for_layer"]}\t{len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])}' \
                       f'\t{len(dic["list_of_edge_for_networkx"])}'
    logger.info(log_message)

    #################################################################
    # [L1]
    # making dataset for "INTER  SAMPLE VS REF  LAYER"  and "INTER SAMPLE" dataset.--------------------------------
    fo_log_key.write("\n\n [L1 ]making dataset for INTER  SAMPLE VS REF  LAYER  and INTER SAMPLE dataset")
    logger.debug("[L1] making dataset for inter sample vs ref layer  and inter sample dataset")
    logger.info("[L1] making dataset for inter sample vs ref layer  and inter sample dataset")

    (list_of_edge_for_networkx_to_show_inter_sample_ref_layer,
     list_of_edge_for_networkx_to_show_inter_sample_layer,
     dic_cluster_total_input_idx_MOD_vs_node_info_inter_sample_ref_layer,
     dic_cluster_total_input_idx_MOD_vs_node_info_inter_sample_layer,
     list_total_input_idx_MOD_inter_sample_ref_layer,
     list_total_input_idx_MOD_inter_sample_layer) =\
        make_edges_and_nodes_inter_layer(list_of_edge_for_networkx_to_show,
                                         dic_cluster_total_input_idx_MOD_vs_node_info,
                                         fo_log_key)

    ########################################################
    # M1
    # [M1] crate networkx graph for each layer
    logger.debug("[multilayer_3d_network_b1/process_3d_network_data]  [M1] crate networkx graph for each layer ")
    for dic_edges_nodes_graph_by_layer in list_dic_edges_nodes_graph_by_layer:
        fo_log.write(f'\n creating FG, whoing edges: {dic_edges_nodes_graph_by_layer["attribute_for_layer"]}'
                     f'{dic_edges_nodes_graph_by_layer["list_of_edge_for_networkx"]}')
        logger.info(f"Creating FG for '{dic_edges_nodes_graph_by_layer['attribute_for_layer']}, "
                    f"{dic_edges_nodes_graph_by_layer['list_of_edge_for_networkx']}'")
        FG = create_networkx_graph(dic_edges_nodes_graph_by_layer["dic_cluster_total_input_idx_MOD_vs_node_info"],
                                   dic_edges_nodes_graph_by_layer["list_of_edge_for_networkx"])
        dic_edges_nodes_graph_by_layer["nx_graph"] = FG

        fo_log.write(f'\nnumber of network (layers )  len of (dic_edges_nodes_graph_by_layer):'
                     f'{len(dic_edges_nodes_graph_by_layer)}\n')
        fo_log.flush()
        logger.info(f'Number of network layers (length of dic_edges_nodes_graph_by_layer): '
                    f'{len(dic_edges_nodes_graph_by_layer)}')

    ################################
    # [O] community detection
    logger.debug("[O]staring community detection ")

    fo_log_key.write("\n\n beginning of community detection [O]")
    log_message = 'Beginning of community detection [O]' \
                  '\nattribute_for_layer\tlength of dic_cluster_total_input_idx_MOD_vs_node_info' \
                  '\tlength of list_of_edge_for_networkx'
    for dic in list_dic_edges_nodes_graph_by_layer:
        fo_log_key.write(
            f'\n{dic["attribute_for_layer"]} : '
            f'dic_cluster_total_input_idx_MOD_vs_node_info: {len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])}'
            f',    list_of_edge_for_networkx: {len(dic["list_of_edge_for_networkx"])}')
        log_message += f'\n{dic["attribute_for_layer"]}\t{len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])}' \
                       f'\t{len(dic["list_of_edge_for_networkx"])}'
    logger.info(log_message)

    n_level_community_detection = dic_config["n_level_community_detection"]

    fo_log_key.write("\n level of community detection (0: no community detection): " + str(n_level_community_detection))
    logger.info(f'Level of community detection (0: no community detection): {n_level_community_detection}')

    if n_level_community_detection > 0:
        (list_dic_edges_nodes_graph_by_layer,
         list_of_edge_for_networkx_to_show_inter_sample_ref_layer,
         list_of_edge_for_networkx_to_show_inter_sample_layer,) =\
            perform_community_detection(n_level_community_detection,
                                        list_dic_edges_nodes_graph_by_layer,
                                        list_of_edge_for_networkx_to_show_inter_sample_ref_layer,
                                        list_of_edge_for_networkx_to_show_inter_sample_layer,
                                        fo_log_key)

    logger.debug("[O3] community detection finished")
    #  end community detection

    fo_log_key.write("\n\n\n after community detection")
    fo_log_key.write("\n list_of_edge_for_networkx_to_show_inter_sample_ref_layer : "
                     f"{len(list_of_edge_for_networkx_to_show_inter_sample_ref_layer)}")
    logger.info(f'After community detection, length of list_of_edge_for_networkx_to_show_inter_sample_ref_layer: '
                f'{len(list_of_edge_for_networkx_to_show_inter_sample_ref_layer)}')

    ############################################################
    # [N1]  SUBGRAPH
    ############################################################
    # SUBGRAPH  (sample)
    #  here you are making list of global_accession of which mass feature show significance quantitative shift.

    fo_log_key.write("\n\n [N1]  SUBGRAPH")
    log_message = f'[N1] SUBGRAPH\nattribute_for_layer\tlength of dic_cluster_total_input_idx_MOD_vs_node_info' \
                  f'\tlist_of_edge_for_networkx'
    for dic in list_dic_edges_nodes_graph_by_layer:
        fo_log_key.write(
            f'\n{dic["attribute_for_layer"]} : '
            f'dic_cluster_total_input_idx_MOD_vs_node_info:{len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])}'
            f',    list_of_edge_for_networkx: {len(dic["list_of_edge_for_networkx"])}')
        log_message += f'\n{dic["attribute_for_layer"]}\t{len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])}' \
                       f'\t{len(dic["list_of_edge_for_networkx"])}'
    logger.info(log_message)
    logger.debug("[multilayer_3d_network_b1/process_3d_network_data]  [N1] creating [SUBGRAPH  (sample)] ")

    fo_log_key.write(f'\n\nconfig_o.stat_value: {dic_config["stat_value"]}   {type(dic_config["stat_value"])}\n')
    fo_log_key.write(f'config_o.quant_value: {dic_config["quant_value"]}\n')
    fo_log_key.write(f'topn: {dic_config["subgraph_num_core_nodes"]}\n')
    fo_log_key.write(f'config_o.subgraph_depth: {dic_config["subgraph_depth"]}\n')

    logger.info(f'config_o.stat_value: {dic_config["stat_value"]}')
    logger.info(f'config_o.quant_value: {dic_config["quant_value"]}')
    logger.info(f'config_o.subgraph_num_core_nodes: {dic_config["subgraph_num_core_nodes"]}')
    logger.info(f'config_o.subgraph_depth: {dic_config["subgraph_depth"]}')
    logger.info(f'config_o.quant_polar: {dic_config["quant_polar"]}')

    ##############################
    # [N1a] Quantitative subgraph
    ##############################
    if dic_global_accession_vs_mass_feature:
        list_dic_edges_nodes_graph_by_layer, l_total_input_idx_mod_sample =\
            create_quantitative_subgraph(dic_config["subgraph_type"],
                                         dic_config["subgraph_depth"],
                                         dic_config["subgraph_num_core_nodes"],
                                         dic_config["quant_polar"],
                                         dic_config["quant_value_type"],
                                         dic_config["quant_value"],
                                         dic_config["stat_value"],
                                         dic_global_accession_vs_mass_feature,
                                         list_dic_edges_nodes_graph_by_layer,
                                         fo_log_key)

    else:
        l_total_input_idx_mod_sample = []

    ##########################################
    # [N1b] user select node-based subgraph
    ##########################################
    fo_log_key.write("\n\n [N1b] user select node-based subgraph")
    log_message = f'[N1b] user select node-based subgraph\nattribute_for_layer' \
                  f'\tlength of dic_cluster_total_input_idx_MOD_vs_node_info\tlength of list_of_edge_for_networkx'
    for dic in list_dic_edges_nodes_graph_by_layer:
        fo_log_key.write(
            f'\n{dic["attribute_for_layer"]} : '
            f'dic_cluster_total_input_idx_MOD_vs_node_info:{len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])}'
            f',    list_of_edge_for_networkx: {len(dic["list_of_edge_for_networkx"])}')
        log_message += f'\n{dic["attribute_for_layer"]}\t{len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])}' \
                       f'\t{len(dic["list_of_edge_for_networkx"])}'
    logger.info(log_message)

    logger.debug("[multilayer_3d_network_b1/process_3d_network_data]  [N1b] user select node-based subgraph")

    # combine all sample data
    FG_all_samples = create_sample_networkx_graph(list_dic_edges_nodes_graph_by_layer,
                                                  list_of_edge_for_networkx_to_show_inter_sample_layer,
                                                  fo_log)
    
    # create subgraph for sample-combined dataset -----------------------------------------------------
    l_global_accession_for_subgraph_sample_user_selected = dic_config["l_global_accession_for_node_select_subgraph"]

    if len(l_global_accession_for_subgraph_sample_user_selected) > 0:
        (list_dic_edges_nodes_graph_by_layer,
         FG_all_samples,
         l_total_input_idx_mod_sample_new) = \
            extract_subgraph_based_on_sample_global_accession(l_global_accession_for_subgraph_sample_user_selected,
                                                              dic_config["node_select_subgraph_depth"],
                                                              list_dic_edges_nodes_graph_by_layer,
                                                              FG_all_samples,
                                                              fo_log)

        l_total_input_idx_mod_sample += l_total_input_idx_mod_sample_new

    fo_log.write("\n\n len of list_of_edge_for_networkx_to_show_inter_sample_ref_layer:"
                 f"{len(list_of_edge_for_networkx_to_show_inter_sample_ref_layer)}")
    fo_log.write("\n\n list_of_edge_for_networkx_to_show_inter_sample_ref_layer"
                 f"{list_of_edge_for_networkx_to_show_inter_sample_ref_layer}")
    logger.info(f'Length of list_of_edge_for_networkx_to_show_inter_sample_ref_layer: '
                f'{len(list_of_edge_for_networkx_to_show_inter_sample_ref_layer)}')
    logger.info(f'list_of_edge_for_networkx_to_show_inter_sample_ref_layer: '
                f'{list_of_edge_for_networkx_to_show_inter_sample_ref_layer}')
    fo_log_key.write("\n\n len of list_of_edge_for_networkx_to_show_inter_sample_ref_layer: "
                     f"{len(list_of_edge_for_networkx_to_show_inter_sample_ref_layer)}")
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
            f'\n{dic["attribute_for_layer"]} : '
            f'dic_cluster_total_input_idx_MOD_vs_node_info:{len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])}'
            f',    list_of_edge_for_networkx: {len(dic["list_of_edge_for_networkx"])}')
        log_message += f'\n{dic["attribute_for_layer"]}\t{len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])}' \
                       f'\t{len(dic["list_of_edge_for_networkx"])}'
    logger.info(log_message)

    logger.debug("[multilayer_3d_network_b1/process_3d_network_data]  [N1r create subgraph for ref layer")
    # !!!!!!!!!!!!!!!! NOTE !!!!!!!!!!!!!
    #  l_total_input_idx_mod_sample is already only contain use-selected nodes, is user use that function. !
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

    fo_log.write(f"\n\nl_total_input_idx_mod_ref_for_subgraph_base\n{l_total_input_idx_mod_ref_for_subgraph_base}")
    fo_log_key.write("\n len of list_of_edge_for_networkx_to_show_inter_sample_ref_layer:"
                     f"{len(list_of_edge_for_networkx_to_show_inter_sample_ref_layer)}\n")
    fo_log_key.write("\n len of l_total_input_idx_mod_ref_for_subgraph_base:"
                     f"{len(l_total_input_idx_mod_ref_for_subgraph_base)}\n")
    fo_log_key.write(f"\n len of l_total_input_idx_mod_sample:{len(l_total_input_idx_mod_sample)}\n")
    fo_log_key.flush()

    logger.info(f'l_total_input_idx_mod_ref_for_subgraph_base: {l_total_input_idx_mod_ref_for_subgraph_base}')
    logger.info(f'Length of list_of_edge_for_networkx_to_show_inter_sample_ref_layer: '
                f'{len(list_of_edge_for_networkx_to_show_inter_sample_ref_layer)}')
    logger.info(f'Length of l_total_input_idx_mod_ref_for_subgraph_base: '
                f'{len(l_total_input_idx_mod_ref_for_subgraph_base)}')
    logger.info(f'Length of l_total_input_idx_mod_sample: {len(l_total_input_idx_mod_sample)}')

    #######################################################
    # preserve nodes, edges only connected to sample layer
    depth_ref_preserve_interlayer = 10
    f_preserve_only_ref_node_inter_layer_connected = 1

    if f_preserve_only_ref_node_inter_layer_connected > 0:
        list_dic_edges_nodes_graph_by_layer =\
            extract_ref_subgraph_based_on_total_input_idx(list_dic_edges_nodes_graph_by_layer,
                                                          list_total_input_idx_MOD_inter_sample_ref_layer,
                                                          depth_ref_preserve_interlayer)

    #########################################################
    # This is for extracting compounds and related, BASED ON EXTERNAL CMPD INFO
    # !!!!!!!!!!!!!!!! tempo
    ref_subgraph_depth = 3
    # if "external compound info" is selected to make subgraph
    if dic_config.get("ref_filter_ext") == "ext_cmpd_info":
        list_dic_edges_nodes_graph_by_layer =\
            extract_ref_subgraph_in_external_compounds(list_dic_edges_nodes_graph_by_layer, ref_subgraph_depth)

    # ----------------------------------------------------------------------
    # ref layer subgraph extraction 2.
    # this time, based on ref layer subgraph connected to "usr specified" node.

    fo_log.write("\n\nlen of l_total_input_idx_mod_ref_for_subgraph_base : "
                 f"{len(l_total_input_idx_mod_ref_for_subgraph_base)}\n")
    fo_log_key.write("\n\nlen of l_global_accession_for_subgraph_sample_user_selected : "
                     f"{len(l_global_accession_for_subgraph_sample_user_selected)}\n")
    fo_log_key.write("\n\nlen of l_total_input_idx_mod_ref_for_subgraph_base : "
                     f"{len(l_total_input_idx_mod_ref_for_subgraph_base)}\n")
    fo_log_key.flush()

    logger.info(f'Length of l_total_input_idx_mod_ref_for_subgraph_base: '
                f'{len(l_total_input_idx_mod_ref_for_subgraph_base)}')
    logger.info(f'Length of l_global_accession_for_subgraph_sample_user_selected'
                f'{len(l_global_accession_for_subgraph_sample_user_selected)}')

    f_show_all_ref = 0
    # doing ref layer extraction only if subgraph mode is ON for sample.
    if len(l_global_accession_for_subgraph_sample_user_selected) > 0:

        fo_log_key.write("\n performing ref layer subgraph 2 creation (user specified node)\n")
        logger.info('Performing ref layer subgraph 2 creation (user specified node)')

        fo_log.write(
            f"\n\n l_total_input_idx_mod_ref_for_subgraph_base:{l_total_input_idx_mod_ref_for_subgraph_base}")
        logger.info(f'l_total_input_idx_mod_ref_for_subgraph_base: {l_total_input_idx_mod_ref_for_subgraph_base}')

        # look for ref layer--------------------------------
        list_dic_edges_nodes_graph_by_layer =\
            extract_ref_subgraph_based_on_total_input_idx(list_dic_edges_nodes_graph_by_layer,
                                                          l_total_input_idx_mod_ref_for_subgraph_base,
                                                          ref_subgraph_depth)

    fo_log.write("\n\nmade ref subgraph" + str(list_dic_edges_nodes_graph_by_layer) + "\n\n\n")
    logger.info(f'Made ref subgraph: {list_dic_edges_nodes_graph_by_layer}')

    ######################################################################################
    # [N1t] UPDATING (SUB)GRAPH info
    # now it is possible some graph by layer has No nodes.
    # for example, after ref compound filtration, organooxygen ref layer has NO spectra.
    # now you want to keep "graph by layer" that actually have spectra(node)

    # output log
    fo_log_key.write("\n\n [N1t] UPDATING (SUB)GRAPH info")
    log_message = f'[N1t] UPDATING (SUB)GRAPH info\n"attribute_for_layer' \
                  f'\tlength of dic_cluster_total_input_idx_MOD_vs_node_info' \
                  f'\tlength of list_of_edge_for_networkx'
    for dic in list_dic_edges_nodes_graph_by_layer:
        fo_log_key.write(
            f'\n{dic["attribute_for_layer"]} : '
            f'dic_cluster_total_input_idx_MOD_vs_node_info:{len(dic["dic_cluster_total_input_idx_MOD_vs_node_info"])}'
            f',    list_of_edge_for_networkx: {len(dic["list_of_edge_for_networkx"])}')
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

    ###############################################
    #  pick up interlayer edges that fits subgraph
    ###############################################
    logger.debug("[multilayer_3d_network_b1/process_3d_network_data]  [N1z] pick up interlayer edges that fits subgraph")

    list_of_edge_for_networkx_to_show_inter_sample_ref_layer, list_of_edge_for_networkx_to_show_inter_sample_layer =\
        update_edges(list_dic_edges_nodes_graph_by_layer,
                     list_of_edge_for_networkx_to_show_inter_sample_ref_layer,
                     list_of_edge_for_networkx_to_show_inter_sample_layer)

    fo_log.write("list_of_edge_for_networkx_to_show_inter_sample_ref_layer \n"
                 f"{list_of_edge_for_networkx_to_show_inter_sample_ref_layer}\n")
    logger.info(f'list_of_edge_for_networkx_to_show_inter_sample_ref_layer: '
                f'{list_of_edge_for_networkx_to_show_inter_sample_ref_layer}')

    ########################
    # combine subgraph  and upper layer
    ########################
    logger.debug("starting [combine subgraph and upper layer]")

    ###################################
    # get layer attributes
    #####################################

    logger.debug("define layer")
    dic_layer_id_vs_attribute_for_layer = define_layers(dic_config["str_key_attribute_to_base_layer"],
                                                        dic_cluster_total_input_idx_MOD_vs_node_info)
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

    ####################################
    # Main [Q] get layout, 2D
    ######################################
    fo_log_key.write("\n\n beginning of main [Q]")
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


