from logging import getLogger
import networkx as nx


logger = getLogger(__name__)


def add_2d_layout(list_dic_edges_nodes_graph_by_layer, dic_config, FG_all_samples):
    """

    Parameters
    ----------
    list_dic_edges_nodes_graph_by_layer
    dic_config
    FG_all_samples

    Returns
    -------
    list_dic_edges_nodes_graph_by_layer_new : list
     Updated list_dic_edges_nodes_graph_by_layer
    """
    # Update list_dic_edges_nodes_graph_by_layer
    list_dic_edges_nodes_graph_by_layer_new = []

    # this layout ( dictionary in which key is total_input_idx_mod, value is coordinate) keeps all sample layout and CAN BE SHARED among all samples.
    layout_2d_all_samples = nx.spring_layout(FG_all_samples, k=0.1, scale=2)
    logger.debug("[multilayer_3d_network_b1/process_3d_network_data]  finished nx spring layout for all samples")

    ############
    # Get layout

    logger.debug("(main) get layout 2D")
    for dic in list_dic_edges_nodes_graph_by_layer:
        f_base_layer = 0
        if dic["attribute_for_layer"].startswith("sample_"):
            f_base_layer = 1

        # for sample data sets (layer). can be more than one.
        ##################
        if f_base_layer == 1:
            # NOTE !!!!  the layout for all samples "layout_2d_all_samples" contain all nodes in multiple sample layers.
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
        else:
            layout_2d = nx.spring_layout(dic["nx_graph"], k=0.01, scale=2)

        dic_new = {}
        for k, v in dic.items():
            dic_new[k] = v
        dic_new["layout_2d"] = layout_2d

        for idxm, o in dic_new["dic_cluster_total_input_idx_MOD_vs_node_info"].items():
            if idxm not in layout_2d:
                logger.warning(f"!!!!!! {idxm} not present !!!!")
                dic_config["all_data_valid"] = False
                dic_config["l_func_invalid_data"].append("read_config_file")

        list_dic_edges_nodes_graph_by_layer_new.append(dic_new)

    return list_dic_edges_nodes_graph_by_layer_new
