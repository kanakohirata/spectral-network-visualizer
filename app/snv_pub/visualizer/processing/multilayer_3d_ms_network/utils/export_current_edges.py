from logging import getLogger
import os

logger = getLogger(__name__)


def export_current_edges(list_dic_edges_nodes_graph_by_layer,
                         list_of_edge_for_networkx_to_show_inter_sample_ref_layer,
                         list_of_edge_for_networkx_to_show_inter_sample_layer):
    fo_edge_info = open(os.path.join(os.getcwd(), 'logs', "current_edge_info.tsv"), "w", encoding="utf-8")

    ########################
    # write edge info
    ##########################

    # note CLUSTERID1 actually means CLUSTERIDX  "X"

    fo_edge_info.write("X_TOTAL_INPUT_IDX_MOD"
                       + "\t" + "Y_TOTAL_INPUT_IDX_MOD"
                       + "\t" + "SPEC_SIM_SCORE"
                       + "\t" + "DELTA_MZ"
                       + "\t" + "EDGE_TYPE"
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
