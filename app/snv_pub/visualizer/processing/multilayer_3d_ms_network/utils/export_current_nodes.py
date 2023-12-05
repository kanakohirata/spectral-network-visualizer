from logging import getLogger
import os

logger = getLogger(__name__)


def export_current_nodes(list_dic_edges_nodes_graph_by_layer):
    fo_node_info = open(os.path.join(os.getcwd(), 'logs', "current_node_info.tsv"), "w", encoding="utf-8")

    str_fo = "TOTAL_INPUT_IDX_MOD" \
             + "\t" + "LAYER" \
             + "\t" + "CLISTER_ID" \
             + "\t" + "GLOBAL_ACCESSION" \
             + "\t" + "DATASET"\
             + "\t" + "DS_SPLIT_KEY"\
             + "\t" + "INCHI"\
             + "\t" + "INCHI_KEY"\
             + "\t" + "ACCESSION_NUMBER"\
             + "\t" + "CMPD_NAME" \
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
                + "\t" + node_info.spec_cluster.represen_spec_uni.accession_number \
                + "\t" + node_info.spec_cluster.represen_spec_uni.name \
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
