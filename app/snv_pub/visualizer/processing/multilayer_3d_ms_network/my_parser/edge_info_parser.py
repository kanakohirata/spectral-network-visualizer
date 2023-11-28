import dask.dataframe as dd
from logging import getLogger
from .. import dedicated_dictionaries


logger = getLogger(__name__)


def read_edge_info_old(dic_config):
    ###########################################################################################
    # section [H]
    # [H] read EDGES info
    ##########################################################################################
    l_edges = []
    list_edge_info = []
    list_cluster_global_accession_to_pass_threshold = []
    """:type: list[Edge_info]"""

    # df_edge_info = pd.read_csv(conf_o.filename_edge_info, delimiter='\t', index_col=False)
    df_edge_info = dd.read_csv(dic_config["filename_edge_info"], sep='\t', dtype='object', encoding='utf-8')
    logger.debug("   --read_edge_info--  reading edgeinfo")

    for key, row in df_edge_info.iterrows():
        # obj_edge_info = edge_info_a1.Edge_info()
        # obj_edge_info = Edge_info()
        dic_edge_info = dedicated_dictionaries.get_initialized_dic_edge_info()
        dic_edge_info["spec_cluster_x_id"] = row['X_CLUSTERID']
        # obj_edge_info.spec_cluster_x_global_accession = row['REP_SPECTRUM_X_GLOBAL_ACCESSION']
        x_global_accession = row['REP_SPECTRUM_X_GLOBAL_ACCESSION']
        dic_edge_info["spec_cluster_x_global_accession"] = x_global_accession.replace('"', '')

        dic_edge_info["spec_cluster_y_id"] = row['Y_CLUSTERID']
        # obj_edge_info.spec_cluster_y_global_accession = row['REP_SPECTRUM_Y_GLOBAL_ACCESSION']
        y_global_accession = row['REP_SPECTRUM_Y_GLOBAL_ACCESSION']
        dic_edge_info["spec_cluster_y_global_accession"] = y_global_accession.replace('"', '')

        l_edges.append((row['X_CLUSTERID'], row['Y_CLUSTERID']))

        dic_edge_info["spec_sim_score"] = float(row['SCORE'])
        dic_edge_info["delta_mz"] = float(row['DELTA_MZ'])

        if dic_edge_info["spec_sim_score"] > dic_config["score_threshold"]:
            list_edge_info.append(dic_edge_info)
            # list_cluster_global_accession_to_pass_threshold.append( dic_edge_info["spec_cluster_x_global_accession"] )
            # list_cluster_global_accession_to_pass_threshold.append(dic_edge_info["spec_cluster_y_global_accession"])
            # make it nr
            # list_cluster_global_accession_to_pass_threshold = list(set(list_cluster_global_accession_to_pass_threshold))

    logger.debug(" -- read_edge_info-- FINISHED reading edgeinfo")

    # fo_log.write("\n[E] reading edge info\n")
    # fo_log.write( "\nlen of list_cluster_global_accession_to_pass_threshold:\n" + str(len(list_cluster_global_accession_to_pass_threshold))  + "\n")
    # fo_log.write("\nlen of list_edge_info:\n" + str(len(list_edge_info)) + "\n")

    return l_edges, list_edge_info


def read_edge_info(path, score_threshold):
    """
    Read edge info file.
    Parameters
    ----------
    path : str
        Path of edge info file.
    score_threshold : float
        Threshold of spectral similarity score.

    Returns
    -------

    """

    l_edges = []
    list_edge_info = []

    df_edge_info = dd.read_csv(path, sep='\t', dtype='object', encoding='utf-8')
    logger.debug('  --read_edge_info--  reading edge info')

    for key, row in df_edge_info.iterrows():
        dic_edge_info = dedicated_dictionaries.get_initialized_dic_edge_info()
        dic_edge_info['spec_cluster_x_id'] = row['X_CLUSTERID']
        dic_edge_info['spec_cluster_x_global_accession'] = row['REP_SPECTRUM_X_GLOBAL_ACCESSION']

        dic_edge_info['spec_cluster_y_id'] = row['Y_CLUSTERID']
        dic_edge_info['spec_cluster_y_global_accession'] = row['REP_SPECTRUM_Y_GLOBAL_ACCESSION']

        l_edges.append((row['X_CLUSTERID'], row['Y_CLUSTERID']))

        dic_edge_info['spec_sim_score'] = float(row['SCORE'])
        dic_edge_info['delta_mz'] = float(row['DELTA_MZ'])

        if dic_edge_info['spec_sim_score'] > score_threshold:
            list_edge_info.append(dic_edge_info)

    logger.debug("  --read_edge_info--  FINISHED reading edge info")

    return l_edges, list_edge_info
