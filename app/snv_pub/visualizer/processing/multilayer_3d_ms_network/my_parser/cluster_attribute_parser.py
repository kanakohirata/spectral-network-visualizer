import ast
import dask.dataframe as dd
from logging import getLogger
import pandas as pd
import re
from .. import dedicated_dictionaries


logger = getLogger(__name__)


def read_cluster_attribute_old(dic_config):
    ###########################################################################################
    # section D
    # [D] read spec cluster info
    ##########################################################################################

    l_cluster_total_input_idx = []

    # this dictionary is inted to hold whole input cluster info.
    # cluster_total_input_idx is simply count of input cluster. NOT cluster id created by generator

    dic_cluster_total_input_idx_vs_cluster_info = {}

    """:type: list[Spec_cluster]"""
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  use conf_o
    # df_clusterinfo = pd.read_csv(conf_o.filename_cluster_info, delimiter='\t', index_col=False)

    df_clusterinfo = dd.read_csv(dic_config["filename_cluster_info"], sep='\t', dtype='object')

    # first scan clusterinfo file just to collect cluster number.
    count = -1
    str_o = ""

    f_present_RETENTION_TIME_IN_SEC = 0
    for col in df_clusterinfo.columns:

        logger.debug(col)
        if col == "RETENTION_TIME_IN_SEC":
            f_present_RETENTION_TIME_IN_SEC = 1
            logger.debug("rt present")

    for key, row in df_clusterinfo.iterrows():
        valid_entry = 1
        count += 1

        l_cluster_total_input_idx.append(str(count))
        ############################################
        # make new object------------------------------
        dic_clst = dedicated_dictionaries.get_initialized_dic_spec_cluster()

        ############################################
        # get info from the file-----------------
        dic_clst["cluster_id"] = row['CLUSTER_ID']

        # !!!!!!!!!!!!!!!!!!!!!!!!!!
        # get rid of "  mark, which cause trouble later on in Dash interface.
        str_global_accession = row['GLOBAL_ACCESSION']
        dic_clst["compound_name"] = row['CMPD_NAME']
        dic_clst["global_accession"] = str_global_accession.replace('"', '')

        dic_clst["dataset"] = row['DATASET']
        dic_clst["tag"] = row['DATASET']

        dic_clst["list_compound_categories"] = []

        if row['INCHI'] is not None:
            dic_clst["inchi"] = row['INCHI']

        if row['INCHI_KEY'] is not None:
            dic_clst["inchi_key"] = str(row['INCHI_KEY'])

        dic_spec_uni = dedicated_dictionaries.get_initialized_dic_spec_uni()
        dic_spec_uni["precursor_mz"] = float((row['PRECURSOR_MZ']))

        if f_present_RETENTION_TIME_IN_SEC == 1:
            if (row['RETENTION_TIME_IN_SEC'] != "None") and (row['RETENTION_TIME_IN_SEC'] is not "None"):
                dic_spec_uni["retention_time_in_sec"] = float((row['RETENTION_TIME_IN_SEC']))

        dic_spec_uni["accession_number"] = (row['ACCESSION_NUMBER'])

        # peak list--------------
        dic_spec_uni["peak_list_mz_int_rel"] = list(ast.literal_eval(row['PEAK_LIST_MZ_INT_REL']))

        dic_spec_uni["list_cmpd_classification_kingdom"] = list(ast.literal_eval(row['CMPD_CLASSIFICATION_KINGDOM']))
        dic_spec_uni["list_cmpd_classification_superclass"] = list(
            ast.literal_eval(row['CMPD_CLASSIFICATION_SUPERCLASS']))
        dic_spec_uni["list_cmpd_classification_class"] = list(ast.literal_eval(row['CMPD_CLASSIFICATION_CLASS']))
        dic_spec_uni["source_filename"] = row['FILE_NAME']
        dic_clst["represen_spec_uni"] = dic_spec_uni
        # !!!!!!!!!not using ???
        dic_clst["list_spec_uni"] = []
        
        if valid_entry:
            # note idx is now converted to str, because it will be modified for layer info later,
            dic_cluster_total_input_idx_vs_cluster_info[str(count)] = dic_clst

        str_o = str_o + "\n" + "cluster_total_input_idx:" + str(count) + "\tglobal accession:" + str(
            dic_clst["global_accession"]) + "\tprecmz:" + str(dic_spec_uni["precursor_mz"]) + \
                dic_spec_uni["source_filename"] + str(dic_spec_uni["list_cmpd_classification_kingdom"])

    return dic_cluster_total_input_idx_vs_cluster_info


def read_cluster_attribute(path):
    """
    Read cluster attribute file.
    Parameters
    ----------
    path: str
        Path of cluster attribute file.

    Returns
    -------
    A dictionary of cluster attribute.
        {'0': {'cluster_id': 0, 'global_accession': 'aaa', 'dataset': 'sample', ...},
         '1': {'cluster_id': 1, 'global_accession': 'bbb', 'dataset': 'sample', ...}, ...}
    """

    # this dictionary is intended to hold whole input cluster info.
    # cluster_total_input_idx is simply count of input cluster. NOT cluster id created by generator

    df_cluster_info = pd.read_csv(path, sep='\t')
    flag_present_retention_time_in_sec = False
    for col in df_cluster_info.columns:

        logger.debug(col)
        if col == "RETENTION_TIME_IN_SEC":
            flag_present_retention_time_in_sec = True
            logger.debug('"RETENTION_TIME_IN_SEC" is in cluster attribute columns.')

    # Rename columns.
    df_cluster_info.rename(inplace=True,
                           columns={'CLUSTER_ID': 'cluster_id',
                                    'GLOBAL_ACCESSION': 'global_accession',
                                    'CMPD_NAME': 'compound_name',
                                    'DATASET': 'dataset',
                                    'INCHI': 'inchi',
                                    'INCHI_KEY': 'inchi_key',
                                    'PRECURSOR_MZ': 'precursor_mz',
                                    'ACCESSION_NUMBER': 'accession_number',
                                    'PEAK_LIST_MZ_INT_REL': 'peak_list_mz_int_rel',
                                    'CMPD_CLASSIFICATION_KINGDOM': 'list_cmpd_classification_kingdom',
                                    'CMPD_CLASSIFICATION_SUPERCLASS': 'list_cmpd_classification_superclass',
                                    'CMPD_CLASSIFICATION_CLASS': 'list_cmpd_classification_class',
                                    'FILE_NAME': 'source_filename'})

    # Correct retention time data.
    if flag_present_retention_time_in_sec:
        df_cluster_info.rename(inplace=True, columns={'RETENTION_TIME_IN_SEC': 'retention_time_in_sec'})
        df_cluster_info['retention_time_in_sec'] = df_cluster_info['retention_time_in_sec'].astype(str)
        df_cluster_info['retention_time_in_sec'] =\
            df_cluster_info['retention_time_in_sec'].apply(lambda x: 0.0 if not re.match(r'\d+\.?\d+', x) else x)
        df_cluster_info['retention_time_in_sec'] = df_cluster_info['retention_time_in_sec'].fillna(0.0)
    else:
        df_cluster_info['retention_time_in_sec'] = 0.0

    df_cluster_info['cluster_id'] = df_cluster_info['cluster_id'].astype(str)

    # Fill compound data.
    df_cluster_info[['compound_name', 'inchi', 'inchi_key']]\
        = df_cluster_info[['compound_name', 'inchi', 'inchi_key']].fillna('')

    # Add a column for first part of inchikey.
    df_cluster_info['inchi_key_first'] = df_cluster_info['inchi_key'].apply(lambda x: x.split('-')[0])

    # Fill precursor m/z and retention time.
    df_cluster_info['precursor_mz'] = df_cluster_info['precursor_mz'].astype(str)
    df_cluster_info['precursor_mz'] =\
        df_cluster_info['precursor_mz'].apply(lambda x: 0.0 if not re.match(r'\d+\.?\d+', x) else x)
    df_cluster_info['precursor_mz'] = df_cluster_info['precursor_mz'].fillna(0.0)
    df_cluster_info[['precursor_mz', 'retention_time_in_sec']]\
        = df_cluster_info[['precursor_mz', 'retention_time_in_sec']].astype({'precursor_mz': float, 'retention_time_in_sec': float})

    # Initialize other columns.
    df_cluster_info['tag'] = df_cluster_info['dataset']
    df_cluster_info['list_compound_categories'] = [[], ] * len(df_cluster_info)
    df_cluster_info['list_spec_uni'] = [[], ] * len(df_cluster_info)
    df_cluster_info['list_COMMON_NAME'] = [[], ] * len(df_cluster_info)
    df_cluster_info['list_external_compound_UNIQUE_ID'] = [[], ] * len(df_cluster_info)
    df_cluster_info['list_pathway_UNIQUE_ID'] = [[], ] * len(df_cluster_info)
    df_cluster_info['list_pathway_COMMON_NAME'] = [[], ] * len(df_cluster_info)
    df_cluster_info[['layer', 'extra_node_color']] = ''
    df_cluster_info['unique_level'] = -1
    df_cluster_info['is_toxic'] = False

    # Split df_cluster_info into new df_cluster_info and df_spectrum
    df_spectrum = df_cluster_info[['accession_number',
                                   'source_filename',
                                   'precursor_mz',
                                   'retention_time_in_sec',
                                   'peak_list_mz_int_rel',
                                   'list_cmpd_classification_kingdom',
                                   'list_cmpd_classification_superclass',
                                   'list_cmpd_classification_class']]

    df_cluster_info = df_cluster_info[['cluster_id',
                                       'global_accession',
                                       'compound_name',
                                       'dataset',
                                       'tag',
                                       'list_compound_categories',
                                       'inchi',
                                       'inchi_key',
                                       'inchi_key_first',
                                       'list_spec_uni',
                                       'list_COMMON_NAME',
                                       'list_external_compound_UNIQUE_ID',
                                       'list_pathway_UNIQUE_ID',
                                       'list_pathway_COMMON_NAME',
                                       'layer',
                                       'extra_node_color',
                                       'unique_level',
                                       'is_toxic']]

    # Initialize other columns of df_spectrum.
    df_spectrum[['authors',
                 'name',
                 'cas_no',
                 'inchi_key',
                 'inchi',
                 'smiles',
                 'pubchem_cid',
                 'KEGG_id',
                 'ms_level',
                 'spectrum_type',
                 'instrument',
                 'instrument_type',
                 'fragmentation_type',
                 'source_filepath']] = ''
    
    df_spectrum[['unique_number',
                 'exact_mass',
                 'retention_time_in_min']] = -1
    
    df_spectrum[['fragmentation_energy',
                 'mass_resolution_class']] = 0
    
    df_spectrum['charge_type'] = [{}, ] * len(df_spectrum)
    df_spectrum['precursor_type'] = [{}, ] * len(df_spectrum)
    df_spectrum['peak_list_mz_int_abs'] = [[], ] * len(df_spectrum)
    df_spectrum['peak_list_mz_annoid'] = [[], ] * len(df_spectrum)
    df_spectrum['dic_annoid_struct'] = [{}, ] * len(df_spectrum)
    df_spectrum['list_cmpd_classification_alternative_parent'] = [[], ] * len(df_spectrum)

    # Convert peaks from str to list
    df_spectrum['peak_list_mz_int_rel'] = df_spectrum['peak_list_mz_int_rel'].apply(_parse_string_spectrum)
    # Convert classification from str to list
    df_spectrum['list_cmpd_classification_kingdom'] = df_spectrum['list_cmpd_classification_kingdom'].apply(_parse_string_classification_list)
    df_spectrum['list_cmpd_classification_superclass'] = df_spectrum['list_cmpd_classification_superclass'].apply(_parse_string_classification_list)
    df_spectrum['list_cmpd_classification_class'] = df_spectrum['list_cmpd_classification_class'].apply(_parse_string_classification_list)

    spectra = df_spectrum.to_dict(orient='records')
    df_cluster_info['represen_spec_uni'] = spectra

    dic_cluster_total_input_idx_vs_cluster_info = df_cluster_info.to_dict(orient='index')

    return dic_cluster_total_input_idx_vs_cluster_info

    
def _parse_string_spectrum(string_spectrum):
    if pd.isna(string_spectrum):
        return []

    string_numbers = re.findall(r'\d+\.\d+', string_spectrum)
    numbers = [float(_x) for _x in string_numbers]
    spectrum = [list(numbers[i:i + 2]) for i in range(0, len(numbers), 2)]
    return spectrum


def _parse_string_classification_list(string_classification_list):
    if pd.isna(string_classification_list):
        return []

    if '"' in string_classification_list:
        classifications = re.findall(r'"(.*?)"', string_classification_list)
    else:
        classifications = re.findall(r"'(.*?)'", string_classification_list)
    return classifications

