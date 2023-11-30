from logging import basicConfig, getLogger, DEBUG

logger = getLogger(__name__)
logger.setLevel(DEBUG)
basicConfig(
    format='[%(asctime)s] %(name)s %(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def get_initialized_dic_spec_cluster():
    """
    Returns
    -------
    dic_spec_cluster : dict
        {
            'cluster_id': '',
            'compound_name': '',
            'global_accession': '',
            'represen_dic_spec_uni': {},
            'list_dic_spec_uni': [],
            'list_COMMON_NAME': [],
            'list_external_compound_UNIQUE_ID': [],
            'list_pathway_UNIQUE_ID': [],
            'list_pathway_COMMON_NAME': [],
            'list_compound_categories': [],
            'inchi': '',
            'inchi_key': '',
            'tag': '',
            'dataset': '',
            'layer': '',
            "unique_level": -1,
            'is_toxic': False,
            'extra_node_color': ''
        }
    """
    dic_spec_cluster = {
        'cluster_id': '',
        'compound_name': '',
        'global_accession': '',
        'represen_dic_spec_uni': {},
        'list_dic_spec_uni': [],
        'list_COMMON_NAME': [],
        'list_external_compound_UNIQUE_ID': [],
        'list_pathway_UNIQUE_ID': [],
        'list_pathway_COMMON_NAME': [],
        'list_compound_categories': [],
        'inchi': '',
        'inchi_key': '',
        'tag': '',
        'dataset': '',
        'layer': '',
        "unique_level": -1,
        'is_toxic': False,
        'extra_node_color': ''
    }

    return dic_spec_cluster


def get_initialized_dic_node_info():
    """
    Returns
    -------
    dic_node_info : dict
        {
            'total_input_idx': 'none',
            'total_input_idx_mod': 'none',
            'dic_spec_cluster': {},
            'layer': '',
            'color': 0,
            'suspect': '',
            'l_source_suspect': [],
            'name': ''
        }
    """

    dic_node_info = {
        'total_input_idx': 'none',
        'total_input_idx_mod': 'none',
        'dic_spec_cluster': {},
        'layer': '',
        'color': 0,
        'suspect': '',
        'l_source_suspect': [],
        'name': ''
    }

    return dic_node_info


def get_initialized_dic_edge_info():
    """
    Returns
    -------
    dic_edge_info : dict
        {
            'spec_cluster_x_total_input_idx': 0,
            'spec_cluster_y_total_input_idx': 0,
            'spec_cluster_x_id': 0,
            'spec_cluster_y_id': 0,
            'spec_cluster_x_global_accession': '',
            'spec_cluster_y_global_accession': '',
            'delta_mz': -1,
            'spec_sim_score': -1,
            'number_matched_peaks': -1,
            'edge_type': '',
            'cmpd_sim_score': -1
        }
    """

    dic_edge_info = {
        'spec_cluster_x_total_input_idx': 0,
        'spec_cluster_y_total_input_idx': 0,
        'spec_cluster_x_id': 0,
        'spec_cluster_y_id': 0,
        'spec_cluster_x_global_accession': '',
        'spec_cluster_y_global_accession': '',
        'delta_mz': -1,
        'spec_sim_score': -1,
        'number_matched_peaks': -1,
        'edge_type': '',
        'cmpd_sim_score': -1
    }

    return dic_edge_info


def get_initialized_dic_spec_uni():
    """
    Returns
    -------
    dic_spec_uni : dict
        {
            'source': '',
            'accession_number': '',
            'unique_number': -1,
            'authors': '',
            'name': '',
            'cas_no': '',
            'inchi_key': '',
            'smiles': '',
            'inchi': '',
            'pubchem_cid': '',
            'KEGG_id': '',
            'ms_level': '',
            'spectrum_type': '',
            'ionization_mode': '',
            'instrument': '',
            'instrument_type': '',
            'fragmentation_type': '',
            'fragmentation_energy': 0,
            'mass_resolution_class': 0,
            'exact_mass': -1,
            'precursor_mz': -1,
            'charge_type': {},
            'precursor_type': {},
            'retention_time_in_sec': -1,
            'retention_time_in_min': -1,
            'peak_list_mz_int_abs': [],
            'peak_list_mz_int_rel': [],
            'peak_list_mz_annoid': [],
            'dic_annoid_struct': {},
            'list_cmpd_classification_kingdom': [],
            'list_cmpd_classification_class': [],
            'list_cmpd_classification_superclass': [],
            'list_cmpd_classification_alternative_parent': [],
            'source_filename': '',
            'source_filepath': ''
        }
    """
    dic_spec_uni = {
        'source': '',
        'accession_number': '',
        'unique_number': -1,
        'authors': '',
        'name': '',
        'cas_no': '',
        'inchi_key': '',
        'smiles': '',
        'inchi': '',
        'pubchem_cid': '',
        'KEGG_id': '',
        'ms_level': '',
        'spectrum_type': '',
        'ionization_mode': '',
        'instrument': '',
        'instrument_type': '',
        'fragmentation_type': '',
        'fragmentation_energy': 0,
        'mass_resolution_class': 0,
        'exact_mass': -1,
        'precursor_mz': -1,
        'charge_type': {},
        'precursor_type': {},
        'retention_time_in_sec': -1,
        'retention_time_in_min': -1,
        'peak_list_mz_int_abs': [],
        'peak_list_mz_int_rel': [],
        'peak_list_mz_annoid': [],
        'dic_annoid_struct': {},
        'list_cmpd_classification_kingdom': [],
        'list_cmpd_classification_class': [],
        'list_cmpd_classification_superclass': [],
        'list_cmpd_classification_alternative_parent': [],
        'source_filename': '',
        'source_filepath': ''
    }

    return dic_spec_uni

