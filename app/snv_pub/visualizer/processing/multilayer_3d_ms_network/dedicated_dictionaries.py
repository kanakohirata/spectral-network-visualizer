from logging import basicConfig, getLogger, DEBUG

logger = getLogger(__name__)
logger.setLevel(DEBUG)
basicConfig(
    format='[%(asctime)s] %(name)s %(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def get_initialized_dic_spec_cluster():

    dic_spec_cluster = {}
    dic_spec_cluster["cluster_id"] = ""
    dic_spec_cluster["compound_name"] = ""
    dic_spec_cluster["global_accession"] = ""
    dic_spec_cluster["represen_dic_spec_uni"] = {}
    dic_spec_cluster["list_dic_spec_uni"] =[]


    dic_spec_cluster["list_COMMON_NAME"] =[]
    dic_spec_cluster["list_external_compound_UNIQUE_ID"] =[]
    dic_spec_cluster["list_pathway_UNIQUE_ID"] = []
    dic_spec_cluster["list_pathway_COMMON_NAME"] = []


    dic_spec_cluster["list_compound_categories"] =[]
    dic_spec_cluster["inchi"] =""
    dic_spec_cluster["inchi_key"] =""
    dic_spec_cluster["tag"] =""
    dic_spec_cluster["dataset"] =""
    dic_spec_cluster["layer"] =""
    dic_spec_cluster["unique_level"] = -1

    dic_spec_cluster['is_toxic'] = False
    dic_spec_cluster['extra_node_color'] = ''

    return dic_spec_cluster





def get_initialized_dic_node_info():

    dic_node_info ={}

    dic_node_info["total_input_idx"] = "none"
    dic_node_info["total_input_idx_mod"] = "none"
    dic_node_info["dic_spec_cluster"] = {}
    dic_node_info["layer"] = ""
    dic_node_info["color"] = 0
    dic_node_info["suspect"] = ""
    dic_node_info["l_source_suspect"] = []
    dic_node_info["name"] = ""


    return dic_node_info


def get_initialized_dic_edge_info():
    dic_edge_info ={}

    dic_edge_info["spec_cluster_x_total_input_idx"] = 0
    dic_edge_info["spec_cluster_y_total_input_idx"] = 0
    dic_edge_info["spec_cluster_x_global_accession"] = ""
    dic_edge_info["spec_cluster_y_global_accession"] = ""
    dic_edge_info["delta_mz"] = -1
    dic_edge_info["spec_sim_score"] = -1
    dic_edge_info["number_matched_peaks"] = -1
    dic_edge_info["edge_type"] = ""
    dic_edge_info["cmpd_sim_score"] = -1

    return dic_edge_info


def get_initialized_dic_spec_uni():
    dic_spec_uni ={}

    dic_spec_uni["source"] = ""
    dic_spec_uni["accession_number"] = ""
    dic_spec_uni["unique_number"] = -1
    dic_spec_uni["authors"] = ""
    dic_spec_uni["name"] = ""
    dic_spec_uni["cas_no"] = ""
    dic_spec_uni["inchi_key"] = ""
    dic_spec_uni["smiles"] = ""
    dic_spec_uni["inchi"] = ""
    dic_spec_uni["pubchem_cid"] =     ""
    dic_spec_uni["KEGG_id"] = ""

    dic_spec_uni["ms_level"] = ""
    dic_spec_uni["spectrum_type"] = ""
    dic_spec_uni["ionization_mode"] = ""

    dic_spec_uni["instrument"] = ""
    dic_spec_uni["instrument_type"] = ""
    dic_spec_uni["fragmentation_type"] = ""
    dic_spec_uni["fragmentation_energy"] = 0
    dic_spec_uni["mass_resolution_class"] = 0



    dic_spec_uni["exact_mass"] = -1
    dic_spec_uni["precursor_mz"] = -1

    dic_spec_uni["charge_type"] = {}
    dic_spec_uni["precursor_type"] = {}


    dic_spec_uni["retention_time_in_sec"] = -1
    dic_spec_uni["retention_time_in_min"] = -1
    dic_spec_uni["peak_list_mz_int_abs"] = []
    dic_spec_uni["peak_list_mz_int_rel"] = []
    dic_spec_uni["peak_list_mz_annoid"] = []
    dic_spec_uni["dic_annoid_struct"] = {}

    dic_spec_uni["list_cmpd_classification_kingdom"] = []
    dic_spec_uni["list_cmpd_classification_class"] = []
    dic_spec_uni["list_cmpd_classification_superclass"] = []
    dic_spec_uni["list_cmpd_classification_alternative_parent"] = []


    dic_spec_uni["source_filename"] = ""
    dic_spec_uni["source_filepath"] = ""

    return dic_spec_uni

