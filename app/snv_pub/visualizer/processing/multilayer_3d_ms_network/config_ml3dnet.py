import configparser
import ast
import sys

from logging import basicConfig, getLogger, DEBUG

logger = getLogger(__name__)
logger.setLevel(DEBUG)
basicConfig(
    format='[%(asctime)s] %(name)s %(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def get_dic_config_initialized(setting=0):
    dic_config = {}

    # status of input/processed data.
    # [if all values are valid]
    dic_config["all_data_valid"] = True
    dic_config["l_func_invalid_data"] = []

    # [filename] ----------------------------------------------------------------------
    dic_config["filename_edge_info"] = "output_test_0\output_.edgeinfo.tsv"
    dic_config["filename_cluster_info"] = "output_test_0\output_.cluster_attribute.tsv"
    dic_config["filename_feature_table"] = ""
    dic_config["foldername_ext_cmpd_info"] = ""

    # [filter]------------------------
    # e.g.     none, list_cmpd_classification_superclass
    dic_config["filter_select_category"] = "none"
    # e.g.   ###   Phenylpropanoids and polyketides    Lipids and lipid-like molecules
    dic_config["filter_select_keyword"] = ""

    # [layer]--------------------------------
    #  you can choose from  : source_filename   , list_cmpd_classification_superclass , list_compound_categories
    #  if you want to make split layers according to superclass classification  :  list_cmpd_classification_superclass
    #  for categories of compound specified in external compound info (foldername_ext_cmpd_info)   :   list_compound_categories
    dic_config["type_attribute_for_layer_separation"] = "list_cmpd_classification_superclass"
    ## here you can define attribute (as string) to define which layer is base layer.
    ##  if you want to spec in file "base.mgf", then = base.mgf
    dic_config["str_key_attribute_to_base_layer"] = "sample"

    dic_config['color_toxic_compound'] = False

    #   0 : no mesh   1: with mesh
    dic_config["mesh_mode"] = 1

    # [mass_range] -------------
    dic_config["mass_lower_limit"] = 100.0
    dic_config["mass_higher_limit"] = 1000.0

    dic_config["dic_suspect_cmpd_vs_mz"] = {}

    # adduct mass for suspect list -----------------------------

    dic_config["list_filename_suspect_mz"] = ["PFAS_37_reg.tsv"]

    # list of mz for adduct calc.   e.g. [1.0072]
    # normally you should use l_adduct_type_for_suspect  below
    dic_config["l_adduct_mass_for_suspect"] = []

    # adduct type for suspect
    # currently not on ini file.
    dic_config["l_adduct_type_for_suspect"] = ["M+H"]

    # mass defect--------------------------------------------------
    dic_config["list_mass_defect"] = []

    # product/fragment ion required
    dic_config["mz_tolerance_for_fragment"] = 0.01
    dic_config["list_product_mz_required"] = []

    # mz tolerance -------
    dic_config["mz_tol"] = 0.01

    # [threshold]--------------
    dic_config["score_threshold"] = 0.75

    # network realated ----------------
    dic_config["subgraph_type"] = "node_quant"
    dic_config["subgraph_num_core_nodes"] = 1000
    dic_config["subgraph_depth"] = 3
    dic_config["stat_value"] = 0.05
    dic_config["quant_polar"] = "qauant_val_both"
    # if 0, quant value is ignored   if 2, ratio  < 1/2 or  > 2 will be taken for subgraph
    dic_config["quant_value"] = 2

    # 'ratio' or 'loading' are available for 'quant_value_type'
    dic_config['quant_value_type'] = 'ratio'

    dic_config["l_global_accession_for_node_select_subgraph"] = []

    dic_config["l_total_input_idx_to_remove"] = []

    dic_config["node_select_subgraph_depth"] = 10  # TODO: Add (K. Hirata, 20221214)

    # [community] --------------------
    # 0: no community detection, 1- : level (round) of community detection
    dic_config["n_level_community_detection"] = 0

    return dic_config


def read_config_file(config_filename):
    # read config

    # config_o = Config_multilayer_3d()

    dic_config = {}
    #####################
    # [B]read config file ----------------------------------------------------------------
    ####################
    inifile = configparser.SafeConfigParser()
    inifile.read("./config.ini")

    # file name ----------------------------------------------------------------------
    dic_config["filename_edge_info"] = inifile.get("filename", "filename_edge_info")
    dic_config["filename_cluster_info"] = inifile.get("filename", "filename_cluster_info")
    dic_config["filename_feature_table"] = inifile.get("filename", "feature_table")
    dic_config["foldername_ext_cmpd_info"] = inifile.get("filename", "foldername_ext_cmpd_info")

    dic_config["foldername_spectra"] = inifile.get("files", "spectra_folder")

    dic_config["type_attribute_for_layer_separation"] = "source_filename"

    dic_config["filter_select_category"] = inifile.get("filter", "select_category")
    dic_config["filter_select_keyword"] = inifile.get("filter", "select_keyword")

    dic_config["type_attribute_for_layer_separation"] = inifile.get("layer", "type_attribute_for_layer_separation")

    dic_config["str_key_attribute_to_base_layer"] = '[]'

    str_from_file = inifile.get("layer", "str_key_attribute_to_base_layer")
    if len(str_from_file) > 1:
        dic_config["str_key_attribute_to_base_layer"] = str_from_file

    dic_config["mesh_mode"] = int(inifile.get("layer", "mesh_mode"))

    # mass range -------------
    dic_config["mass_lower_limit"] = float(inifile.get("mass_range", "mass_lower_limit"))
    dic_config["mass_higher_limit"] = float(inifile.get("mass_range", "mass_higher_limit"))

    # MODIFIED20220714-------------
    #  dictionary suspect compound vs mz ----------------------------
    line_v = inifile.get("mass_range", "list_suspect_mz")
    # dic_as_str = (line_v.strip().rstrip('\n').split(","))
    try:
        dic_config["dic_suspect_cmpd_vs_mz"] = ast.literal_eval(line_v)

    except ValueError:

        print(" something wrong with -mass_range-dic_suspect_cmpd_vs_mz-")
        # print (dic_as_str)
        dic_config["all_data_valid"] = False
        dic_config["l_func_invalid_data"].append("read_config_file")

    # altearnatively, accept list of file supect---------------------------------------------
    line_v = inifile.get("mass_range", "list_filename_suspect_mz")
    list_as_str = (line_v.strip().rstrip('\n').split(","))
    list_as_str = [e.strip() for e in list_as_str]
    print(list_as_str)
    list_filename_suspect = []
    for e in list_as_str:
        # avoid adding " " as filename

        list_filename_suspect.append(e)
        print("appended")
    dic_config["list_filename_suspect_mz"] = list_filename_suspect

    # adduct mass for suspect list -----------------------------
    line_v = inifile.get("mass_range", "list_adduct_mass_for_suspect")

    # only when mass defect parameter is specified.  set to -1 or " " if mass defect is not used
    if not line_v == "none":
        list_as_str = (line_v.strip().rstrip('\n').split(","))

        print(list_as_str)
        if len(list_as_str) > 0:
            try:
                dic_config["l_adduct_mass_for_suspect"] = [float(val) for val in list_as_str]
            except:
                print(" something wrong with -mass_range-list_adduct_mass_for_suspect-")
                dic_config["all_data_valid"] = False
                dic_config["l_func_invalid_data"].append("read_config_file")

    # mass defect--------------------------------------------------
    line_v = inifile.get("mass_range", "list_mass_defect")

    dic_config["list_mass_defect"] = []

    # only when mass defect parameter is specified.  set to -1 or " " if mass defect is not used
    if not line_v == "none":
        list_as_str = (line_v.strip().rstrip('\n').split(","))

        print(list_as_str)
        if len(list_as_str) > 0:
            try:
                dic_config["list_mass_defect"] = [float(val) for val in list_as_str]
            except ValueError:
                print(" something wrong with -mass_range-list_mass_defect-")
                dic_config["all_data_valid"] = False
                dic_config["l_func_invalid_data"].append("read_config_file")

    # product/fragment ion required
    dic_config["list_product_mz_required"] = []

    line_v = inifile.get("mass_range", "list_product_mz_required")

    if not line_v == "none":
        list_as_str = (line_v.strip().rstrip('\n').split(","))
        try:
            dic_config["list_product_mz_required"] = [float(val) for val in list_as_str]
        except ValueError:
            print(" something wrong with -mass_range-list_product_mz_required-")
            dic_config["all_data_valid"] = False
            dic_config["l_func_invalid_data"].append("read_config_file")

    # mz tolerance -------
    dic_config["mz_tol"] = float(inifile.get("mass_range", "mz_tol"))

    # ---------MODIFIED20220714

    # network realated ----------------
    dic_config["subgraph_type"] = inifile.get("subgraph", "subgraph_type")
    dic_config["subgraph_num_core_nodes"] = int(inifile.get("subgraph", "num_core_nodes"))
    dic_config["subgraph_depth"] = int(inifile.get("subgraph", "depth"))
    dic_config["stat_value"] = float(inifile.get("subgraph", "stat_value"))
    dic_config["quant_ploar"] = "qauant_val_both"
    dic_config["quant_value"] = float(inifile.get("subgraph", "quant_value"))

    dic_config["n_level_community_detection"] = int(inifile.get("community", "level_community_detection"))

    dic_config["l_global_accession_for_node_select_subgraph"] = []

    line_v = inifile.get("subgraph", "l_total_input_idx_to_remove")
    if not line_v == "none":
        list_as_str = (line_v.strip().rstrip('\n').split(","))
        try:
            dic_config["l_total_input_idx_to_remove"] = [str(val) for val in list_as_str]
        except ValueError:
            print(" something wrong with -subgraph-l_total_input_idx_to_remove-")
            dic_config["all_data_valid"] = False
            dic_config["l_func_invalid_data"].append("read_config_file")

    # threshold ----
    dic_config["score_threshold"] = float(inifile.get("threshold", "score_threshold"))

    #####################
    #    output
    #################

    # output file names
    dic_config["output_folder_name"] = inifile.get("output", "output_folder_name")
    dic_config["output_file_name"] = inifile.get("output", "output_file_name")

    line_v = inifile.get("output", "list_id_for_graphical_scoring_details")
    list_as_str = (line_v.strip().rstrip('\n').split(","))
    list_as_str = [val.rstrip().lstrip() for val in list_as_str]
    dic_config["list_id_for_graphical_scoring_details"] = [int(val) for val in list_as_str]

    fo_log = open("log.txt", "w")
    str_o = ""
    str_o = "score threshold is " + str(dic_config["score_threshold"])
    fo_log.write(str_o)
    str_o = ""

    return dic_config
