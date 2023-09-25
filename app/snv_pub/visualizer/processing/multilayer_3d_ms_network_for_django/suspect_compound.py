import pandas as pd
import sys
from . import adduct_calc
from . import formula

from logging import getLogger
logger = getLogger(__name__)

# MODIFIED20220714----------------
def read_file_suspect_mz( dic_conf ):
    dic_filename_vs_dic_name_vs_suspect_mz = {}
    for filename_suspect_mz in dic_conf["list_filename_suspect_mz"]:
        logger.debug("filename_suspect_mz", filename_suspect_mz)

        #filename_suspect_mz = dic_conf.filename_suspect_mz
        try:
            df_list_suspect_mz = pd.read_csv(filename_suspect_mz, delimiter='\t', index_col=False)
        except ValueError as e:
            logger.error(e)
            logger.error(" read_list_suspect_mz did not work. check  filename_suspect_mz\n  \
               check if file contain extra data out side of headear specified region \n \
                or presence of extra row")
            dic_conf["all_data_valid"] = False
            dic_conf["l_func_invalid_data"].append("read_file_suspect_mz")

        df_list_suspect_mz = df_list_suspect_mz.fillna("nd")

        colname_for_mz = ""
        colname_for_name = ""
        colname_for_formula = ""
        for colname in df_list_suspect_mz.columns:
            if colname in ["mz", "MZ", "m/z"]:
                colname_for_mz = colname
            if colname in ["Name","name", "NAME", "Chemical Name", "description", "annotation"]:
                colname_for_name = colname
            if colname in ["Molecular Formula", "molecular formula", "Formula", "formula", "Elemental Composition", "elemental composition", ]:
                colname_for_formula = colname

        # terminate if no column name for compound is not find.
        if colname_for_name == "":
            logger.error ("ERROR no column name for comound name was found in suspect file")
            dic_conf["all_data_valid"] = False
            dic_conf["l_func_invalid_data"].append("read_file_suspect_mz")

        # terminate if colum name for BOTH formula and mz were not foundd
        if colname_for_mz == ""  and colname_for_formula == ""  :
            logger.error ("ERROR no column name for comound name was found in suspect file")
            dic_conf["all_data_valid"] = False
            dic_conf["l_func_invalid_data"].append("read_file_suspect_mz")


        dic_name_vs_suspect_mz = {}

        # if table contain formula info -----------------------------
        if len(colname_for_formula) > 0:
            for key, row in df_list_suspect_mz.iterrows():
                # try to read name colum only when it exists
                if len(colname_for_name) > 0:
                    exact_mass = 0
                    suspect_formula = row[colname_for_formula]
                    suspect_name = row[colname_for_name]
                    try :
                        exact_mass = formula.get_monoisotopic_mass(suspect_formula)
                    except :
                        logger.error(" Error at read_file_suspect_mz function. something wrong with mass calculation from formula")
                        dic_conf["all_data_valid"] = False
                        dic_conf["l_func_invalid_data"].append("read_file_suspect_mz")


                    # add to suspect list if the mass is more than 50u
                    if exact_mass > 50.0 :
                        for adduct_mass in dic_conf["l_adduct_mass_for_suspect"] :
                            dic_name_vs_suspect_mz[suspect_name + "_adduct_" + str(adduct_mass)] = exact_mass + adduct_mass

                        for adduct_type in dic_conf["l_adduct_type_for_suspect"]:
                            mz = adduct_calc.get_mz_based_on_exactmw_and_adduct(exact_mass, adduct_type)
                            dic_name_vs_suspect_mz[suspect_name + "_adduct_" + adduct_type] = mz

        # if table has mz info,  but  NOT formula info -----------------------------
        if len(colname_for_mz) > 0  and len(colname_for_formula) == 0:
            for key, row in df_list_suspect_mz.iterrows():
                suspect_mz = float(row[colname_for_mz])

                # try to read name colum only when it exists
                if len(colname_for_name) > 0:
                    suspect_name = row[colname_for_name]
                    dic_name_vs_suspect_mz[suspect_name] = suspect_mz

        dic_filename_vs_dic_name_vs_suspect_mz[filename_suspect_mz] = dic_name_vs_suspect_mz

    return dic_filename_vs_dic_name_vs_suspect_mz

# ----------------MODIFIED20220714




