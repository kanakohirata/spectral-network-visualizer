import os.path
import pandas as pd
from . import adduct_calc
from . import formula

from logging import getLogger

logger = getLogger(__name__)


def read_file_suspect_mz(dic_conf):
    dic_filename_vs_dic_name_vs_suspect_mz = {}
    for path_suspect_mz in dic_conf["list_filename_suspect_mz"]:
        filename_suspect_mz = os.path.basename(path_suspect_mz)
        logger.debug(f"filename_suspect_mz: {filename_suspect_mz}", )

        try:
            df_list_suspect_mz = pd.read_csv(path_suspect_mz, delimiter='\t', index_col=False)
        except ValueError as e:
            logger.error(e)
            logger.error(" read_list_suspect_mz did not work. check  filename_suspect_mz\n"
                         "check if file contain extra data out side of header specified region \n"
                         "or presence of extra row")
            dic_conf["all_data_valid"] = False
            dic_conf["l_func_invalid_data"].append("read_file_suspect_mz")
            return {}

        # Rename columns: replace ' ' with '_' and lowercase.
        columns_to_rename = {}
        for _column in df_list_suspect_mz.columns:
            columns_to_rename[_column] = _column.replace(' ', '_').lower()
        df_list_suspect_mz.rename(inplace=True, columns=columns_to_rename)

        colname_for_mz = ""
        colname_for_name = ""
        colname_for_formula = ""
        for colname in df_list_suspect_mz.columns:
            if colname in ["mz", "m/z"]:
                colname_for_mz = colname
            if colname in ["name", "chemical_name", "description", "annotation", "compound_name"]:
                colname_for_name = colname
            if colname in ["molecular_formula", "formula", "elemental_composition"]:
                colname_for_formula = colname

        # terminate if no column name for compound is not find.
        if colname_for_name == "":
            logger.error('ERROR: no column name "name" was found in suspect file')
            dic_conf["all_data_valid"] = False
            dic_conf["l_func_invalid_data"].append("read_file_suspect_mz")
            return {}

        # terminate if colum name for BOTH formula and mz were not founded
        if colname_for_mz == "" and colname_for_formula == "":
            logger.error('ERROR: "formula" and "m/z" were found in suspect file columns')
            dic_conf["all_data_valid"] = False
            dic_conf["l_func_invalid_data"].append("read_file_suspect_mz")
            return {}

        df_list_suspect_mz = df_list_suspect_mz[~pd.isna(df_list_suspect_mz[colname_for_name])]
        df_list_suspect_mz.drop_duplicates(subset=[colname_for_name], inplace=True)
        dic_name_vs_suspect_mz = {}

        # if table contain formula info -----------------------------
        if colname_for_formula:
            if not (dic_conf["l_adduct_mass_for_suspect"] or dic_conf["l_adduct_type_for_suspect"]):
                dic_name_vs_suspect_mz = {}
            else:
                df_list_suspect_mz = df_list_suspect_mz[~pd.isna(df_list_suspect_mz[colname_for_formula])]

                df_formula = pd.DataFrame(df_list_suspect_mz[colname_for_formula])
                # Remove duplicated formula
                df_formula.drop_duplicates(inplace=True)
                try:
                    df_formula["exact_mass"] = df_formula[colname_for_formula].apply(formula.get_monoisotopic_mass)
                except:
                    logger.error(
                        "Error at read_file_suspect_mz function. something wrong with mass calculation from formula")
                    dic_conf["all_data_valid"] = False
                    dic_conf["l_func_invalid_data"].append("read_file_suspect_mz")
                    return {}

                # Add to suspect list if the mass is more than 50u
                df_formula = df_formula[df_formula["exact_mass"] > 50]
                for _adduct_mass in dic_conf["l_adduct_mass_for_suspect"]:
                    df_formula[str(_adduct_mass)] = df_formula["exact_mass"].apply(lambda x: x + _adduct_mass)

                for _adduct_type in dic_conf["l_adduct_type_for_suspect"]:
                    df_formula[_adduct_type] = df_formula["exact_mass"].apply(
                        lambda x: adduct_calc.get_mz_based_on_exactmw_and_adduct(x, _adduct_type))

                df_list_suspect_mz = pd.merge(df_list_suspect_mz, df_formula, on=colname_for_formula, how='left')
                for adduct_mass in dic_conf["l_adduct_mass_for_suspect"]:
                    indexes = df_list_suspect_mz[colname_for_name] + f"_adduct_{adduct_mass}"
                    ser_suspect = pd.Series(df_list_suspect_mz[str(adduct_mass)].to_list(), index=indexes)
                    dic_name_vs_suspect_mz = {**dic_name_vs_suspect_mz, **ser_suspect.to_dict()}

                for adduct_type in dic_conf["l_adduct_type_for_suspect"]:
                    indexes = df_list_suspect_mz[colname_for_name] + f"_adduct_{adduct_type}"
                    ser_suspect = pd.Series(df_list_suspect_mz[adduct_type].to_list(), index=indexes)
                    dic_name_vs_suspect_mz = {**dic_name_vs_suspect_mz, **ser_suspect.to_dict()}

        # if table has mz info, but NOT formula info -----------------------------
        elif colname_for_mz:
            df_list_suspect_mz = df_list_suspect_mz[~pd.isna(df_list_suspect_mz[colname_for_mz])]
            df_list_suspect_mz[colname_for_mz] = df_list_suspect_mz[colname_for_mz].astype(float)
            indexes = df_list_suspect_mz[colname_for_name] + "_m/z_" + df_list_suspect_mz[colname_for_mz].astype(str)
            ser_suspect = pd.Series(df_list_suspect_mz[colname_for_mz].to_list(), index=indexes)
            dic_name_vs_suspect_mz = ser_suspect.to_dict()

        dic_filename_vs_dic_name_vs_suspect_mz[filename_suspect_mz] = dic_name_vs_suspect_mz

    return dic_filename_vs_dic_name_vs_suspect_mz
