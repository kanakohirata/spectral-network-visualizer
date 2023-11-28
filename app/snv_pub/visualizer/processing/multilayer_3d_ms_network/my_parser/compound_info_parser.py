from glob import glob
from django.core.cache import cache
from logging import getLogger
import os
import pandas as pd
import sys
import zipfile
from .. import read_t3db_a1

logger = getLogger(__name__)


def read_external_cmpd_info_old(dic_config, dic_cluster_total_input_idx_vs_cluster_info):
    logger.info(f'Start {sys._getframe().f_code.co_name}()')
    fo = open("read_external_cmpd_info.txt", "w")
    ###########################################
    #  [F] Read external info for node/cluster
    #############################################

    #####################
    # inchi key match mode    0: exact match ,  1:  XXXX-YYY-ZZZ match to XXXX-AAA-BBB
    inchi_key_match_mode = 1

    #########################
    # external compound CLASSYFIRE classification file.
    #############################

    filename_ext_classyfire_compound_table = "visualizer/data/cmpd_classification_classyfire.tsv"

    dic_inchikey_vs_dic_classyfire_classification = {}
    dic_inchikey_1st_part_vs_dic_classyfire_classification = {}

    df_ext_classyfire_cmpd_table = pd.read_csv(filename_ext_classyfire_compound_table, delimiter='\t', index_col=False)
    df_ext_classyfire_cmpd_table = df_ext_classyfire_cmpd_table.fillna("no_classification")

    path_t3db_xml = 'visualizer/data/t3db_xml/'  # TODO: Temporary T3DB XML directory
    if not os.path.exists('./visualizer/data/t3db_xml/t3db.xml'):
        with zipfile.ZipFile('./visualizer/data/t3db_xml/t3db.zip') as zf:
            zf.extractall('./visualizer/data/t3db_xml')

    list_t3db_inchikey_complete = []
    list_t3db_inchikey_head = []

    # Use cache of list_t3db_obj
    list_t3db_obj = cache.get('list_t3db_obj')
    if not list_t3db_obj:
        list_t3db_obj = read_t3db_a1.read_t3db_a1_interparse(path_t3db_xml)
        cache.set('list_t3db_obj', list_t3db_obj, 60 * 10)

    for t3db_obj in list_t3db_obj:
        list_t3db_inchikey_complete.append(t3db_obj.inchi_key)
        list_t3db_inchikey_head.append(t3db_obj.inchi_key.split('-')[0])

    count = -1
    for key, row in df_ext_classyfire_cmpd_table.iterrows():
        valid_entry = 1
        count += 1
        inchi_key = row['INCHI_KEY']

        # l_cluster_total_input_idx.append(str(count))

        ############################################
        # make new entry------------------------------
        dic_for_entry = {}
        # obj_clst = spec_cluster_a1.Spec_cluster()

        ############################################
        # get info from the file-----------------
        dic_for_entry['CMPD_CLASSIFICATION_KINGDOM'] = [row['CMPD_CLASSIFICATION_KINGDOM']]
        dic_for_entry['CMPD_CLASSIFICATION_SUPERCLASS'] = [row['CMPD_CLASSIFICATION_SUPERCLASS']]
        dic_for_entry['CMPD_CLASSIFICATION_CLASS'] = [row['CMPD_CLASSIFICATION_CLASS']]
        dic_for_entry['CMPD_CLASSIFICATION_SUBCLASS'] = [row['CMPD_CLASSIFICATION_SUBCLASS']]

        dic_inchikey_vs_dic_classyfire_classification[inchi_key] = dic_for_entry

        dic_inchikey_1st_part_vs_dic_classyfire_classification[inchi_key.split("-")[0]] = dic_for_entry

    logger.debug("  ---read_external_cmpd_info---- finished reading external CLASS file")

    fo_noclass = open("noclass_inchikey_inchi.txt", "w")

    # complement classification info for input clusters.
    # iterate cluster/node
    for input_idx, cluster_info in dic_cluster_total_input_idx_vs_cluster_info.items():

        # print "#", input_idx, cluster_info.inchi_key
        # if the current cluster (actually its represen uni) does not have superclass info
        l_tag_noclassification = ['noclassification', 'no_classification', 'NO_CLASSIFICATION']

        # do this only for spec cluster whose identity and inchi key is known
        if len(cluster_info["inchi_key"]) > 4:
            # if the spectrum has classification info
            if len(cluster_info["represen_spec_uni"]["list_cmpd_classification_superclass"]) > 0:
                if cluster_info["represen_spec_uni"]["list_cmpd_classification_superclass"][
                    0] in l_tag_noclassification:
                    # print cluster_info.inchi_key, "class info lacks", cluster_info.represen_spec_uni.list_cmpd_classification_superclass
                    dic_classyfire_classification = {}

                    # inchi key match mode : 1  only match 1st part of inchikey
                    if inchi_key_match_mode == 1:
                        # if inchikey is present in the dictionary of external info.
                        if cluster_info["inchi_key"].split("-")[
                            0] in dic_inchikey_1st_part_vs_dic_classyfire_classification:
                            dic_classyfire_classification = dic_inchikey_1st_part_vs_dic_classyfire_classification[
                                cluster_info["inchi_key"].split("-")[0]]
                            # print "dic_classyfire_classification", dic_classyfire_classification
                            cluster_info["represen_spec_uni"]["list_cmpd_classification_kingdom"] = \
                                dic_classyfire_classification['CMPD_CLASSIFICATION_KINGDOM']
                            cluster_info["represen_spec_uni"]["list_cmpd_classification_superclass"] = \
                                dic_classyfire_classification['CMPD_CLASSIFICATION_SUPERCLASS']
                            cluster_info["represen_spec_uni"]["list_cmpd_classification_class"] = \
                                dic_classyfire_classification['CMPD_CLASSIFICATION_CLASS']
                            cluster_info["represen_spec_uni"]["list_cmpd_classification_subclass"] = \
                                dic_classyfire_classification['CMPD_CLASSIFICATION_SUBCLASS']

                        if cluster_info["inchi_key"].split("-")[0] in list_t3db_inchikey_head:
                            cluster_info['is_toxic'] = True

                        else:
                            fo_noclass.write(cluster_info["inchi_key"] + "\t" + cluster_info["inchi"] + "\n")
                            logger.info(f'No classification data: {cluster_info["inchi_key"]}, {cluster_info["inchi"]}')

                    else:
                        if cluster_info["inchi_key"] in list_t3db_inchikey_complete:
                            cluster_info['is_toxic'] = True

                    # print "complemented sperclass", cluster_info.represen_spec_uni.list_cmpd_classification_superclass
                """
                for t3db_obj in list_t3db_obj:
                    f_match = 0
                    if inchi_key_match_mode == 0 :
                        if cluster_info.inchi_key == t3db_obj.inchi_key:
                """

    fo_noclass.close()

    path = dic_config["foldername_ext_cmpd_info"]

    ##############################
    # NON T3DB mode
    if path != "t3db_xml":
        # list_t3db_obj = read_t3db_a1.read_t3db_a1_interparse(path)

        # key is inchikey and value is dictionary where all external compound info is stored. key is header. value is value specified in tsv
        dic_inchikey_vs_dic_external_compound_info = {}

        for file in glob(path + '\*.tsv'):

            df_ext_compound_table = pd.read_csv(file, delimiter='\t', index_col=False)
            # df_ext_classyfire_cmpd_table = df_ext_classyfire_cmpd_table.fillna("no_classification")
            count = -1
            for key, row in df_ext_compound_table.iterrows():
                valid_entry = 1
                count += 1
                inchi_key = row['INCHI_KEY']
                skin_reactivity = row['SKIN_REACT']

                if not inchi_key in dic_inchikey_vs_dic_external_compound_info:
                    dic_inchikey_vs_dic_external_compound_info[inchi_key] = {}
                    dic_inchikey_vs_dic_external_compound_info[inchi_key]['SKIN_REACT'] = skin_reactivity

                if inchi_key in dic_inchikey_vs_dic_external_compound_info:
                    dic_inchikey_vs_dic_external_compound_info[inchi_key]['SKIN_REACT'] = skin_reactivity

        for input_idx, cluster_info in dic_cluster_total_input_idx_vs_cluster_info.items():
            l_compound_category_for_this_node = []

            # cluster_info.inchi_key

            for inchikey_ext, dic_external_compound_info in dic_inchikey_vs_dic_external_compound_info.items():

                f_match = 0
                if inchi_key_match_mode == 0:
                    if cluster_info.inchi_key == inchikey_ext:
                        f_match = 1

                if inchi_key_match_mode == 1:
                    if cluster_info.inchi_key.split("-")[0] == inchikey_ext.split("-")[0]:
                        f_match = 1

                if f_match == 1:
                    # count_t3db_found_compounds = count_t3db_found_compounds + 1

                    for k, v in dic_external_compound_info.items():
                        # l_toxin_category_for_this_node = l_toxin_category_for_this_node + t3db_obj.list_categories
                        l_compound_category_for_this_node.append(k)
            cluster_info["list_compound_categories"] = list(set(l_compound_category_for_this_node))

    ##############################
    # T3DB mode
    if path.endswith("t3db_xml/"):  # TODO: if path == "t3db_xml" -> if path.endswith("t3db_xml")
        # list_t3db_obj = read_t3db_a1.read_t3db_a1_interparse(path)

        fo.write("length of list_t3db_obj )" + str(len(list_t3db_obj)))
        logger.info(f'Length of list_t3db_obj: {len(list_t3db_obj)}')

        # iterate cluster/node
        count_t3db_found_compounds = 0
        for input_idx, cluster_info in dic_cluster_total_input_idx_vs_cluster_info.items():
            l_toxin_category_for_this_node = []
            for t3db_obj in list_t3db_obj:
                f_match = 0
                if inchi_key_match_mode == 0:
                    if cluster_info['inchi_key'] == t3db_obj.inchi_key:
                        f_match = 1

                if inchi_key_match_mode == 1:
                    if cluster_info['inchi_key'].split("-")[0] == t3db_obj.inchi_key.split("-")[0]:
                        f_match = 1

                if f_match == 1:
                    count_t3db_found_compounds = count_t3db_found_compounds + 1
                    l_toxin_category_for_this_node = l_toxin_category_for_this_node + t3db_obj.list_categories

            cluster_info['list_compound_categories'] = list(set(l_toxin_category_for_this_node))

        fo.write("count_t3db_found_compounds )" + str(count_t3db_found_compounds))
        logger.warning(f'count_t3db_found_compounds: {count_t3db_found_compounds}')
        fo.flush()

    logger.debug("--read_external_cmpd_info--  finished")


def read_classyfire_table(path):
    df = pd.read_csv(path, delimiter='\t', index_col=False)

    # Remove row if it has no 'INCHI_KEY'
    df = df[~pd.isna(df['INCHI_KEY'])]
    df = df.fillna("no_classification")
    df['CMPD_CLASSIFICATION_KINGDOM'] = df['CMPD_CLASSIFICATION_KINGDOM'].apply(lambda x: [x])
    df['CMPD_CLASSIFICATION_SUPERCLASS'] = df['CMPD_CLASSIFICATION_SUPERCLASS'].apply(lambda x: [x])
    df['CMPD_CLASSIFICATION_CLASS'] = df['CMPD_CLASSIFICATION_CLASS'].apply(lambda x: [x])
    df['CMPD_CLASSIFICATION_SUBCLASS'] = df['CMPD_CLASSIFICATION_SUBCLASS'].apply(lambda x: [x])
    # Remove duplicated inchikeys.
    df.drop_duplicates(inplace=True, subset=['INCHI_KEY'])
    # Set 'INCHI_KEY' to index
    df.set_index('INCHI_KEY')
    # Convert DataFrame to dict
    dic_inchikey_vs_dic_classyfire_classification = df.to_dict(orient='index')

    df.reset_index(inplace=True)
    # Create column for first part of inchikey.
    df['inchikey_first'] = df['INCHI_KEY'].apply(lambda x: x.split('-')[0])
    # Remove duplicated inchikey first parts.
    df.drop_duplicates(inplace=True, subset=['inchikey_first'])
    # Set 'inchikey_first' to index
    df.set_index('inchikey_first')
    df = df[['CMPD_CLASSIFICATION_KINGDOM',
             'CMPD_CLASSIFICATION_SUPERCLASS',
             'CMPD_CLASSIFICATION_CLASS',
             'CMPD_CLASSIFICATION_SUBCLASS']]
    # Convert DataFrame to dict
    dic_inchikey_1st_part_vs_dic_classyfire_classification = df.to_dict(orient='index')

    return (dic_inchikey_vs_dic_classyfire_classification,
            dic_inchikey_1st_part_vs_dic_classyfire_classification)


def _read_t3db_xml():
    path_t3db_xml = 'visualizer/data/t3db_xml/'  # TODO: Temporary T3DB XML directory
    if not os.path.exists('./visualizer/data/t3db_xml/t3db.xml'):
        with zipfile.ZipFile('./visualizer/data/t3db_xml/t3db.zip') as zf:
            zf.extractall('./visualizer/data/t3db_xml')

    dic_inchikey_vs_t3db_obj = {}
    dic_inchikey_1st_part_vs_t3db_obj = {}

    # Use cache of list_t3db_obj
    list_t3db_obj = cache.get('list_t3db_obj')
    if not list_t3db_obj:
        list_t3db_obj = read_t3db_a1.read_t3db_a1_interparse(path_t3db_xml)
        cache.set('list_t3db_obj', list_t3db_obj, 60 * 10)

    for t3db_obj in list_t3db_obj:
        # TODO: Check duplicated inchikeys
        dic_inchikey_vs_t3db_obj[t3db_obj.inchi_key] = t3db_obj
        dic_inchikey_1st_part_vs_t3db_obj[t3db_obj.inchi_key.split('-')[0]] = t3db_obj

    return dic_inchikey_vs_t3db_obj, dic_inchikey_1st_part_vs_t3db_obj


def add_external_cmpd_info(dic_cluster_total_input_idx_vs_cluster_info, folder_path_ext_cmpd_info):
    logger.info(f'Start {sys._getframe().f_code.co_name}()')
    fo = open("read_external_cmpd_info.txt", "w")

    #####################
    # inchi key match mode
    # 0: exact match
    # 1: XXXX-YYY-ZZZ match to XXXX-AAA-BBB
    inchi_key_match_mode = 1

    ######################################################
    # external compound CLASSYFIRE classification file.
    ######################################################
    filepath_ext_classyfire_compound_table = "visualizer/data/cmpd_classification_classyfire.tsv"
    (dic_inchikey_vs_dic_classyfire_classification,
     dic_inchikey_1st_part_vs_dic_classyfire_classification) = read_classyfire_table(
        filepath_ext_classyfire_compound_table)

    dic_inchikey_vs_t3db_obj, dic_inchikey_1st_part_vs_t3db_obj = _read_t3db_xml()

    logger.debug("  ---read_external_cmpd_info--- finished reading external CLASS file")

    fo_noclass = open("noclass_inchikey_inchi.txt", "w")

    # complement classification info for input clusters.
    # iterate cluster/node
    l_tag_noclassification = ['noclassification', 'no_classification', 'NO_CLASSIFICATION']
    for input_idx, cluster_info in dic_cluster_total_input_idx_vs_cluster_info.items():
        # do this only for spec cluster whose identity and inchi key is known
        if len(cluster_info['inchi_key']) < 14:
            continue

        # if the spectrum has classification info
        if (len(cluster_info['represen_spec_uni']['list_cmpd_classification_superclass']) == 0
                or cluster_info['represen_spec_uni']['list_cmpd_classification_superclass'][0] in l_tag_noclassification):

            # inchi key match mode : 1 (only match 1st part of inchikey)
            if inchi_key_match_mode == 1:
                # if inchikey is present in the dictionary of external info.
                if cluster_info['inchi_key_first'] in dic_inchikey_1st_part_vs_dic_classyfire_classification:
                    dic_classyfire_classification = \
                        dic_inchikey_1st_part_vs_dic_classyfire_classification[cluster_info['inchi_key_first']]

                    # Fill kingdom
                    cluster_info['represen_spec_uni']['list_cmpd_classification_kingdom'] = \
                        dic_classyfire_classification['CMPD_CLASSIFICATION_KINGDOM']
                    # Fill superclass
                    cluster_info['represen_spec_uni']['list_cmpd_classification_superclass'] = \
                        dic_classyfire_classification['CMPD_CLASSIFICATION_SUPERCLASS']
                    # Fill class
                    cluster_info['represen_spec_uni']['list_cmpd_classification_class'] = \
                        dic_classyfire_classification['CMPD_CLASSIFICATION_CLASS']
                    # Fill subclass
                    cluster_info['represen_spec_uni']['list_cmpd_classification_subclass'] = \
                        dic_classyfire_classification['CMPD_CLASSIFICATION_SUBCLASS']
                else:
                    fo_noclass.write(cluster_info['inchi_key'] + '\t' + cluster_info['inchi'] + '\n')
                    logger.info(f'No classification data: {cluster_info["inchi_key"]}, {cluster_info["inchi"]}')

                if cluster_info['inchi_key_first'] in dic_inchikey_1st_part_vs_t3db_obj:
                    cluster_info['is_toxic'] = True

            # inchi key match mode : 0 (exact match)
            else:
                if cluster_info['inchi_key'] in dic_inchikey_vs_dic_classyfire_classification:
                    dic_classyfire_classification = dic_inchikey_vs_dic_classyfire_classification[cluster_info['inchi_key']]

                    # Fill kingdom
                    cluster_info['represen_spec_uni']['list_cmpd_classification_kingdom'] = \
                        dic_classyfire_classification['CMPD_CLASSIFICATION_KINGDOM']
                    # Fill superclass
                    cluster_info['represen_spec_uni']['list_cmpd_classification_superclass'] = \
                        dic_classyfire_classification['CMPD_CLASSIFICATION_SUPERCLASS']
                    # Fill class
                    cluster_info['represen_spec_uni']['list_cmpd_classification_class'] = \
                        dic_classyfire_classification['CMPD_CLASSIFICATION_CLASS']
                    # Fill subclass
                    cluster_info['represen_spec_uni']['list_cmpd_classification_subclass'] = \
                        dic_classyfire_classification['CMPD_CLASSIFICATION_SUBCLASS']

                else:
                    fo_noclass.write(cluster_info['inchi_key'] + '\t' + cluster_info['inchi'] + '\n')
                    logger.info(f'No classification data: {cluster_info["inchi_key"]}, {cluster_info["inchi"]}')

                if cluster_info['inchi_key'] in dic_inchikey_vs_t3db_obj:
                    cluster_info['is_toxic'] = True

    fo_noclass.close()

    # Add list_compound_categories data to cluster info
    if not folder_path_ext_cmpd_info:
        return

    ##############################
    # NON T3DB mode
    if not folder_path_ext_cmpd_info.endswith('t3db_xml/'):
        # key is inchikey and value is dictionary where all external compound info is stored.
        # key is header. value is value specified in tsv
        dic_inchikey_vs_dic_external_compound_info = {}
        dic_inchikey_1st_part_vs_dic_external_compound_info = {}
        paths = []
        for f in os.listdir(folder_path_ext_cmpd_info):
            p = os.path.join(folder_path_ext_cmpd_info, f)
            if os.path.isfile(p) and f.endswith('.tsv'):
                paths.append(p)

        for path in paths:
            df_ext_compound_table = pd.read_csv(path, delimiter='\t', index_col=False)

            if 'INCHI_KEY' not in df_ext_compound_table.columns or 'SKIN_REACT' not in df_ext_compound_table.columns:
                continue

            for key, row in df_ext_compound_table.iterrows():
                inchi_key = row['INCHI_KEY']
                inchi_key_first = row['INCHI_KEY'].split('-')[0]
                skin_reactivity = row['SKIN_REACT']

                if inchi_key in dic_inchikey_vs_dic_external_compound_info:
                    dic_inchikey_vs_dic_external_compound_info[inchi_key]['SKIN_REACT'] = skin_reactivity
                else:
                    dic_inchikey_vs_dic_external_compound_info[inchi_key] = {'SKIN_REACT': skin_reactivity}

                if inchi_key_first in dic_inchikey_1st_part_vs_dic_external_compound_info:
                    dic_inchikey_1st_part_vs_dic_external_compound_info[inchi_key_first]['SKIN_REACT'] = skin_reactivity
                else:
                    dic_inchikey_1st_part_vs_dic_external_compound_info[inchi_key_first] = {'SKIN_REACT': skin_reactivity}

        for input_idx, cluster_info in dic_cluster_total_input_idx_vs_cluster_info.items():
            if inchi_key_match_mode == 0:  # Exact match
                if cluster_info['inchi_key'] in dic_inchikey_vs_dic_external_compound_info:
                    dic_external_compound_info = dic_inchikey_vs_dic_external_compound_info[cluster_info['inchi_key']]
                    list_compound_categories = list(set(cluster_info["list_compound_categories"]
                                                        + list(dic_external_compound_info.keys())))

                    cluster_info["list_compound_categories"] = list_compound_categories

            else:
                if cluster_info['inchi_key_first'] in dic_inchikey_1st_part_vs_dic_external_compound_info:
                    dic_external_compound_info = dic_inchikey_1st_part_vs_dic_external_compound_info[cluster_info['inchi_key_first']]
                    list_compound_categories = list(set(cluster_info["list_compound_categories"]
                                                        + list(dic_external_compound_info.keys())))
                    cluster_info["list_compound_categories"] = list_compound_categories

    ##############################
    # T3DB mode
    else:
        fo.write("length of list_t3db_obj )" + str(len(dic_inchikey_vs_t3db_obj)))
        logger.info(f'Length of list_t3db_obj: {len(dic_inchikey_vs_t3db_obj)}')

        # iterate cluster/node
        count_t3db_found_compounds = 0
        for input_idx, cluster_info in dic_cluster_total_input_idx_vs_cluster_info.items():
            l_toxin_category_for_this_node = []

            if inchi_key_match_mode == 0:  # Exact match
                if cluster_info['inchi_key'] in dic_inchikey_vs_t3db_obj:
                    count_t3db_found_compounds += 1
                    t3db_obj = dic_inchikey_vs_t3db_obj[cluster_info['inchi_key']]
                    cluster_info['list_compound_categories'] = t3db_obj.list_categories
            else:
                if cluster_info['inchi_key_first'] in dic_inchikey_1st_part_vs_t3db_obj:
                    count_t3db_found_compounds += 1
                    t3db_obj = dic_inchikey_1st_part_vs_t3db_obj[cluster_info['inchi_key_first']]
                    cluster_info['list_compound_categories'] = t3db_obj.list_categories

        fo.write("count_t3db_found_compounds )" + str(count_t3db_found_compounds))
        logger.warning(f'count_t3db_found_compounds: {count_t3db_found_compounds}')
        fo.flush()

    logger.debug("  --read_external_cmpd_info--  finished")
