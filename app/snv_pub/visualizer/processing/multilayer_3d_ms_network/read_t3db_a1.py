# last modified 20200614



#######
# usage
"""
# path where (only) t3db xml is located.
path = 't3db_source'
list_hmdb_obj = read_t3db_a1.read_t3db_a1_interparse(path)
print "--output---"
for hmdb_obj in (list_hmdb_obj):
    print "\n---\n HMDB accession" , hmdb_obj.id
    print "hmdb_name, ", hmdb_obj.name
    print "inchi, ", hmdb_obj.inchi
    print "inchikey, ", hmdb_obj.inchikey
    print "categories, " , hmdb_obj.list_categories
    print "types, " , hmdb_obj.list_types

"""

import glob
import re
from logging import getLogger
import math
import csv
import copy
import sys

import numpy as np
import os
import csv
import matplotlib.pyplot as plt
# encoding: utf-8
from xml.etree import ElementTree


# a3 version ---------------------------------------------------
# now it reads XML file and returns list of compound info

logger = getLogger(__name__)


"""

class std_cmpd_class:
    def __init__(self):
            self.name    =   ""
            self.id    =    0
            self.mix    =    ""
            self.cas_num    =    ""
            self.ele_composition    =   ""
            self.exact_mass    =    ""
            self.mz    =    0
            self.lod = 0
            self.mix =""
            self.hmdb_id =""
            self.smiles =""
#   id	mix	casno	name	exact mass	mz
#   1	1	100-52	serotonin	100.2	101.2
"""


class db_toxin_cmpd_class:
    def __init__(self):
            self.name    =   ""
            self.toxin_class = ""
            self.id    =    0
            self.mix    =    ""
            self.cas_num    =    ""
            self.pubchem_cid    =    ""
            self.chemical_formula    =   ""
            self.monoiso_mol_weight    =   0
            self.mz    =    0
            self.lod = 0
            self.mix =""
            self.hmdb_id =""
            self.smiles =""
            self.inchi = ""
            self.inchi_key = ""
            self.exact_mass = 0.0
            self.tax_direct_parent=""
            self.tax_kingdom =""
            self.tax_super_class =""
            self.tax_class =""
            self.tax_sub_class =""
            self.tax_molecular_framework =""


            self.list_categories = []
            self.list_types = []


            self.list_tax_alternative_parent = []
            self.list_tax_substituent = []
            self.list_origins = []
            self.flag_endogenous = 0

            self.experimental_properties_water_solubility = ""

            self.predicted_properties_solubility = ""

            self.list_binary_vector_nominal_1k = []




list_spec_obj =[]

list_db_cmpd_info_obj=[]


##############
#  this interparse version is for large XML file
###############
def read_t3db_a1_interparse(xml_path):

    """

    :type xml_path: object
    """
    list_hmdb_obj =[]
    logger.warning(f'@@@ File in xml_path: {os.listdir(xml_path)}')
    files = [os.path.join(xml_path, filename) for filename in os.listdir(xml_path) if filename.endswith('.xml')]
    logger.warning(f"@@@ {files}")
    for file in files:  # TODO: glob() is not working, use os.listdir()
        logger.warning(f'read_t3db_a1_interparse: {file}')
        f = open(file,"r",encoding="utf-8")

        myregex_XML_ext = re.compile(".XML"  , re.IGNORECASE)

        match =  myregex_XML_ext.search(file)

        print(file)
        if not (match) :
            print("from readhmdb_a3. exiting because reading non XML file")
            sys.exit()

        if (match):
            print("--read_t3db_a1_interparse-- reading hmdb file")

            #print match.group(1)




            context = ElementTree.iterparse(f, events=('start', 'end'))
            #tree = ElementTree.parse(f)

            #root = context.getroot()
            _, root = next(context)

            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # define namespace
            # """""""""""""""""""""""""""""
            #namespace = r'{http://www.hmdb.ca}'
            namespace = ""

            count = 0

            print("--read_t3db_a1_interparse--  starting iterate for metabolites")



            for event, elem in context:
                #print "elem-TTTAGG" , elem.tag
                if event == 'end' and elem.tag == namespace+'compound':
                #if event == 'end' and elem.tag == 'compound':
                    count = count +1

                    if( count%1000 == 0): print(count, end=' ')
                    #if (count  == 5): sys.exit()


                    current_toxin_cmpd = db_toxin_cmpd_class()
                    list_tax = []
                    for name in elem.findall( namespace + 'common_name'):
                        #print ""
                        current_toxin_cmpd.name = name.text
                        #print "[common_name]", current_toxin_cmpd.name

                    for toxin_class in elem.findall(namespace + 'class'):
                        current_toxin_cmpd.toxin_class = toxin_class.text
                        #print "[class]",current_toxin_cmpd.toxin_class


                    for categories in elem.findall(namespace + 'categories'):
                        for x in categories.findall(namespace + 'category'):


                            # current_toxin_cmpd.list_tax_alternative_parent.append(x.text)
                            current_toxin_cmpd.list_categories.append(x.text)

                    for types in elem.findall(namespace + 'types'):
                        for x in types.findall(namespace + 'type'):

                            # current_toxin_cmpd.list_tax_alternative_parent.append(x.text)
                            current_toxin_cmpd.list_types.append(x.text)


                    for accession in elem.findall( namespace + 'accession'):
                        current_toxin_cmpd.id = accession.text
                        #current_toxin_cmpd.id
                        #print current_toxin_cmpd.id
                    for cas in elem.findall( namespace + 'cas_registry_number'):

                        if cas.text != None :
                            current_toxin_cmpd.cas_num = cas.text


                    for monoiso_mol_weight in elem.findall( namespace + 'monisotopic_molecular_weight'):
                        #print monoiso_mol_weight  , type(monoiso_mol_weight)
                        if monoiso_mol_weight.text is not None:
                            current_toxin_cmpd.monoiso_mol_weight = float(monoiso_mol_weight.text)

                    for chemical_formula in elem.findall( namespace + 'chemical_formula'):
                        current_toxin_cmpd.chemical_formula = chemical_formula.text

                    for smiles in elem.findall(namespace +'smiles'):

                        str_smiles = smiles.text
                        current_toxin_cmpd.smiles =str_smiles

                    for inchi in elem.findall(namespace +'inchi'):
                        current_toxin_cmpd.inchi = inchi.text

                    for inchikey in elem.findall(namespace +'inchikey'):
                        #print "inchikey.text" , inchikey.text
                        if inchikey.text != None :
                            if inchikey.text.startswith("InChIKey=") :
                                current_toxin_cmpd.inchi_key = inchikey.text.split("=")[1]



                        """
                        pattern = r"="
                        repatter = re.compile(pattern)


                        matchOB = repatter.search(inchikey.text)
                        if matchOB :
                            current_hmdb_cmpd.inchikey = inchikey.text.split("=")[1]

                        if not matchOB :
                            current_hmdb_cmpd.inchikey = inchikey.text
                        """

                    for mycid in elem.findall(namespace + 'pubchem_compound_id'):
                        str_pubchemcid  = mycid.text
                        current_toxin_cmpd.pubchem_cid = str_pubchemcid



                    """
                    # taxonomy ---------------------------------
                    for tax in elem.findall(namespace + 'taxonomy'):

                        for tax_kngdm in tax.findall(namespace + 'kingdom'):

                            list_tax.append(tax_kngdm.text)
                            if tax_kngdm.text is not None:
                                current_toxin_cmpd.tax_kingdom = tax_kngdm.text

                            #print tax_kngdm.text
                        for tax_spclass in tax.findall(namespace + 'super_class'):
                            #print "super_class:" , tax_spclass.text
                            list_tax.append(tax_spclass.text)
                            if tax_spclass.text is not None:
                                current_toxin_cmpd.tax_super_class = tax_spclass.text

                        for tax_class in tax.findall(namespace +'class'):
                            #print "class:" , tax_class.text
                            list_tax.append(tax_class.text)
                            if tax_class.text is not None:
                                current_toxin_cmpd.tax_class = tax_class.text

                        for tax_sub_class in tax.findall(namespace +'sub_class'):
                            #print "class:" , tax_class.text

                            if tax_sub_class.text is not  None :
                                current_toxin_cmpd.tax_sub_class = tax_sub_class.text

                        for tax_directparent in tax.findall(namespace +'direct_parent'):
                            #print "direct_parent:" , tax_directparent.text
                            list_tax.append(tax_directparent.text)
                            current_toxin_cmpd.tax_direct_parent = tax_directparent.text
                        for mf in tax.findall(namespace +'molecular_framework'):
                            current_toxin_cmpd.tax_molecular_framework = mf.text
                        for tax_substituent in tax.findall(namespace +'substituents'):
                            for x in tax_substituent.findall(namespace +'substituent'):
                                current_toxin_cmpd.list_tax_substituent.append(x.text)

                        for tax_alternative_parents in tax.findall(namespace + 'alternative_parents'):
                            for x in tax_alternative_parents.findall(namespace + 'alternative_parent'):
                                current_toxin_cmpd.list_tax_alternative_parent.append(x.text)
                    #-------------------------------------------------

                    """

                    """
                    # experimental_properties ---------------------------------
                    for experimental_properties in elem.findall(namespace + 'experimental_properties'):

                        for property in experimental_properties.findall(namespace + 'property'):
                            f_ws = 0
                            for kind in property.findall(namespace + 'kind'):
                                if kind.text == "water_solubility":
                                    f_ws = 1
                                    #print "WATER SOLLLLLLLLLLLLLL"
                                    #print property.text
                            for value in property.findall(namespace + 'value'):
                                if f_ws == 1 :
                                    #print "WATER S:" , value.text

                                    current_toxin_cmpd.experimental_properties_water_solubility = value.text
                    """
                    """
                    # predicted_properties ---------------------------------
                    for predicted_properties in elem.findall(namespace + 'predicted_properties'):

                        for property in predicted_properties.findall(namespace + 'property'):
                            f_ws = 0
                            for kind in property.findall(namespace + 'kind'):
                                if kind.text == "solubility":
                                    f_ws = 1
                                    #print "WATER SOLLLLLLLLLLLLLL"
                                    #print property.text
                            for value in property.findall(namespace + 'value'):
                                if f_ws == 1 :
                                    #print "WATER PREDICTED S:" , value.text

                                    current_toxin_cmpd.predicted_properties_solubility = value.text

                    """


                    """
                    #Substituent -----------------------------------------
                    list_substituent =[]
                    for substituent in elem.findall('taxonomy/substituents/substituent'):
                        #print "substituent:" , substituent.text
                        list_substituent.append(substituent.text)
                    current_hmdb_cmpd.list_tax_substituent = list_substituent
                    """
                    """
                    # alternative parent--------------------------
                    list_alternative_parent =[]
                    for alternative_parent in elem.findall('taxonomy/alternative_parents/alternative_parent'):
                        print "sualternative_parent:" , alternative_parent.text
                        list_alternative_parent.append(alternative_parent.text)
                    current_hmdb_cmpd.list_tax_alternative_parent = list_alternative_parent
                    ###########
                    """

                    list_origins = []
                    flag_endogenous = 0
                    for origin in elem.findall(namespace +'ontology/origins/origin'):


                        list_origins.append(origin.text)
                        if( origin.text == 'Endogenous'  or origin.text == 'endogenous'    ):
                            flag_endogenous =1

                            current_toxin_cmpd.flag_endogenous = 1
                            current_toxin_cmpd.list_origins = list_origins

                    #tax
                    #current_hmdb_cmpd.cas_num = get_cas_number_from_hmdb_xml(file)
                    #current_hmdb_cmpd.pubchem_cid    =  get_pubchemcid_from_hmdb_xml(file)
                    #current_hmdb_cmpd.hmdb_id =  get_hmdb_id_from_hmdb_xml(file)
                    #current_hmdb_cmpd.tax_super_class =    get_chem_taxonomy_info_list_from_hmdb_xml(file)[1]
                    #current_hmdb_cmpd.tax_class =    get_chem_taxonomy_info_list_from_hmdb_xml(file)[2]
                    #current_hmdb_cmpd.smiles  = get_smiles_from_hmdb_xml(file)

                    list_hmdb_obj.append(current_toxin_cmpd)
                    root.clear()
    return list_hmdb_obj














def readhmdb_a3(xml_path):

    list_hmdb_obj =[]
    #for file in glob.glob('E:\work\hmdb_metabolites\*.xml'):
    for file in glob.glob(xml_path + '\*.xml'):
    #for file in glob.glob('C:\\Users\\esk\\Dropbox_eisuke.hayakawaATgmail\\Dropbox\\Programming_drpbx\\python\\mass_spec_related\\chemclass_vs_nominal_delta\\xml\\*.xml'):
        f = open(file,"r")



        myregex_XML_ext = re.compile(".XML"  , re.IGNORECASE)

        match =  myregex_XML_ext.search(file)

        print(file)
        if not (match) :
            print("from readhmdb_a3. exiting because reading non XML file")
            sys.exit()

        if (match):
            print("reading hmdb file")

            #print match.group(1)


            print(file)
            print(f)
            print("parsing ElementTree")
            tree = ElementTree.parse(f)
            print("getting root")
            root = tree.getroot()


            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # define namespace
            # """""""""""""""""""""""""""""
            namespace = r'{http://www.hmdb.ca}'


            print("nama" , namespace)





            count = 0

            print("starting iterate for metabolites")

            for metabolite in root.findall(namespace + 'metabolite', ):


                count = count +1
                if( count%100 == 0): print(count, end=' ')

                current_hmdb_cmpd = db_cmpd_class()
                list_tax = []

                for accession in metabolite.findall( namespace + 'accession'):

                    current_hmdb_cmpd.id = accession.text
                for cas in metabolite.findall( namespace + 'cas_registry_number'):
                    current_hmdb_cmpd.cas_num = cas.text


                for monoiso_mol_weight in metabolite.findall( namespace + 'monisotopic_moleculate_weight'):
                    current_hmdb_cmpd.monoiso_mol_weight = float(monoiso_mol_weight.text)

                for smiles in metabolite.findall(namespace +'smiles'):

                    str_smiles = smiles.text
                    current_hmdb_cmpd.smiles =str_smiles

                for inchi in metabolite.findall(namespace +'inchi'):

                    current_hmdb_cmpd.inchi = inchi.text

                for mycid in metabolite.findall(namespace + 'pubchem_compound_id'):
                    str_pubchemcid  = mycid.text
                    current_hmdb_cmpd.pubchem_cid = str_pubchemcid

                # taxonomy ---------------------------------
                for tax in metabolite.findall(namespace + 'taxonomy'):


                    for tax_kngdm in tax.findall(namespace + 'kingdom'):

                        list_tax.append(tax_kngdm.text)
                        current_hmdb_cmpd.tax_kingdom = tax_kngdm.text

                        print(tax_kngdm.text)
                    for tax_spclass in tax.findall(namespace + 'super_class'):
                        #print "super_class:" , tax_spclass.text
                        list_tax.append(tax_spclass.text)
                        current_hmdb_cmpd.tax_super_class = tax_spclass.text

                    for tax_class in tax.findall(namespace +'tclass'):
                        #print "class:" , tax_class.text
                        list_tax.append(tax_class.text)
                        current_hmdb_cmpd.tax_class = tax_class.text

                    for tax_directparent in tax.findall(namespace +'direct_parent'):
                        #print "direct_parent:" , tax_directparent.text
                        list_tax.append(tax_directparent.text)
                        current_hmdb_cmpd.tax_direct_parent = tax_directparent.text

                list_origins = []
                flag_endogenous = 0
                for origin in metabolite.findall(namespace +'ontology/origins/origin'):


                    list_origins.append(origin.text)
                    if( origin.text == 'Endogenous'  or origin.text == 'endogenous'    ):
                        flag_endogenous =1

                        current_hmdb_cmpd.flag_endogenous = 1
                current_hmdb_cmpd.list_origins = list_origins

                #-------------------------------------------------


                #Substituent -----------------------------------------
                list_substituents =[]
                for substituent in root.findall('taxonomy/substituents/substituent'):
                    #print "substituent:" , substituent.text
                    list_substituents.append(substituent.text)
                current_hmdb_cmpd.list_substituents = list_substituents

                #tax
                #current_hmdb_cmpd.cas_num = get_cas_number_from_hmdb_xml(file)
                #current_hmdb_cmpd.pubchem_cid    =  get_pubchemcid_from_hmdb_xml(file)
                #current_hmdb_cmpd.hmdb_id =  get_hmdb_id_from_hmdb_xml(file)
                #current_hmdb_cmpd.tax_super_class =    get_chem_taxonomy_info_list_from_hmdb_xml(file)[1]
                #current_hmdb_cmpd.tax_class =    get_chem_taxonomy_info_list_from_hmdb_xml(file)[2]
                #current_hmdb_cmpd.smiles  = get_smiles_from_hmdb_xml(file)

                #print current_hmdb_cmpd.hmdb_id
                #print current_hmdb_cmpd.smiles
                list_hmdb_obj.append(current_hmdb_cmpd)

    return list_hmdb_obj
















##################################################


def get_smiles_from_hmdb_xml(str_file_name):
    str_smiles = " "
    page = {}
#    f = open("test.xml", "r")

    #with open('hmdb_xml\HMDB00001.XML', 'rt') as f:
    with open(str_file_name, 'rt') as f:

        tree = ElementTree.parse(f)
        root = tree.getroot()
        #print "roottag" , root.tag

        for smiles in root.findall('smiles'):
            str_smiles = smiles.text

            #print smiles.tag
            #print smiles.text

    print("str_smiles : " , str_smiles)
    if( str_smiles == None):

        str_smiles = ""
    return str_smiles


def get_cas_number_from_hmdb_xml(str_file_name):

    page = {}
#    f = open("test.xml", "r")
    str_cas = ""
    #with open('hmdb_xml\HMDB00001.XML', 'rt') as f:
    with open(str_file_name, 'rt') as f:

        tree = ElementTree.parse(f)
        root = tree.getroot()
        #print "roottag" , root.tag

        for cas in root.findall('cas_registry_number'):
            str_cas = cas.text

            #print smiles.tag
            #print smiles.text
    return str_cas


def get_pubchemcid_from_hmdb_xml(str_file_name):
    str_pubchemcid = ""
    with open(str_file_name, 'rt') as f:

        tree = ElementTree.parse(f)
        root = tree.getroot()

        for mycid in root.findall('pubchem_compound_id'):
            str_pubchemcid  = mycid.text

            #print smiles.tag
            #print smiles.text
    return str_pubchemcid


def get_hmdb_id_from_hmdb_xml(str_file_name):
    str_hmdb_id = ""
    with open(str_file_name, 'rt') as f:

        tree = ElementTree.parse(f)
        root = tree.getroot()

        for myid in root.findall('accession'):
            str_hmdb_id  = myid.text


    return str_hmdb_id



def get_chem_taxonomy_info_list_from_hmdb_xml(str_file_name):

    page = {}
#    f = open("test.xml", "r")
    list_tax = []
    #with open('hmdb_xml\HMDB00001.XML', 'rt') as f:
    with open(str_file_name, 'rt') as f:

        tree = ElementTree.parse(f)
        root = tree.getroot()

        for tax_kngdm in root.findall('taxonomy/kingdom'):
            #print "kingdom:" , tax_kngdm.text
            list_tax.append(tax_kngdm.text)

        for tax_spclass in root.findall('taxonomy/super_class'):
            #print "super_class:" , tax_spclass.text
            list_tax.append(tax_spclass.text)
        for tax_class in root.findall('taxonomy/class'):
            #print "class:" , tax_class.text
            list_tax.append(tax_class.text)
        for tax_directparent in root.findall('taxonomy/direct_parent'):
            #print "direct_parent:" , tax_directparent.text
            list_tax.append(tax_directparent.text)


            #print smiles.tag
            #print smiles.text
    return list_tax



def get_all_substituents_from_hmdb_xml(str_file_name):
    str_smiles =""

    with open(str_file_name, 'rt') as f:

        tree = ElementTree.parse(f)
        root = tree.getroot()

        for substituent in root.findall('taxonomy/substituents/substituent'):
            print("substituent:" , substituent.text)



    return str_smiles




# main --------------------------------------



#read hmdb all files-=========================================
#import generate_theo_fragment_a1
#import random
#import delta_generator_a1


myregex_hmdb_id = re.compile("(HMDB[0-9]+).")


list_spec_obj =[]

list_db_cmpd_info_obj=[]

def readhmdb_a2(xml_path):

    list_hmdb_obj =[]
    #for file in glob.glob('E:\work\hmdb_metabolites\*.xml'):
    for file in glob.glob(xml_path + '\*.xml'):
    #for file in glob.glob('C:\\Users\\esk\\Dropbox_eisuke.hayakawaATgmail\\Dropbox\\Programming_drpbx\\python\\mass_spec_related\\chemclass_vs_nominal_delta\\xml\\*.xml'):
        f = open(file,"r")

        match =myregex_hmdb_id.search(file)
        str_hmdb_id = ""
        if (match):
            print("reading hmdb file")
            current_hmdb_cmpd = db_cmpd_class()
            #print match.group(1)
            str_hmdb_id = match.group(1)

            #print file

            tree = ElementTree.parse(f)
            root = tree.getroot()
            #print "roottag" , root.tag


            list_tax = []
            for cas in root.findall('cas_registry_number'):
                current_hmdb_cmpd.cas_num = cas.text

            for monoiso_mol_weight in root.findall('monisotopic_moleculate_weight'):
                current_hmdb_cmpd.monoiso_mol_weight = float(monoiso_mol_weight.text)

            for smiles in root.findall('smiles'):
                str_smiles = smiles.text
                current_hmdb_cmpd.smiles =str_smiles

            for mycid in root.findall('pubchem_compound_id'):
                str_pubchemcid  = mycid.text
                current_hmdb_cmpd.pubchem_cid = str_pubchemcid

            # taxonomy ---------------------------------
            for tax_kngdm in root.findall('taxonomy/kingdom'):
                #print "kingdom:" , tax_kngdm.text
                list_tax.append(tax_kngdm.text)
                current_hmdb_cmpd.tax_kingdom = tax_kngdm.text
            for tax_spclass in root.findall('taxonomy/super_class'):
                #print "super_class:" , tax_spclass.text
                list_tax.append(tax_spclass.text)
                current_hmdb_cmpd.tax_super_class = tax_spclass.text

            for tax_class in root.findall('taxonomy/class'):
                #print "class:" , tax_class.text
                list_tax.append(tax_class.text)
                current_hmdb_cmpd.tax_class = tax_class.text

            for tax_directparent in root.findall('taxonomy/direct_parent'):
                #print "direct_parent:" , tax_directparent.text
                list_tax.append(tax_directparent.text)
                current_hmdb_cmpd.tax_direct_parent = tax_directparent.text

            list_origins = []
            flag_endogenous = 0
            for origin in root.findall('ontology/origins/origin'):

                #print origin.text
                list_origins.append(origin.text)
                if( origin.text == 'Endogenous'  or origin.text == 'endogenous'    ):
                    flag_endogenous =1
                    print("endogenous !!!!")
                    current_hmdb_cmpd.flag_endogenous = 1
            current_hmdb_cmpd.list_origins = list_origins

            #-------------------------------------------------


            #Substituent -----------------------------------------
            list_substituents =[]
            for substituent in root.findall('taxonomy/substituents/substituent'):
                #print "substituent:" , substituent.text
                list_substituents.append(substituent.text)
            current_hmdb_cmpd.list_substituents = list_substituents

            #tax
            #current_hmdb_cmpd.cas_num = get_cas_number_from_hmdb_xml(file)
            #current_hmdb_cmpd.pubchem_cid    =  get_pubchemcid_from_hmdb_xml(file)
            #current_hmdb_cmpd.hmdb_id =  get_hmdb_id_from_hmdb_xml(file)
            #current_hmdb_cmpd.tax_super_class =    get_chem_taxonomy_info_list_from_hmdb_xml(file)[1]
            #current_hmdb_cmpd.tax_class =    get_chem_taxonomy_info_list_from_hmdb_xml(file)[2]
            #current_hmdb_cmpd.smiles  = get_smiles_from_hmdb_xml(file)

            #print current_hmdb_cmpd.hmdb_id
            #print current_hmdb_cmpd.smiles
            list_hmdb_obj.append(current_hmdb_cmpd)

    return list_hmdb_obj




#------------------------------------------------------






