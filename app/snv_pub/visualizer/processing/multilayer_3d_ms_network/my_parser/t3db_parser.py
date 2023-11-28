from logging import getLogger
import os
import re
from xml.etree import ElementTree

logger = getLogger(__name__)


class ToxinCompound:
    def __init__(self):
        self.name = ''
        self.toxin_class = ''
        self.id = 0
        self.mix = ''
        self.cas_num = ''
        self.pubchem_cid = ''
        self.chemical_formula = ''
        self.monoiso_mol_weight = 0
        self.mz = 0
        self.lod = 0
        self.mix = ''
        self.hmdb_id = ''
        self.smiles = ''
        self.inchi = ''
        self.inchi_key = ''
        self.exact_mass = 0.0
        self.tax_direct_parent = ''
        self.tax_kingdom = ''
        self.tax_super_class = ''
        self.tax_class = ''
        self.tax_sub_class = ''
        self.tax_molecular_framework = ''

        self.list_categories = []
        self.list_types = []

        self.list_tax_alternative_parent = []
        self.list_tax_substituent = []
        self.list_origins = []
        self.flag_endogenous = 0

        self.experimental_properties_water_solubility = ''

        self.predicted_properties_solubility = ''

        self.list_binary_vector_nominal_1k = []


def read_t3db_xml_files(dir_path):
    """
    Read T3DB XML files.
    Parameters
    ----------
    dir_path : str
        A folder path including T3DB XML files.

    Returns
    -------
    A list of ToxinCompound instances.
    """
    paths = [os.path.join(dir_path, filename) for filename in os.listdir(dir_path)
             if filename.endswith('.xml') or filename.endswith('.XML')]

    list_toxin_compound = []
    for path in paths:
        list_toxin_compound += read_t3db_interparse(path)

    return list_toxin_compound


def read_t3db_interparse(path):
    logger.debug('-- read_t3db_interparse --')
    logger.debug(f'Read a T3DB XML file: {path}')
    f = open(path, 'r', encoding='utf-8')

    my_regex_xml_ext = re.compile('.XML', re.IGNORECASE)
    match = my_regex_xml_ext.search(path)
    if not match:
        raise ValueError('The T3BD file must be XML format.')

    context = ElementTree.iterparse(f, events=('start', 'end'))
    _, root = next(context)

    list_toxin_compound = []
    namespace = ''
    count = 0
    for event, elem in context:
        if event == 'end' and elem.tag == namespace + 'compound':
            count = count + 1

            if count % 1000 == 0:
                print(count, end=' ')

            current_toxin_cmpd = ToxinCompound()
            for name in elem.findall(namespace + 'common_name'):
                current_toxin_cmpd.name = name.text

            for toxin_class in elem.findall(namespace + 'class'):
                current_toxin_cmpd.toxin_class = toxin_class.text

            for categories in elem.findall(namespace + 'categories'):
                for x in categories.findall(namespace + 'category'):
                    current_toxin_cmpd.list_categories.append(x.text)

            for types in elem.findall(namespace + 'types'):
                for x in types.findall(namespace + 'type'):
                    current_toxin_cmpd.list_types.append(x.text)

            for accession in elem.findall(namespace + 'accession'):
                current_toxin_cmpd.id = accession.text
                
            for cas in elem.findall(namespace + 'cas_registry_number'):

                if cas.text is not None:
                    current_toxin_cmpd.cas_num = cas.text

            for monoiso_mol_weight in elem.findall(namespace + 'monisotopic_molecular_weight'):
                # print monoiso_mol_weight  , type(monoiso_mol_weight)
                if monoiso_mol_weight.text is not None:
                    current_toxin_cmpd.monoiso_mol_weight = float(monoiso_mol_weight.text)

            for chemical_formula in elem.findall(namespace + 'chemical_formula'):
                current_toxin_cmpd.chemical_formula = chemical_formula.text

            for smiles in elem.findall(namespace + 'smiles'):
                str_smiles = smiles.text
                current_toxin_cmpd.smiles = str_smiles

            for inchi in elem.findall(namespace + 'inchi'):
                current_toxin_cmpd.inchi = inchi.text

            for inchikey in elem.findall(namespace + 'inchikey'):
                # print "inchikey.text" , inchikey.text
                if inchikey.text is not None:
                    if inchikey.text.startswith("InChIKey="):
                        current_toxin_cmpd.inchi_key = inchikey.text.split("=")[1]

            for mycid in elem.findall(namespace + 'pubchem_compound_id'):
                str_pubchemcid = mycid.text
                current_toxin_cmpd.pubchem_cid = str_pubchemcid

            list_origins = []
            for origin in elem.findall(namespace + 'ontology/origins/origin'):

                list_origins.append(origin.text)
                if origin.text == 'Endogenous' or origin.text == 'endogenous':

                    current_toxin_cmpd.flag_endogenous = 1
                    current_toxin_cmpd.list_origins = list_origins

            list_toxin_compound.append(current_toxin_cmpd)
            root.clear()
    
    return list_toxin_compound
