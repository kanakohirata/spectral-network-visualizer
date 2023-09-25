from . import spectrum_uni_a2
from rdkit import Chem

class Spec_cluster:
    def __init__(self):
        self.cluster_id = 0
        self.compound_name  = ""
        self.global_accession = ""
        self.represen_spec_uni = spectrum_uni_a2.spectrum_uni_class
        self.list_spec_uni =[]

        #compound identity related
        self.list_COMMON_NAME = []
        self.list_external_compound_UNIQUE_ID = []
        self.list_pathway_UNIQUE_ID = []
        self.list_pathway_COMMON_NAME =[]

        self.list_compound_categories =[]

        self.inchi = ""
        self.inchi_key = ""
        #self.rdkit_mol = Chem
        self.tag = ""

        self.dataset =""
        self.layer =""
        self.unique_level = -1

        #
        #    mol = Chem.MolFromInchi(OriginalInChI)
        #    my_.rdkit_mol  = mol