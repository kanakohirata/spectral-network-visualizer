from django.core.exceptions import ValidationError
from django.db import models
from django.utils import timezone
from django.utils.translation import gettext_lazy
from logging import getLogger
import numexpr as ne
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import re
import unicodedata

from .processing.get_structure_data import get_mol_structure_2dsvg, calc_ion_mass

logger = getLogger(__name__)


class BaseManager(models.Manager):
    def get_or_none(self, **kwargs):
        try:
            return self.get_queryset().get(**kwargs)
        except self.model.MultipleObjectsReturned as e:
            print(e)
        except self.model.DoesNotExist:
            return None


class CustomQuerySet(models.QuerySet):
    def update(self, **kwargs) -> int:
        logger.debug("Update 'updated_at' field.")
        kwargs['updated_at'] = timezone.now()
        return super().update(**kwargs)


def classification_default():
    return {'kingdom': '', 'superclass': '', 'class': '', 'subclass': '', 'direct_parent': ''}


class CompoundTemplate(models.Model):
    objects = BaseManager.from_queryset(CustomQuerySet)()

    class Meta:
        abstract = True

    cas_rn = models.CharField(blank=True, null=True, max_length=12, verbose_name='CAS Registry Number')
    compound_name = models.CharField(blank=True, null=True, max_length=255, verbose_name='Compound Name')
    synonyms = models.JSONField(default=list, blank=True, null=True, verbose_name='Synonyms')

    inchi = models.TextField(unique=True, blank=True, verbose_name='InChI')
    smiles = models.TextField(default="", blank=True, verbose_name='SMILES')

    # The following fields are automatically set to their values at save(); they have editable=False.
    exact_mass = models.FloatField(default=0, editable=False, verbose_name='Exact Mass')
    inchikey = models.CharField(unique=True, default="", max_length=27, editable=False, verbose_name='InChIKey')
    formula = models.CharField(default="", max_length=255, editable=False, verbose_name='Molecular Formula')
    icon_2dsvg = models.TextField(default="", editable=False, verbose_name='SVG Icon')
    structure_2dsvg = models.TextField(default="", editable=False, verbose_name='2D Chemical Structure')
    created_at = models.DateTimeField(blank=True, null=True, auto_now_add=True,
                                      editable=False, verbose_name='Creation Date')
    updated_at = models.DateTimeField(blank=True, null=True, auto_now=True, editable=False, verbose_name='Update Date')

    def __str_(self):
        return self.inchikey

    def clean(self):
        super().clean()

        # Start custom validation.
        if not (self.inchi or self.smiles):
            raise ValidationError(gettext_lazy('Either InChI or SMILES is required.'))
        else:
            mol = None
            if self.inchi:
                mol = Chem.MolFromInchi(self.inchi)
                if not mol:
                    raise ValidationError(gettext_lazy('RDKit cannot generate a Mol object.'))
            if self.smiles:
                mol = Chem.MolFromSmiles(self.smiles)
                if not mol:
                    raise ValidationError(gettext_lazy('RDKit cannot generate a Mol object.'))

            if self.inchi and self.smiles:
                mol_f_inchi = Chem.MolFromInchi(self.inchi)
                mol_f_smiles = Chem.MolFromSmiles(self.smiles)
                inchikey_f_inchi = Chem.MolToInchiKey(mol_f_inchi)
                inchikey_f_smiles = Chem.MolToInchiKey(mol_f_smiles)

                if inchikey_f_inchi != inchikey_f_smiles:
                    raise ValidationError(
                        gettext_lazy('The InChIKey generated from the InChI does not match that from the SMILES.'))
                # Finish custom validation.

            # Set values of the following fields.
            if self.inchi:
                mol = Chem.MolFromInchi(self.inchi)
            else:
                mol = Chem.MolFromSmiles(self.smiles)
                self.inchi = Chem.MolToInchi(mol)

            if not self.smiles:
                self.smiles = Chem.MolToSmiles(mol)

            self.exact_mass = rdMolDescriptors.CalcExactMolWt(mol)
            self.inchikey = Chem.MolToInchiKey(mol)
            self.formula = rdMolDescriptors.CalcMolFormula(mol)
            self.icon_2dsvg = get_mol_structure_2dsvg(mol, 200, 167)
            self.structure_2dsvg = get_mol_structure_2dsvg(mol, 400, 400)

    def save(self, *args, **kwargs):
        self.full_clean()
        super().save(*args, **kwargs)

    def calc_adduct_mass(self, adduct):
        equation = adduct.mass_equation.replace('M', str(self.exact_mass))
        adduct_mass = float(ne.evaluate(equation))

        return adduct_mass


class Compound(CompoundTemplate):
    objects = BaseManager.from_queryset(CustomQuerySet)()

    class Meta:
        verbose_name = 'Compound'
        verbose_name_plural = 'Compound'

    prefix = 'CMP'

    @classmethod
    def get_id_prefixed(cls, record_id):
        return f'{cls.prefix}{str(record_id).zfill(5)}'

    @property
    def id_prefixed(self):
        return f'{self.get_id_prefixed(self.id)}'

    def __str_(self):
        return self.cas_rn

    def clean(self):
        # Start custom validation.
        if self.cas_rn:
            cas_matching = re.fullmatch(r'\d{1,7}-\d{2}-\d', self.cas_rn)
            if not cas_matching:
                raise ValidationError(gettext_lazy('Invalid CAS RN.'))

        if self.compound_name:
            for char in self.compound_name:
                ret = unicodedata.east_asian_width(char)
                if ret != 'Na':
                    raise ValidationError(gettext_lazy(f'Invalid character in name: {char}'))

        super().clean()


class Classyfire(models.Model):
    objects = BaseManager.from_queryset(CustomQuerySet)()

    compound = models.OneToOneField('Compound', on_delete=models.CASCADE, blank=True, null=True,
                                    verbose_name='Compound ID', related_name='classyfire')

    kingdom = models.CharField(blank=True, null=True, max_length=20, verbose_name='Kingdom')
    superclass = models.CharField(blank=True, null=True, max_length=20, verbose_name='Superclass')
    chem_class = models.CharField(blank=True, null=True, max_length=20, verbose_name='Class')
    subclass = models.CharField(blank=True, null=True, max_length=20, verbose_name='Subclass')
    direct_parent = models.CharField(blank=True, null=True, max_length=20, verbose_name='Direct Parent')
