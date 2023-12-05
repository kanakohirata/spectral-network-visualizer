from django.core.files.storage import default_storage
from django.core.validators import FileExtensionValidator, ValidationError
from django import forms
from django.template import loader
from django.utils.safestring import mark_safe
from django.utils.translation import gettext_lazy
from logging import getLogger
import re


logger = getLogger(__name__)


def convert_form_errors_to_list(form, error_level='error'):
    error_list = []
    error_dict = form.errors.get_json_data()
    for field_name, messages in error_dict.items():
        logger.debug(field_name)
        field_id = form[field_name].auto_id
        error_list.append({'fieldId': field_id, 'fieldName': field_name, 'messages': messages,
                           'level': error_level})

    return error_list


class CustomParameterForm(forms.Form):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        for field in self.fields.values():
            field.required = False
            field.label_suffix = ''
            if field.widget.attrs.get('class'):
                field.widget.attrs['class'] += ' network-parameter'
            else:
                field.widget.attrs.update({'class': 'network-parameter'})


class CustomCheckboxSelectMultiple(forms.CheckboxSelectMultiple):
    template_name = 'visualizer/widgets/custom_checkbox.html'
    option_template_name = 'visualizer/widgets/custom_checkbox_option.html'

    def __init__(self, attrs=None):
        super().__init__(attrs)
        if 'class' in self.attrs:
            self.attrs['class'] += ' custom-checkbox'
        else:
            self.attrs['class'] = 'custom-checkbox'


class InputFileForm(CustomParameterForm):
    edge_file = forms.FileField(label='Edge File', widget=forms.FileInput())
    attribute_file = forms.FileField(label='Attribute File', widget=forms.FileInput())

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        for field in self.fields.values():
            field.required = False
            field.validators.append(FileExtensionValidator(['tsv']))


class BasicParameterForm(CustomParameterForm):
    score_threshold = forms.FloatField(label='Edge score threshold', min_value=0, max_value=1,
                                       widget=forms.NumberInput(attrs={'step': 0.0001, 'placeholder': 'default=0.75'}))

    subgraph_depth = forms.IntegerField(label='Subgraph depth', min_value=0,
                                        widget=forms.NumberInput(attrs={'step': 1, 'placeholder': 'default=3'}))

    mass_lower_limit = forms.FloatField(label='Mass lower threshold', min_value=0,
                                        widget=forms.NumberInput(attrs={'step': 0.0001, 'placeholder': 'default=100'}))

    mass_higher_limit = forms.FloatField(label='Mass higher threshold', min_value=0,
                                         widget=forms.NumberInput(attrs={'step': 0.0001, 'placeholder': 'default=1000'}))

    # community_detection_level is supposed to be an integer.
    # However, since only the greedy_modularity_communities method is currently used for community detection,
    # a Boolean value is instead used for that.
    # When the girvan_newman method becomes available in the future, an integer value should be obtained.
    community_detection_level = forms.BooleanField(label='Perform community detection')


class QuantitativeParameterForm(CustomParameterForm):
    feature_file = forms.FileField(label='Feature File', widget=forms.FileInput(),
                                   validators=[FileExtensionValidator(['csv', 'tsv'])])

    quant_polar_choices = (('', '--------------------'),
                           ('quant_val_both', 'Both'),
                           ('quant_val_increase', 'Increase'),
                           ('quant_val_decrease', 'Decrease'))

    quant_polar = forms.ChoiceField(label='Extraction type', choices=quant_polar_choices)

    subgraph_num_core_nodes = forms.IntegerField(label='Top N Node', min_value=1,
                                                 widget=forms.NumberInput(
                                                     attrs={'step': 1, 'placeholder': '1000'}
                                                 ))

    quant_value_type_choices = (('ratio', 'Ratio ( > 0 )'),
                                ('loading', 'Factor loading'))
    quant_value_type = forms.ChoiceField(label='Type of quantitative value', choices=quant_value_type_choices)

    quant_value = forms.FloatField(label='Threshold of quantitative value',
                                   widget=forms.NumberInput(attrs={'step': 0.0001, 'placeholder': '2'}))

    stat_value = forms.FloatField(label='P-value', min_value=0,
                                  widget=forms.NumberInput(attrs={'step': 0.0001, 'placeholder': '0.05'}))

    quant_subgraph_depth = forms.IntegerField(label='Extraction depth', min_value=1,
                                              widget=forms.NumberInput(attrs={'step': 1, 'placeholder': '10'}))


class GlobalAccessionParameterForm(CustomParameterForm):
    l_global_accession_for_node_select_subgraph = forms.CharField(
        label='Global Accession to Extract',
        help_text='Enter global accessions delimited with a semicolon.')

    node_select_subgraph_depth = forms.IntegerField(label='Extraction depth', min_value=1,
                                                    widget=forms.NumberInput(attrs={'step': 1, 'placeholder': '10'}))

    # node_select_subgraph_mode_choices = (('', '--------------------'),
    #                                      ('quant_val_both', 'Both'),
    #                                      ('quant_val_increase', 'Increase'),
    #                                      ('quant_val_decrease', 'Decrease'))
    # node_select_subgraph_mode = forms.ChoiceField(label='Extraction Mode', choices=node_select_subgraph_mode_choices)


class ReferenceLayerParameterForm(CustomParameterForm):
    filter_select_category_choices = (('none', 'Chemical Taxonomy (Superclass level)'),
                                      ('list_cmpd_classification_superclass', 'Chemical Taxonomy (Class level)'),
                                      ('list_compound_categories', 'Toxicity (T3DB)'))
    filter_select_category = forms.ChoiceField(label='Category', choices=filter_select_category_choices)

    filter_select_keyword_choices = (
        ('Alkaloids and derivatives', 'Alkaloids and derivatives'),
        ('Benzenoids', 'Benzenoids'),
        ('Lipids and lipid-like molecules', 'Lipids and lipid-like molecules'),
        ('Nucleosides, nucleotides, and analogues', 'Nucleosides, nucleotides, and analogues'),
        ('Organic acids and derivatives', 'Organic acids and derivatives'),
        ('Organic oxygen compounds', 'Organic oxygen compounds'),
        ('Organic Polymers', 'Organic Polymers'),
        ('Organoheterocyclic compounds', 'Organoheterocyclic compounds'),
        ('Phenylpropanoids and polyketides', 'Phenylpropanoids and polyketides')
    )
    filter_select_keyword = forms.ChoiceField(label='Filter', choices=filter_select_keyword_choices, disabled=True)

    color_toxic_compound = forms.BooleanField(label='Add color to T3DB compounds')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.fields['filter_select_keyword'].widget.attrs['id'] = 'chemical-superclass-filer'
        self.fields['filter_select_keyword'].widget.attrs['style'] = 'display: none'


class AnnotationParameterFrom(CustomParameterForm):
    mz_tolerance = forms.FloatField(
        label='<i>m/z</i> Tolerance (Da)', min_value=0,
        widget=forms.NumberInput(attrs={'step': 0.0001, 'placeholder': '0.01'})
    )

# ALL_ADDUCTS = ('[M+H]+', '[M+Na]+', '[M+3H]3+', '[M+2H+Na]3+', '[M+H+2Na]3+', '[M+3Na]3+', '[M+2H]2+',
#               '[M+H+NH4]2+', '[M+H+Na]2+', '[M+H+K]2+', '[M+ACN+2H]2+', '[M+2Na]2+', '[M+2ACN+2H]2+', '[M+3ACN+2H]2+',
#               '[M+NH4]+', '[M+CH3OH+H]+', '[M+K]+', '[M+ACN+H]+', '[M+2Na-H]+', '[M+IsoProp+H]+', '[M+ACN+Na]+',
#               '[M+2K-H]+', '[M+DMSO+H]+', '[M+2ACN+H]+', '[M+IsoProp+Na+H]+', '[2M+H]+', '[2M+NH4]+', '[2M+Na]+',
#               '[2M+K]+', '[2M+ACN+H]+', '[2M+ACN+Na]+', '[M-H]-', '[M-3H]3-', '[M-2H]2-', '[M-H2O-H]-', '[M+Na-2H]-',
#               '[M+Cl]-', '[M+K-2H]-', '[M+FA-H]-', '[M+Hac-H]-', '[M+Br]-', '[M+TFA-H]-', '[2M-H]-', '[2M+FA-H]-',
#               '[2M+Hac-H]-','[3M-H]-')


class SuspectCompoundParameterFrom(CustomParameterForm):
    suspect_file = forms.FileField(label='Suspect File',
                                   widget=forms.FileInput(attrs={'multiple': True, 'class': 'suspect-parameter'}),
                                   validators=[FileExtensionValidator(['tsv'])])

    ion_mode_choices = (('pos', 'Positive'), ('neg', 'Negative'))
    ion_mode = forms.ChoiceField(label='Ion Mode', choices=ion_mode_choices,
                                 widget=forms.Select(attrs={'class': 'suspect-parameter'}))

    adducts_for_suspect_compound_choices = (
        ('positive', (
            ('[M+H]+', '[M+H]+'),
            ('[M+Na]+', '[M+Na]+'),
            ('[M+3H]3+', '[M+3H]3+'),
            ('[M+2H+Na]3+', '[M+2H+Na]3+'),
            ('[M+H+2Na]3+', '[M+H+2Na]3+'),
            ('[M+3Na]3+', '[M+3Na]3+'),
            ('[M+2H]2+', '[M+2H]2+'),
            ('[M+H+NH4]2+', '[M+H+NH4]2+'),
            ('[M+H+Na]2+', '[M+H+Na]2+'),
            ('[M+H+K]2+', '[M+H+K]2+'),
            ('[M+ACN+2H]2+', '[M+ACN+2H]2+'),
            ('[M+2Na]2+', '[M+2Na]2+'),
            ('[M+2ACN+2H]2+', '[M+2ACN+2H]2+'),
            ('[M+3ACN+2H]2+', '[M+3ACN+2H]2+'),
            ('[M+NH4]+', '[M+NH4]+'),
            ('[M+CH3OH+H]+', '[M+CH3OH+H]+'),
            ('[M+K]+', '[M+K]+'),
            ('[M+ACN+H]+', '[M+ACN+H]+'),
            ('[M+2Na-H]+', '[M+2Na-H]+'),
            ('[M+IsoProp+H]+', '[M+IsoProp+H]+'),
            ('[M+ACN+Na]+', '[M+ACN+Na]+'),
            ('[M+2K-H]+', '[M+2K-H]+'),
            ('[M+DMSO+H]+', '[M+DMSO+H]+'),
            ('[M+2ACN+H]+', '[M+2ACN+H]+'),
            ('[M+IsoProp+Na+H]+', '[M+IsoProp+Na+H]+'),
            ('[2M+H]+', '[2M+H]+'),
            ('[2M+NH4]+', '[2M+NH4]+'),
            ('[2M+Na]+', '[2M+Na]+'),
            ('[2M+K]+', '[2M+K]+'),
            ('[2M+ACN+H]+', '[2M+ACN+H]+'),
            ('[2M+ACN+Na]+', '[2M+ACN+Na]+')
        )),
        ('negative', (
            ('[M-H]-', '[M-H]-'),
            ('[M-3H]3-', '[M-3H]3-'),
            ('[M-2H]2-', '[M-2H]2-'),
            ('[M-H2O-H]-', '[M-H2O-H]-'),
            ('[M+Na-2H]-', '[M+Na-2H]-'),
            ('[M+Cl]-', '[M+Cl]-'),
            ('[M+K-2H]-', '[M+K-2H]-'),
            ('[M+FA-H]-', '[M+FA-H]-'),
            ('[M+Hac-H]-', '[M+Hac-H]-'),
            ('[M+Br]-', '[M+Br]-'),
            ('[M+TFA-H]-', '[M+TFA-H]-'),
            ('[2M-H]-', '[2M-H]-'),
            ('[2M+FA-H]-', '[2M+FA-H]-'),
            ('[2M+Hac-H]-', '[2M+Hac-H]-'),
            ('[3M-H]-', '[3M-H]-')
        ))
    )
    adducts_for_suspect_compound = forms.MultipleChoiceField(
        choices=adducts_for_suspect_compound_choices,
        widget=CustomCheckboxSelectMultiple(attrs={'class': 'suspect-parameter'}),
    )


def validate_mass_defect(str_mass_defects):
    if str_mass_defects:
        try:
            str_mass_defects = str_mass_defects.rstrip(',')
            for mass_defect in str_mass_defects.split(','):
                mass_defect = mass_defect.strip()
                if mass_defect:
                    _ = float(mass_defect)
        except ValueError:
            raise ValidationError(
                gettext_lazy('Enter numbers delimited with a comma.')
            )


class MassDefectParameterForm(CustomParameterForm):
    str_mass_defects = forms.CharField(label='Mass Defects', validators=[validate_mass_defect],
                                       help_text='Enter numbers delimited with a comma.')


def validate_fragment_mz(str_fragment_mz_values: str):
    if str_fragment_mz_values:
        str_fragment_mz_values = str_fragment_mz_values.strip()
        if not re.search(r'\d+\.?\d+', str_fragment_mz_values):
            raise ValidationError(
                    gettext_lazy('Enter numbers delimited with a comma, space or newline.')
                )


class FragmentMzParameterForm(CustomParameterForm):
    mz_tolerance_for_fragment = forms.FloatField(
        label='<i>m/z</i> Tolerance (Da)', min_value=0,
        widget=forms.NumberInput(attrs={'step': 0.0001, 'placeholder': '0.01'})
    )

    str_fragment_mz_values = forms.CharField(label='Fragment m/z (OR Search)', validators=[validate_fragment_mz],
                                             widget=forms.Textarea(attrs={'rows': '5'}),
                                             help_text='Enter numbers delimited with a comma, space or newline.')
