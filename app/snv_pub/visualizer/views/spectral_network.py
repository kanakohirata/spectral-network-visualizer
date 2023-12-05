from django.conf import settings
from django.core.cache import cache
from django.http import JsonResponse
from django.utils.decorators import method_decorator
from django.views import generic
from django.views.decorators.csrf import csrf_exempt, csrf_protect
from django.shortcuts import render
from logging import getLogger
import os
import pandas as pd
import re
import time
from ..files.uploadhandler import TemporaryDirectoryUploadHandler
from ..forms import *
from ..processing.multilayer_3d_ms_network import config_ml3dnet
from ..processing.multilayer_3d_ms_network import multilayer_3d_network
from ..processing.multilayer_3d_ms_network import reconstruct_graph_data
from ..processing.multilayer_3d_ms_network import make_colormap

logger = getLogger(__name__)
BASE_DIR = settings.BASE_DIR


class SpectralNetworkView(generic.TemplateView):
    template_name = 'visualizer/spectral_network.html'
    input_file_form_class = InputFileForm
    basic_parameter_form_class = BasicParameterForm
    quantitative_parameter_form_class = QuantitativeParameterForm
    global_accession_parameter_form_class = GlobalAccessionParameterForm
    reference_layer_parameter_form_class = ReferenceLayerParameterForm
    annotation_parameter_form_class = AnnotationParameterFrom
    suspect_compound_parameter_form_class = SuspectCompoundParameterFrom
    mass_defect_parameter_form_class = MassDefectParameterForm
    fragment_mz_parameter_from_class = FragmentMzParameterForm

    def get(self, request, *args, **kwargs):
        logger.debug(request.META.get('REMOTE_ADDR'))
        return super().get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        adduct_dict = cache.get('adduct_dict')
        if not adduct_dict:
            logger.debug("'adduct_dict' is not in cache. Read adduct_calculation.tsv")
            adduct_filepath = 'visualizer/data/adduct_calculation.tsv'
            adduct_tabel = pd.read_table(adduct_filepath, index_col='Ion_Name')
            adduct_dict = adduct_tabel.to_dict(orient='index')
            cache.set('adduct_dict', adduct_dict, 60 * 10)  # 10 min

        logger.debug(adduct_dict)

        adduct_names_pos = []
        adduct_names_neg = []
        for adduct_name, adduct_info in adduct_dict.items():
            if adduct_info['Ion_Mode'] == 'positive':
                adduct_names_pos.append(adduct_name)
            else:
                adduct_names_neg.append(adduct_name)

        context['adduct_names_pos'] = adduct_names_pos
        context['adduct_names_neg'] = adduct_names_neg

        if not kwargs.get('input_file_form'):
            context['input_file_form'] = self.input_file_form_class()
        if not kwargs.get('basic_parameter_form'):
            context['basic_parameter_form'] = self.basic_parameter_form_class()
        if not kwargs.get('quantitative_parameter_form'):
            context['quantitative_parameter_form'] = self.quantitative_parameter_form_class()
        if not kwargs.get('global_accession_parameter_form'):
            context['global_accession_parameter_form'] = self.global_accession_parameter_form_class()
        if not kwargs.get('reference_layer_parameter_form'):
            context['reference_layer_parameter_form'] = self.reference_layer_parameter_form_class()
        if not kwargs.get('annotation_parameter_form'):
            context['annotation_parameter_form'] = self.annotation_parameter_form_class()
        if not kwargs.get('suspect_compound_parameter_form'):
            context['suspect_compound_parameter_form'] = self.suspect_compound_parameter_form_class()
        if not kwargs.get('mass_defect_parameter_form'):
            context['mass_defect_parameter_form'] = self.mass_defect_parameter_form_class()
        if not kwargs.get('fragment_mz_parameter_from'):
            context['fragment_mz_parameter_from'] = self.fragment_mz_parameter_from_class()

        return context


@method_decorator(csrf_exempt, name='dispatch')
class GetNetworkData(generic.TemplateView):
    template_name = 'visualizer/spectral_network.html'
    all_form_classes = (InputFileForm, BasicParameterForm, QuantitativeParameterForm, GlobalAccessionParameterForm,
                        ReferenceLayerParameterForm, AnnotationParameterFrom, SuspectCompoundParameterFrom,
                        MassDefectParameterForm, FragmentMzParameterForm)

    def setup(self, request, *args, **kwargs):
        request.upload_handlers = [TemporaryDirectoryUploadHandler(request)]
        super().setup(request, *args, **kwargs)

    @method_decorator(csrf_protect)
    def post(self, request):
        logger.debug(request.POST.getlist('adducts_for_suspect_compound'))
        all_forms_is_valid = True
        field_errors = []
        for form_class in self.all_form_classes:
            form = form_class(request.POST, request.FILES or None)
            if not form.is_valid():
                all_forms_is_valid = False
                field_errors.extend(convert_form_errors_to_list(form))

        if all_forms_is_valid:
            pass
        else:
            logger.debug('Form is invalid!!!!!!!!!!!!!!!!!!!!!!!')
            for field_error in field_errors:
                logger.debug(field_error)
            return JsonResponse({'formValidation': {'isValid': False},
                                 'fieldErrors': field_errors})

        adduct_dict = cache.get('adduct_dict')
        if not adduct_dict:
            logger.debug("'adduct_dict' is not in cache. Read adduct_calculation.tsv")
            adduct_filepath = 'visualizer/data/adduct_calculation.tsv'
            adduct_tabel = pd.read_table(adduct_filepath, index_col='Ion_Name')
            adduct_dict = adduct_tabel.to_dict(orient='index')
            cache.set('adduct_dict', adduct_dict, 60 * 10)  # 10 min

        time_begin = time.time()
        logger.debug(request)
        logger.debug(request.FILES)
        edge_file = request.FILES.get('edge_file')
        attribute_file = request.FILES.get('attribute_file')
        feature_file = request.FILES.get('feature_file')
        compound_files = request.FILES.getlist('compound_file')
        suspect_files = request.FILES.getlist('suspect_file')
        create_or_update = request.POST.get("create_or_update")
        logger.debug(f'Create / Update: {create_or_update}')
        logger.debug(f'Suspect files: {suspect_files}')

        if not edge_file:
            edge_filepath = os.path.join(
                BASE_DIR,
                'visualizer/processing/multilayer_3d_ms_network/output_test_0/output_.edgeinfo.tsv'
            )
        else:
            edge_filepath = edge_file.temporary_file_path()

        if not attribute_file:
            attribute_filepath = os.path.join(
                BASE_DIR,
                'visualizer/processing/multilayer_3d_ms_network/output_test_0/output_.cluster_attribute.tsv'
            )
        else:
            attribute_filepath = attribute_file.temporary_file_path()

        if feature_file:
            feature_filepath = feature_file.temporary_file_path()
        else:
            feature_filepath = ''

        list_compound_filepath = []
        if compound_files:
            for compound_file in compound_files:
                list_compound_filepath.append(compound_file.temporary_file_path())

        list_suspect_filepath = []
        if suspect_files:
            for suspect_file in suspect_files:
                list_suspect_filepath.append(suspect_file.temporary_file_path())

        logger.debug(f'Path of attribute_file: {attribute_filepath}')

        logger.debug(f'request.POST: {request.POST}')
        logger.debug(f'request.session: {request.session}')
        logger.debug(f'request.FILES: {request.FILES}')

        # Create a new config dictionary and set parameters.
        config = config_ml3dnet.get_dic_config_initialized()
        config['filename_edge_info'] = edge_filepath
        config['filename_cluster_info'] = attribute_filepath
        config['filename_feature_table'] = feature_filepath
        config['list_filename_suspect_mz'] = list_suspect_filepath

        # Basic parameters
        if request.POST.get('score_threshold'):
            config['score_threshold'] = float(request.POST.get('score_threshold'))
        if request.POST.get('subgraph_depth'):
            config['subgraph_depth'] = int(request.POST.get('subgraph_depth'))
        if request.POST.get('mass_lower_limit'):
            config['mass_lower_limit'] = float(request.POST.get('mass_lower_limit'))
        if request.POST.get('mass_higher_limit'):
            config['mass_higher_limit'] = float(request.POST.get('mass_higher_limit'))
        if request.POST.get('community_detection_level') == 'true':
            config['n_level_community_detection'] = 1
        else:
            config['n_level_community_detection'] = 0
        if request.POST.get('filter_select_category'):
            config['filter_select_category'] = request.POST.get('filter_select_category')
        if request.POST.get('filter_select_keyword'):
            config['filter_select_keyword'] = request.POST.get('filter_select_keyword')

        if config.get('filter_select_category') == 'none':
            config['type_attribute_for_layer_separation'] = 'list_cmpd_classification_superclass'
        elif config.get('filter_select_category') == 'list_cmpd_classification_superclass':
            config['type_attribute_for_layer_separation'] = 'list_cmpd_classification_class'
        elif config.get('filter_select_category') == 'list_compound_categories':  # TODO: Temporary
            config['type_attribute_for_layer_separation'] = 'list_compound_categories'
            config['foldername_ext_cmpd_info'] = os.path.join(BASE_DIR,
                                                              'visualizer/processing/multilayer_3d_ms_network/t3db_xml/')  # TODO: Temporary
            logger.warning(f'config["foldername_ext_cmpd_info"]: {config["foldername_ext_cmpd_info"]}')

        if request.POST.get('color_toxic_compound') == 'true':
            config['color_toxic_compound'] = True
            config['foldername_ext_cmpd_info'] = os.path.join(BASE_DIR,
                                                              'visualizer/processing/multilayer_3d_ms_network/t3db_xml/')  # TODO: Temporary

        # Quantitative parameters
        if feature_file:
            config['quant_polar'] = request.POST.get('quant_polar')

            if request.POST.get('subgraph_num_core_nodes'):
                config['subgraph_num_core_nodes'] = int(request.POST.get('subgraph_num_core_nodes'))
            if request.POST.get('quant_value_type'):
                config['quant_value_type'] = request.POST.get('quant_value_type')
            if request.POST.get('quant_value'):
                config['quant_value'] = float(request.POST.get('quant_value'))
            if request.POST.get('stat_value'):
                config['stat_value'] = float(request.POST.get('stat_value'))
            if request.POST.get('quant_subgraph_depth'):
                config['subgraph_depth'] = int(request.POST.get('quant_subgraph_depth'))

        # Parameters for extraction by global accession
        str_global_accessions_for_node_select_subgraph = request.POST.get('l_global_accession_for_node_select_subgraph',
                                                                          '')
        if len(str_global_accessions_for_node_select_subgraph) >= 3:
            l_global_accession_for_node_select_subgraph = \
                [global_accession.strip() for global_accession in
                 str_global_accessions_for_node_select_subgraph.split(';')]

            config['l_global_accession_for_node_select_subgraph'] = l_global_accession_for_node_select_subgraph

            if request.POST.get('node_select_subgraph_depth'):
                config['node_select_subgraph_depth'] = int(request.POST.get('node_select_subgraph_depth'))

            config['node_select_subgraph_mode'] = request.POST.get('node_select_subgraph_mode')

        # Parameters for suspect compound
        ion_mode = request.POST.get('ion_mode')
        list_adduct_for_suspect = request.POST.getlist('adducts_for_suspect_compound', [])
        config['l_adduct_type_for_suspect'] = list_adduct_for_suspect
        for adduct_name in list_adduct_for_suspect:
            logger.debug(f"{adduct_name}: {adduct_dict[adduct_name]['Adduct_Mass']}")
            config['l_adduct_mass_for_suspect'].append(adduct_dict[adduct_name]['Adduct_Mass'])

        logger.debug(list_adduct_for_suspect)
        logger.debug(f"l_adduct_mass_for_suspect: {config['l_adduct_mass_for_suspect']}")

        str_mass_defects = request.POST.get('str_mass_defects', '')
        list_mass_defect = []
        if str_mass_defects:
            str_mass_defects = str_mass_defects.rstrip(',')
            for mass_defect in str_mass_defects.split(','):
                mass_defect = mass_defect.strip()
                if mass_defect:
                    list_mass_defect.append(float(mass_defect))
        config['list_mass_defect'] = list_mass_defect
        logger.debug(f"list_mass_defect: {config['list_mass_defect']}")

        if config['l_adduct_mass_for_suspect'] or config['list_mass_defect']:
            if request.POST.get('mz_tolerance'):
                config['mz_tol'] = float(request.POST.get('mz_tolerance'))

        # Parameters for fragment
        str_fragment_mz_values = request.POST.get('str_fragment_mz_values', '')
        list_fragment_mz = []
        if re.search(r'\d+\.?\d+', str_fragment_mz_values):
            list_fragment_mz = [float(n) for n in re.findall(r'\d+\.?\d+', str_fragment_mz_values)]
        config['list_product_mz_required'] = list_fragment_mz

        if config['list_product_mz_required'] and request.POST.get('mz_tolerance_for_fragment'):
            config['mz_tolerance_for_fragment'] = float(request.POST.get('mz_tolerance_for_fragment'))

        # List of total_input_idx_MOD to remove
        logger.debug(request.POST.get('l_total_input_idx_to_remove'))
        logger.debug(type(request.POST.get('l_total_input_idx_to_remove')))
        if request.POST.getlist('l_total_input_idx_to_remove'):
            config['l_total_input_idx_to_remove'] = request.POST.getlist('l_total_input_idx_to_remove')
        # logger.debug(f"l_total_input_idx_to_remove: {config['l_total_input_idx_to_remove']}")

        # Create or Update
        dic_source_data = None
        if create_or_update == 'update':
            logger.debug('Use session data.')
            dic_source_data = request.session['dic_source_data']

            if not dic_source_data:
                logger.debug('Session data is empty.')
                create_or_update = 'create'

        if create_or_update == 'create':
            dic_source_data = multilayer_3d_network.read_data_for_multilayer_3d_network(config)
            request.session['dic_source_data'] = dic_source_data

        if not dic_source_data:
            dic_source_data = multilayer_3d_network.read_data_for_multilayer_3d_network(config)
            request.session['dic_source_data'] = dic_source_data

        time_before_process_3d_network_data = time.time()
        dic_processed_data = multilayer_3d_network.process_3d_network_data(dic_source_data, config)
        time_after_process_3d_network_data = time.time()
        dic_dataset_for_3d_network_visualization = \
            multilayer_3d_network.create_data_for_3d_visualization(dic_processed_data, config)
        time_after_create_data_for_3d_visualization = time.time()

        # temp = dic_dataset_for_3d_network_visualization['list_dic_edges_nodes_graph_by_layer'][0]
        # logger.warning(f'temp: {temp}')

        (list_layer_name_vs_node_trace,
         list_layer_name_vs_edge_trace,
         edge_trace_between_sample_and_ref,
         edge_trace_inter_sample,
         list_layer_id_vs_layer_trace) \
            = reconstruct_graph_data(
            dic_dataset_for_3d_network_visualization['list_dic_edges_nodes_graph_by_layer'],
            dic_processed_data['list_of_edge_for_networkx_to_show_inter_sample_ref_layer'],
            dic_processed_data['list_of_edge_for_networkx_to_show_inter_sample_layer'],
            # dic_processed_data['dic_layer_id_vs_attribute_for_layer'],  # TODO: Changed to dic_dataset_for_3d_network_visualization['dic_layer_id_vs_attribute_for_layer']
            dic_dataset_for_3d_network_visualization['dic_layer_id_vs_attribute_for_layer'],
            dic_dataset_for_3d_network_visualization['l_dic_mesh3d_data']
        )
        time_after_reconstruct_graph_data_v2 = time.time()

        list_layer_name_vs_trace = list_layer_id_vs_layer_trace \
                                   + list_layer_name_vs_edge_trace \
                                   + list_layer_name_vs_node_trace

        colorscale = make_colormap.split_colorscale('Hot')

        logger.debug(f"l_dic_annotations: {dic_dataset_for_3d_network_visualization['l_dic_annotations']}")

        l_dic_annotations = dic_dataset_for_3d_network_visualization['l_dic_annotations']
        for dic_annotatoin in l_dic_annotations:
            dic_annotatoin['font'] = {'color': 'blue'}
            dic_annotatoin['arrowcolor'] = 'blue'

        layout_updatemenus = [
            {
                'active': 1, 'x': 1, 'xanchor': 'right', 'y': 1, 'yanchor': 'top',
                'buttons': [
                    {'label': 'Display layer attribute', 'method': 'relayout',
                     'args': ['scene.annotations', l_dic_annotations],
                     'font': {'color': 'blue'}},
                    {'label': 'Hide layer attribute', 'method': 'relayout',
                     'args': ['scene.annotations', []]}],
            },
        ]

        # layout = {'coloraxis': {'colorscale': hot, 'showscale': True, 'colorbar': {'x': 3}}}
        layout = {'coloraxis': {'colorscale': colorscale, 'colorbar': {'orientation': 'h'}},
                  'legend': {'groupclick': 'toggleitem'}, 'updatemenus': layout_updatemenus, }

        logger.debug(f'All processing time: {time.time() - time_begin} s')

        logger.debug(f'Time of process_3d_network_data: '
                     f'{time_after_process_3d_network_data - time_before_process_3d_network_data} s')

        logger.debug(f'Time of create_data_for_3d_visualization: '
                     f'{time_after_create_data_for_3d_visualization - time_after_process_3d_network_data} s')

        logger.debug(f'Time of reconstruct_graph_data_v2: '
                     f'{time_after_reconstruct_graph_data_v2 - time_after_create_data_for_3d_visualization} s')

        for layer_id, layer_trace in list_layer_id_vs_layer_trace:
            logger.warning(f'layer_id: {layer_id}\n{layer_trace}')
        
        return JsonResponse({'formValidation': {'isValid': True},
                             'listLayerNameVsTrace': list_layer_name_vs_trace,
                             # 'edgeTraceBetweenSampleAndRef': edge_trace_between_sample_and_ref,
                             'edgeTraceInterSample': edge_trace_inter_sample,
                             'listMeshAnnotation': dic_dataset_for_3d_network_visualization['l_dic_annotations'],
                             'layout': layout})
