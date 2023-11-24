import copy
from logging import getLogger
import networkx as nx
import numpy as np
import plotly.colors as pcolors
from rdkit import Chem
import sys
from visualizer.processing.multilayer_3d_ms_network_for_django import multilayer_3d_network
from visualizer.processing.multilayer_3d_ms_network_for_django import make_colormap
from visualizer.processing.get_structure_data import get_mol_structure_2dsvg_base64

logger = getLogger(__name__)


def get_text_color(rgb):
    brightness = max(rgb) / 255
    if brightness > 0.65:
        color = '#000000'
    else:
        color = '#ffffff'

    return color


def reconstruct_graph_data_v2(list_dic_edges_nodes_graph_by_layer,
                              list_of_edge_for_networkx_to_show_inter_sample_ref_layer,
                              list_of_edge_for_networkx_to_show_inter_sample_layer,
                              dic_layer_id_vs_attribute_for_layer,
                              l_dic_mesh3d_data):

    dic_attribute_for_layer_vs_layer_id = {}
    for layer_id, attribute_for_layer in dic_layer_id_vs_attribute_for_layer.items():
        dic_attribute_for_layer_vs_layer_id[attribute_for_layer] = layer_id

    logger.warning(f'dic_layer_id_vs_attribute_for_layer: {dic_layer_id_vs_attribute_for_layer}')
    logger.warning(f'dic_attribute_for_layer_vs_layer_id: {dic_attribute_for_layer_vs_layer_id}')

    list_layer_name_vs_node_trace = []
    list_layer_name_vs_edge_trace = []

    l_node_color_TOTAL = []
    layt_TOTAL = {}
    for dic in list_dic_edges_nodes_graph_by_layer:
        logger.debug(f'@@@ Keys: {dic.keys()}')

        dic_cluster_total_input_idx_MOD_vs_node_info = dic["dic_cluster_total_input_idx_MOD_vs_node_info"]

        f_base_layer = 0
        layer_name = dic["attribute_for_layer"]
        layer_id = dic_attribute_for_layer_vs_layer_id[layer_name]

        if layer_name in ["sample", "sample1", "sample2", "sample3", "base", "base1", "base2", "base3", "base4"]:
            f_base_layer = 1
            print('f_base_layer = 1')

        l_nodes = list(dic_cluster_total_input_idx_MOD_vs_node_info.keys())
        layt = dic["layout_3d"]

        # add layout of this layer to TOTAL layout
        for k, v in layt.items():
            layt_TOTAL[k] = v

        # get  edge color variation -------------
        # l_node_color = []
        # l_node_symbol = []
        # for k, node_info in dic["dic_cluster_total_input_idx_MOD_vs_node_info"].items():
        #     l_node_color.append(node_info.color)
        #     l_node_symbol.append(node_info.symbol)

        l_node_color = dic['l_node_color']
        l_node_symbol = dic['l_node_symbols']

        l_node_color_TOTAL.append(l_node_color)

        list_edge_nx = dic["list_of_edge_for_networkx"]

        l_nodes_str = []
        l_2d_svg_base64 = []
        for node in l_nodes:
            node_info = dic_cluster_total_input_idx_MOD_vs_node_info[node]
            #  2.7   node_info_str = "node:" + node.decode('utf-8') + ";\n" + node_info.spec_cluster.global_accession.decode('utf-8') + ";\n" + node_info.spec_cluster.compound_name.decode('utf-8')
            node_info_str = f'Node: {node}<br>Precursor m/z: {node_info["spec_cluster"]["represen_spec_uni"]["precursor_mz"]}<br>' \
                            f'Retention time (sec): {node_info["spec_cluster"]["represen_spec_uni"]["retention_time_in_sec"]}<br>' \
                            f'ID: {node_info["total_input_idx"]}<br>ID MOD: {node_info["total_input_idx_mod"]}<br>' \
                            f'Compound name: {node_info["spec_cluster"]["compound_name"]}<br>' \
                            f'Global accession: {node_info["spec_cluster"]["global_accession"]}<br>' \
                            f'Suspect name: {node_info["name"]}'

            inchi = node_info["spec_cluster"]["inchi"]
            svg_base64 = ''

            if inchi != '' and not isinstance(inchi, float):
                node_info_str += f'<br>InChI: {inchi}'
                mol = Chem.MolFromInchi(inchi)
                if mol:
                    svg_base64 = get_mol_structure_2dsvg_base64(mol, 200, 200)

            l_2d_svg_base64.append([svg_base64])

            if node_info["spec_cluster"]["inchi_key"]:
                node_info_str += f'<br>InChIKey: {node_info["spec_cluster"]["inchi_key"]}'

            l_nodes_str.append(node_info_str)

        # Xn = [layt[k][0] for k in l_nodes]
        # Yn = [layt[k][1] for k in l_nodes]
        # Zn = [layt[k][2] for k in l_nodes]

        Xn = dic['l_node_Xn_coordinates']
        Yn = dic['l_node_Yn_coordinates']
        Zn = dic['l_node_Zn_coordinates']

        trace = {
            'x': Xn, 'y': Yn, 'z': Zn, 'type': 'scatter3d', 'mode': 'markers', 'name': f'Compounds of {layer_name}',
            'marker': {'symbol': l_node_symbol, 'size': 6, 'color': l_node_color, 'coloraxis': "coloraxis",
                       'line': {'color': 'rgb(50,50,50)', 'width': 0.5}},
            'text': l_nodes_str, 'hoverinfo': 'text', 'hoverlabel': {'align': 'left'},
            'legendgroup': layer_id,
            'customdata': l_2d_svg_base64
        }

        list_layer_name_vs_node_trace.append([layer_name, trace])

        l_edge_Xe_coordinates = []
        l_edge_Ye_coordinates = []
        l_edge_Ze_coordinates = []
        l_edge_color = []
        l_edge_str = []
        edge_custom_data = []

        for n in range(len(list_edge_nx)):
            l_edge_Xe_coordinates += [layt[list_edge_nx[n][0]][0], layt[list_edge_nx[n][1]][0], None]
            l_edge_Ye_coordinates += [layt[list_edge_nx[n][0]][1], layt[list_edge_nx[n][1]][1], None]
            l_edge_Ze_coordinates += [layt[list_edge_nx[n][0]][2], layt[list_edge_nx[n][1]][2], None]

            l_edge_color.append(list_edge_nx[n][2]["color"])
            l_edge_color.append(list_edge_nx[n][2]["color"])
            l_edge_color.append(list_edge_nx[n][2]["color"])

            edge_info_str = f'Score: {list_edge_nx[n][2]["spec_sim_score"]}'
            l_edge_str.extend([edge_info_str, edge_info_str, None])
            edge_custom_data.extend([
                [list_edge_nx[n][2]["color"], ],
                [list_edge_nx[n][2]["color"], ],
                [None, ]
            ])

        trace = {
            'x': l_edge_Xe_coordinates, 'y': l_edge_Ye_coordinates, 'z': l_edge_Ze_coordinates,
            'type': 'scatter3d', 'mode': 'lines', 'name': f'Similarities between spectra of {layer_name}',
            'line': {'color': l_edge_color, 'width': 2},
            'text': l_edge_str, 'hoverinfo': 'text',
            'legendgroup': layer_id,
            'customdata': edge_custom_data
        }

        list_layer_name_vs_edge_trace.append([layer_name, trace])

    color_max = 0
    color_min = 0

    for l_node_color in l_node_color_TOTAL:
        # logger.warning(f'l_node_color: {l_node_color}')
        l_node_color_num = [node_color for node_color in l_node_color
                            if not isinstance(node_color, str)]
        if not l_node_color_num:
            continue

        _color_max = max(l_node_color_num)
        _color_min = min(l_node_color_num)

        if _color_max > color_max:
            color_max = _color_max

        if _color_min < color_min:
            color_min = _color_min

    if color_min == 0 and color_max == 0:
        color_max = 0.5
        color_min = -0.5

    # logger.warning(f'color_max: {color_max}, color_min: {color_min}')
    cmap = make_colormap.from_named_colorscale('Hot')
    for layer_name, trace in list_layer_name_vs_node_trace:
        l_idx_vs_node_color_str = []
        l_idx_vs_node_color_num = []

        for idx, value in enumerate(trace['marker']['color']):
            if isinstance(value, str):
                l_idx_vs_node_color_str.append((idx, value))
            else:
                l_idx_vs_node_color_num.append((idx, value))

        l_text_color = [''] * len(trace['marker']['color'])
        l_node_rgb_225_str = [''] * len(trace['marker']['color'])
        if l_idx_vs_node_color_num:
            l_scaled_color = [make_colormap.scale_0to1(value, color_max, color_min) for idx, value in l_idx_vs_node_color_num]
            # logger.warning(f'l_scaled_color: {l_scaled_color}')
            l_node_rgb_0to1 = list(map(cmap, l_scaled_color))
            l_node_rgb_225 = list(map(pcolors.convert_to_RGB_255, l_node_rgb_0to1))

            l_text_color_temp = list(map(get_text_color, l_node_rgb_225))
            l_node_rgb_225_str_temp = list(map(pcolors.label_rgb, l_node_rgb_225))

            for i, (text_color, node_rgb_225_str) in enumerate(zip(l_text_color_temp, l_node_rgb_225_str_temp)):
                idx = l_idx_vs_node_color_num[i][0]
                l_text_color[idx] = text_color
                l_node_rgb_225_str[idx] = node_rgb_225_str

        if l_idx_vs_node_color_str:
            l_node_rgb_225 = list(map(lambda xcolor: pcolors.hex_to_rgb(xcolor[1]), l_idx_vs_node_color_str))
            l_text_color_temp = list(map(get_text_color, l_node_rgb_225))
            l_node_rgb_225_str_temp = list(map(pcolors.label_rgb, l_node_rgb_225))

            for i, (text_color, node_rgb_225_str) in enumerate(zip(l_text_color_temp, l_node_rgb_225_str_temp)):
                idx = l_idx_vs_node_color_str[i][0]
                l_text_color[idx] = text_color
                l_node_rgb_225_str[idx] = node_rgb_225_str

        arr_node_rgb_225_str = np.array([[node_rgb_225_str, text_color] for node_rgb_225_str, text_color
                                         in zip(l_node_rgb_225_str, l_text_color)])
        arr_current_customdata = np.array(trace['customdata'])
        # logger.warning(f'arr_node_rgb_225_str: {layer_name}\n{arr_node_rgb_225_str}')
        customdata = np.concatenate((arr_node_rgb_225_str, arr_current_customdata), axis=1).tolist()
        trace['customdata'] = customdata

        if trace['legendgroup'].startswith('r'):
            trace['marker']['color'] = l_node_rgb_225_str  # TODO: Update

        ### INTER SAMPLE-REF EDGES --------------
    l_inter_sample_ref_e = list_of_edge_for_networkx_to_show_inter_sample_ref_layer
    l_INTER_SAMPLE_REF_LAYER_edge_Xe_coordinates = []
    l_INTER_SAMPLE_REF_LAYER_edge_Ye_coordinates = []
    l_INTER_SAMPLE_REF_LAYER_edge_Ze_coordinates = []

    l_INTER_SAMPLE_REF_LAYER_edge_color = []
    l_INTER_SAMPLE_REF_LAYER_custom_data = []

    # Split edge traces between sample and references by reference layer.
    d_layer_name_vs_edge_trace_between_sample_and_ref = {}

    for n in range(len(l_inter_sample_ref_e)):
        # note   networkx "layout"  is a dictionary where key is node id and value is array of geometry of node.
        Xe_coordinates = [layt_TOTAL[l_inter_sample_ref_e[n][0]][0],
                          layt_TOTAL[l_inter_sample_ref_e[n][1]][0], None]

        Ye_coordinates = [layt_TOTAL[l_inter_sample_ref_e[n][0]][1],
                          layt_TOTAL[l_inter_sample_ref_e[n][1]][1], None]

        Ze_coordinates = [layt_TOTAL[l_inter_sample_ref_e[n][0]][2],
                          layt_TOTAL[l_inter_sample_ref_e[n][1]][2], None]

        color = l_inter_sample_ref_e[n][2]["color"]

        for x, y, z in zip(Xe_coordinates, Ye_coordinates, Ze_coordinates):
            layer_id = ''
            if z is None or z < 0:
                continue
            else:
                for dic_mesh3d_data in l_dic_mesh3d_data:
                    x_min = min(dic_mesh3d_data['l_x']) - 0.1
                    x_max = max(dic_mesh3d_data['l_x']) + 0.1
                    y_min = min(dic_mesh3d_data['l_y']) - 0.1
                    y_max = max(dic_mesh3d_data['l_y']) + 0.1

                    z_layer = dic_mesh3d_data['l_z'][0]

                    if z == z_layer and (x_min < x < x_max) and (y_min < y < y_max):
                        layer_id = dic_mesh3d_data['text'].split('layer_id_')[1]
                        attribute_for_layer = dic_layer_id_vs_attribute_for_layer[layer_id]
                        layer_name = f'Edges between sample and {attribute_for_layer}'
                        break

            if not layer_id:
                logger.warning(f'Invalid edge: This edge does not belong to any reference layer.\n'
                               f'x: {x}, y: {y}, z: {z}')
                continue
            else:
                if layer_name in d_layer_name_vs_edge_trace_between_sample_and_ref:
                    trace = d_layer_name_vs_edge_trace_between_sample_and_ref[layer_name]
                    trace['x'] += Xe_coordinates
                    trace['y'] += Ye_coordinates
                    trace['z'] += Ze_coordinates
                    trace['line']['color'].extend([color, color, color])
                    trace['customdata'].append(color)
                else:
                    trace = {
                        'x': Xe_coordinates,
                        'y': Ye_coordinates,
                        'z': Ze_coordinates,
                        'type': 'scatter3d', 'mode': 'lines', 'name': layer_name,
                        'line': {'color': [color, color, color], 'width': 1}, 'hoverinfo': 'skip',
                        'customdata': [color,], 'legendgroup': layer_id
                    }
                    d_layer_name_vs_edge_trace_between_sample_and_ref[layer_name] = trace

        l_INTER_SAMPLE_REF_LAYER_edge_Xe_coordinates += Xe_coordinates
        l_INTER_SAMPLE_REF_LAYER_edge_Ye_coordinates += Ye_coordinates
        l_INTER_SAMPLE_REF_LAYER_edge_Ze_coordinates += Ze_coordinates
        l_INTER_SAMPLE_REF_LAYER_edge_color.append(l_inter_sample_ref_e[n][2]["color"])
        l_INTER_SAMPLE_REF_LAYER_edge_color.append(l_inter_sample_ref_e[n][2]["color"])
        l_INTER_SAMPLE_REF_LAYER_edge_color.append(l_inter_sample_ref_e[n][2]["color"])

        l_INTER_SAMPLE_REF_LAYER_custom_data.append([l_inter_sample_ref_e[n][2]["color"]])

    for layer_name, edge_trace in d_layer_name_vs_edge_trace_between_sample_and_ref.items():
        list_layer_name_vs_edge_trace.append((layer_name, edge_trace))

    edge_trace_between_sample_and_ref = {
        'x': l_INTER_SAMPLE_REF_LAYER_edge_Xe_coordinates,
        'y': l_INTER_SAMPLE_REF_LAYER_edge_Ye_coordinates,
        'z': l_INTER_SAMPLE_REF_LAYER_edge_Ze_coordinates,
        'type': 'scatter3d', 'mode': 'lines', 'name': 'INTER SAMPLE REF LAYER EDGES',
        'line': {'color': l_INTER_SAMPLE_REF_LAYER_edge_color, 'width': 1}, 'hoverinfo': 'skip',
        'customdata': l_INTER_SAMPLE_REF_LAYER_custom_data,
    }

    # inter sample layers ------------------
    l_inter_sample_e = list_of_edge_for_networkx_to_show_inter_sample_layer
    l_INTER_SAMPLE_LAYER_edge_Xe_coordinates = []
    l_INTER_SAMPLE_LAYER_edge_Ye_coordinates = []
    l_INTER_SAMPLE_LAYER_edge_Ze_coordinates = []

    l_INTER_SAMPLE_LAYER_edge_color = []

    for n in range(len(l_inter_sample_e)):
        # note   networkx "layout"  is a dictionary where key is node id and value is array of geometry of node.

        Xe_coordinates = [layt_TOTAL[l_inter_sample_e[n][0]][0],
                          layt_TOTAL[l_inter_sample_e[n][1]][0], None]

        Ye_coordinates = [layt_TOTAL[l_inter_sample_e[n][0]][1],
                          layt_TOTAL[l_inter_sample_e[n][1]][1], None]

        Ze_coordinates = [layt_TOTAL[l_inter_sample_e[n][0]][2],
                          layt_TOTAL[l_inter_sample_e[n][1]][2], None]

        l_INTER_SAMPLE_LAYER_edge_Xe_coordinates += Xe_coordinates
        l_INTER_SAMPLE_LAYER_edge_Ye_coordinates += Ye_coordinates
        l_INTER_SAMPLE_LAYER_edge_Ze_coordinates += Ze_coordinates
        l_INTER_SAMPLE_LAYER_edge_color.append(l_inter_sample_e[n][2]["color"])
        l_INTER_SAMPLE_LAYER_edge_color.append(l_inter_sample_e[n][2]["color"])
        l_INTER_SAMPLE_LAYER_edge_color.append(l_inter_sample_e[n][2]["color"])

        l_INTER_SAMPLE_LAYER_edge_color.append([l_inter_sample_e[n][2]["color"]])

    edge_trace_inter_sample = {
        'x': l_INTER_SAMPLE_LAYER_edge_Xe_coordinates,
        'y': l_INTER_SAMPLE_LAYER_edge_Ye_coordinates,
        'z': l_INTER_SAMPLE_LAYER_edge_Ze_coordinates,
        'type': 'scatter3d', 'mode': 'lines', 'name': 'INTER SAMPLE LAYER EDGES',
        'line': {'color': l_INTER_SAMPLE_LAYER_edge_color, 'width': 3}, 'hoverinfo': 'skip',
        'customdata': l_INTER_SAMPLE_LAYER_edge_color
    }

    #######################################################################################################
    #######################################################################################################

    list_layer_id_vs_layer_trace = []
    for dic_mesh3d_data in l_dic_mesh3d_data:
        logger.warning(f'dic_mesh3d_data: {dic_mesh3d_data}')
        layer_id = dic_mesh3d_data['text'].split('layer_id_')[1]
        logger.warning(f"dic_mesh3d_data['text']: {dic_mesh3d_data['text']}, layer_id: {layer_id}")

        trace = {'x': dic_mesh3d_data['l_x'], 'y': dic_mesh3d_data['l_y'], 'z': dic_mesh3d_data['l_z'],
                 'type': 'mesh3d', 'text': f'layer_id_{layer_id}',
                 'color': dic_mesh3d_data['color'], 'opacity': 0.20, 'hoverinfo': 'skip',
                 'legendgroup': layer_id}

        list_layer_id_vs_layer_trace.append([layer_id, trace])

    return (list_layer_name_vs_node_trace, list_layer_name_vs_edge_trace,
            edge_trace_between_sample_and_ref, edge_trace_inter_sample, list_layer_id_vs_layer_trace)
