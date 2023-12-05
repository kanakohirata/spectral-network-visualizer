from logging import getLogger


logger = getLogger(__name__)


def define_layers(key_attribute_to_base_layer: str,
                  dic_cluster_total_input_idx_MOD_vs_node_info: dict):
    logger.debug('starting [define_layers]')

    ###################################################
    # iterate list of node and add node with attribute
    ###################################################
    list_attribute_for_layer = []
    for total_input_idx_mod, node_info in dic_cluster_total_input_idx_MOD_vs_node_info.items():
        list_attribute_for_layer.append(node_info['layer'])

    list_attribute_for_layer = list(set(list_attribute_for_layer))

    logger.debug(f"list_attribute_for_layer {list_attribute_for_layer}")

    ######################################################
    # [V]Now  list_attribute_for_layer  contains all attribute for LAYER SEPARATION.
    #  define layer based on the attributes

    # but first you want to define base layer
    # Tell which attribute is used to specify base layer.
    #  if  your base layer spectra has not CMPD_CLASSIFICATION_SUPERCLASS description, use    str_key_attribute_to_base_layer = '[]'

    # you want to remove base layer key from list.

    if key_attribute_to_base_layer in list_attribute_for_layer:
        list_attribute_for_layer.remove(key_attribute_to_base_layer)

    # create dictionary where key is str of attribute ['organic compound'] and value is layer_id
    # !!!!!!!!!!!!!!!!   here you are making "sample" layer....
    dic_attribute_for_layer_vs_layer_id = {key_attribute_to_base_layer: 0}

    for n in range(0, len(list_attribute_for_layer)):
        # it needs to be n+1 since 0 is for base layer you already created
        dic_attribute_for_layer_vs_layer_id[list_attribute_for_layer[n]] = n + 1

    # make switched version
    dic_layer_id_vs_attribute_for_layer = {}

    for k, v in dic_attribute_for_layer_vs_layer_id.items():
        dic_layer_id_vs_attribute_for_layer[v] = k

    logger.warning(f'dic_layer_id_vs_attribute_for_layer: {dic_layer_id_vs_attribute_for_layer}')
    logger.debug(f"list_attribute_for_layer {list_attribute_for_layer}")
    logger.debug("finishing [define_layer]")

    return dic_layer_id_vs_attribute_for_layer
