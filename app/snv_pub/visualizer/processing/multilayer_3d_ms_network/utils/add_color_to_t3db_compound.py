from logging import getLogger


logger = getLogger(__name__)


def add_color_to_t3db_compound(is_colored, dic_cluster_total_input_idx_vs_cluster_info):
    logger.info('Start add_color_to_t3db_compound()')
    if not is_colored:
        for input_idx, cluster_info in dic_cluster_total_input_idx_vs_cluster_info.items():
            cluster_info['extra_node_color'] = ''
    else:
        for input_idx, cluster_info in dic_cluster_total_input_idx_vs_cluster_info.items():
            if cluster_info['is_toxic']:
                cluster_info['extra_node_color'] = '#00e08e'
