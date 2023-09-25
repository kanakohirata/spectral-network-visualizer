import plotly
import sys


from logging import getLogger
logger = getLogger(__name__)


# stretching as function -----------------

def stretch_transform_list_float(l_o, new_min, new_max):

    # making as list of tuple
    l_o_tuple = []

    for i in range(0, len(l_o)):
        l_o_tuple.append((i, l_o[i]))

    l_o_tuple_sorted = sorted(l_o_tuple, key=lambda tup: tup[1], reverse=True)

    max_v_o = l_o_tuple_sorted[0][1]

    l_o_tuple_sorted = sorted(l_o_tuple, key=lambda tup: tup[1])
    min_v_o = l_o_tuple_sorted[0][1]

    l_o_tuple_sorted_flat = []

    for tup in l_o_tuple_sorted:
        l_o_tuple_sorted_flat.append((tup[0], tup[1] - min_v_o))

    l_tuple_idx_dif_r = []
    for i in range(0, len(l_o_tuple_sorted_flat)):
        if i == 0:
            # put 0 for origin, it wont do anything...
            dif_r = 0
        if i != 0:
            dif_r = float(l_o_tuple_sorted_flat[i][1]) / float(l_o_tuple_sorted_flat[-1][1])
        l_tuple_idx_dif_r.append((l_o_tuple_sorted_flat[i][0], dif_r))

    new_len_range = new_max - new_min

    l_transformed_tuple = []
    for tup_idx_dif_r in l_tuple_idx_dif_r:
        transformed = new_min + tup_idx_dif_r[1] * new_len_range
        l_transformed_tuple.append((tup_idx_dif_r[0], transformed))

    l_transformed_tuple_sorted = sorted(l_transformed_tuple, key=lambda tup: tup[0])

    l_stretch_transformed = []
    for tup in l_transformed_tuple_sorted:
        l_stretch_transformed.append(tup[1])

    return l_stretch_transformed


# now you have to transform

def rearrange_networkx_2d_layout_make_3d_zstat(layt_ori, x_min, x_max, y_min, y_max, z_stat):

    # networkx layout is dictionary
    l_id = []
    l_x = []
    l_y = []
    for key, l_x_y in layt_ori.items():
        l_id.append(key)
        l_x.append(l_x_y[0])
        l_y.append(l_x_y[1])

    l_x_new = stretch_transform_list_float(l_x, x_min, x_max)
    l_y_new = stretch_transform_list_float(l_y, y_min, y_max)
    l_z_new = [z_stat for l in layt_ori]

    l_layt_new = {}
    for n in range(0, len(l_x_new)):
        l_layt_new[l_id[n]] = ([l_x_new[n], l_y_new[n], l_z_new[n]])
    return l_layt_new




# This is for igraph layout !!!!!
# this function accept 2d network layout, then rearrange within given aerea, adding static z value
# this will end up with 2d network on mesh in 3d world.

def rearrange_2d_layout_make_3d_zstat(layt_ori,   x_min , x_max , y_min, y_max, z_stat ):
    l_x = [  l[0]      for l in layt_ori ]
    l_y = [  l[1]      for l in layt_ori ]

    l_x_new = stretch_transform_list_float(    l_x      , x_min,x_max)
    l_y_new = stretch_transform_list_float(  l_y       , y_min,y_max)
    l_z_new = [  z_stat     for l in layt_ori ]

    l_layt_new = []
    for n in range(0, len(l_x_new )):
        l_layt_new.append(   [ l_x_new[n] ,  l_y_new[n]  ,  l_z_new[n]]        )
    return l_layt_new


# this is for networkX

# making coordinates for nodes and edges BASED ON GIVEN LAYOUT. (this function wont change layout)
# making Xn,Yn,Zn for node coordinates,   Xe, Ye, Ze for edge coordinate (origin->target coordinate)



#def get_multi_traces_3d_network_x(layt, l_nodes, l_edges, l_node_color=[],  l_edges_color = [], dic_total_input_idx_mod_vs_node_obj = {}, l_edge_trace_idx = [] ):

def get_multi_traces_3d_network_x(list_dic_edges_nodes_graph_by_layer , list_of_edge_for_networkx_to_show_inter_sample_ref_layer, list_of_edge_for_networkx_to_show_inter_sample_layer):
    # l_node_color is simply a list of value of the node.  [0.1, 1.0, 1.2]
    #  l_edges_color  is simply a list of int that indicate color  [ 0, 1, 0, 1, 1]

    # note iterateing list of nodes, which consists of node id (could be name).
    # other function for igraph wont work same way...

    # !!!!
    # if there is input regarding changing node color based on quantification etc, takeover "node_color" arameter
    logger.debug("starting   multilayer_3d_rescale_functions_a2,get_multi_traces_3d_network_x")
    fo_log = open("log_multiple_traces.txt", "w")
    fo_log_layt = open("log_multiple_traces_layt.txt", "w")


    traces_all  = []
    dic_cluster_total_input_idx_MOD_vs_Xn_TOTAL = {}
    dic_cluster_total_input_idx_MOD_vs_Yn_TOTAL = {}
    dic_cluster_total_input_idx_MOD_vs_Zn_TOTAL = {}



    # this layout keeps all layout from all layers. This will be used for inter layer edges,
    layt_TOTAL = {}
    for dic in list_dic_edges_nodes_graph_by_layer:

        dic_cluster_total_input_idx_MOD_vs_node_info = dic["dic_cluster_total_input_idx_MOD_vs_node_info"]

        ##
        # !!!!!!!!!! checking if this is sample laers or not
        f_base_layer =0
        if dic["attribute_for_layer"] in ["sample", "sample1", "sample2", "sample3", "base", "base1", "base2", "base3",  "base4"]:
            f_base_layer = 1

        fo_log.write("\ndic[attribute_for_layer]:" + str(dic["attribute_for_layer"]) )
        fo_log.write("\ndic[attribute_for_layer]:" + str(dic ["layout_3d"]))

        l_nodes = list( dic["dic_cluster_total_input_idx_MOD_vs_node_info"].keys())
        layt = dic ["layout_3d"]


        fo_log_layt.write( "\n\n" +   str(layt))
        fo_log_layt.flush()


        # add layout of this layer to TOTAL layout
        for k,v in layt.items():
            layt_TOTAL[k]=v

        # get node color -------------

        l_node_color =  [ dic_cluster_total_input_idx_MOD_vs_node_info[ cluster_total_input_idx_MOD]   for  cluster_total_input_idx_MOD   in l_nodes         ]


        # get  edge color variation -------------
        l_edge_color = []
        l_node_color =  []
        for k, node_info in  dic["dic_cluster_total_input_idx_MOD_vs_node_info"].items():
            l_node_color.append( node_info.color )





        #MODIFIED20220714 ---------------------
        # change symbols if node match to supect mz list --------------------

        # first you have to figure out howmany suspect files were used----
        #  for example:     "my_suspect_file_A" : "x" , "my_suspect_file_B" : "triangle-up"
        dic_filename_suspect_vs_symbol = {}
        for k, node_info in  dic["dic_cluster_total_input_idx_MOD_vs_node_info"].items():
            if len(node_info.suspect) > 0 :
                dic_filename_suspect_vs_symbol[node_info.suspect] = ""

        l_plotly_symbols_for_suspect_file = ["x","diamond","cross","square"]

        # assign symbol to suspect filename ---------
        count = 0

        for k , v in dic_filename_suspect_vs_symbol.items():

            v = l_plotly_symbols_for_suspect_file[count]
            dic_filename_suspect_vs_symbol[k] = l_plotly_symbols_for_suspect_file[count]
            count = count + 1

        l_node_symbols = []

        for k, node_info in  dic["dic_cluster_total_input_idx_MOD_vs_node_info"].items():
            my_symbol = "circle"
            if len( node_info.suspect )> 0 :
                #my_symbol = "x"
                my_symbol = dic_filename_suspect_vs_symbol[node_info.suspect]
            l_node_symbols.append( my_symbol)

        # ------------------MODIFIED20220714



        list_edge_nx = dic["list_of_edge_for_networkx"]

        Xn = [layt[k][0] for k in l_nodes]
        Yn = [layt[k][1] for k in l_nodes]
        Zn = [layt[k][2] for k in l_nodes]


        # record laytout of all nodes. This will be used for INTER LAYER edges.
        for cluster_total_input_idx_MOD in l_nodes:
            dic_cluster_total_input_idx_MOD_vs_Xn_TOTAL[cluster_total_input_idx_MOD] = layt[cluster_total_input_idx_MOD][0]
            dic_cluster_total_input_idx_MOD_vs_Yn_TOTAL[cluster_total_input_idx_MOD] = layt[cluster_total_input_idx_MOD][1]
            dic_cluster_total_input_idx_MOD_vs_Zn_TOTAL[cluster_total_input_idx_MOD] = layt[cluster_total_input_idx_MOD][2]

        Xe = []
        Ye = []
        Ze = []

        dic_edge_tarce_idx_vs_l_Xe_coordinates = {}
        dic_edge_tarce_idx_vs_l_Ye_coordinates = {}
        dic_edge_tarce_idx_vs_l_Ze_coordinates = {}

        # make dictionary for each edge colors.   like dic of red edges, dic of black edge
        # if "l_edges_in_use_color"   is  [0,1,0,1]  ,  the color of edges with idx = 0,2 is  0 (black), and   edge with index = 1 and 3 is 1 (red)

        l_edge_Xe_coordinates = []
        l_edge_Ye_coordinates = []
        l_edge_Ze_coordinates = []
        l_edge_color = []
        for n in range(len(list_edge_nx)):
            # note   networkx "layout"  is a dictionary where key is node id and value is array of geometry of node.

            l_edge_Xe_coordinates += [   layt[list_edge_nx[n][0]][0], layt[list_edge_nx[n][1]][0], None   ]
            l_edge_Ye_coordinates += [   layt[list_edge_nx[n][0]][1], layt[list_edge_nx[n][1]][1], None   ]
            l_edge_Ze_coordinates += [   layt[list_edge_nx[n][0]][2], layt[list_edge_nx[n][1]][2], None   ]

            l_edge_color.append(list_edge_nx[n][2]["color"])
            l_edge_color.append(list_edge_nx[n][2]["color"])
            l_edge_color.append(list_edge_nx[n][2]["color"])
        l_nodes_str = [str(node) for node in l_nodes]


        l_nodes_str =  []

        for node in l_nodes :
            node_info =  dic_cluster_total_input_idx_MOD_vs_node_info[node]
            #  2.7   node_info_str = "node:" + node.decode('utf-8') + ";\n" + node_info.spec_cluster.global_accession.decode('utf-8') + ";\n" + node_info.spec_cluster.compound_name.decode('utf-8')
            node_info_str = "{" +   r"'NODE':" + node +  \
                            r",'PREC_MZ':" + str(node_info.spec_cluster.represen_spec_uni.precursor_mz) + \
                            r",'RT_SEC':" + str(node_info.spec_cluster.represen_spec_uni.retention_time_in_sec) + \
                            r",'TOTAL_INPUT_IDX_MOD':"  +  node_info.total_input_idx_mod + \
                            ";'SUSPECT_NAME':" + node_info.name + \
                            ";'COMPOUND_NAME':" + node_info.spec_cluster.compound_name + \
                            r",'GLOBAL_ACCESSION':" + node_info.spec_cluster.global_accession +\
                            "}"

            l_nodes_str.append(node_info_str)

        # making traces of edges..
        # making as trace (trace 1 is node, 2 is edges)

        l_edge_traces = []
        dic_id_vs_color = {  0 : "black" , 1: "red" , 2:"purple" , 2:"green"}


        #for n in range (len(list_edge_nx)):
        edge_trace = plotly.graph_objs.Scatter3d(x=l_edge_Xe_coordinates,
                                             y=l_edge_Ye_coordinates,
                                             z=l_edge_Ze_coordinates,
                                             mode='lines',
                                             line=dict(color=  l_edge_color, width=2),
                                             hoverinfo='none'
                                             )
        l_edge_traces.append( edge_trace)


        # https://stackoverflow.com/questions/58335021/plotly-how-to-set-width-to-specific-line

        # trace for nodes ==============================


        # setting node size
        if f_base_layer == 0 :
            node_size = 3
        if f_base_layer == 1 :
            node_size = 6

        l_node_trace = plotly.graph_objs.Scatter3d(x=Xn,
                                             y=Yn,
                                             z=Zn,
                                             mode='markers',
                                             name='COMPOUNDS',
                                             marker=dict(
                                                        # MODIFIED20220714-----
                                                        #symbol='circle',
                                                        symbol=l_node_symbols,
                                                        # ----- MODIFIED20220714
                                                         #size= node_size,
                                                         size=6,
                                                         # color=["red", "blue", "green", "gold", "purple"],

                                                        ## if color coded, node_color is actually a list.
                                                         color= l_node_color,
                                                         # see color map  https://awesomeopensource.com/project/bhaskarvk/colormap
                                                         #  https://community.plot.ly/t/what-colorscales-are-available-in-plotly-and-which-are-the-default/2079
                                                         # useful colorscale   'Viridis'  'YlGnBu', 'Portland'  'Rainbow','Hot'
                                                         # colorscale='Viridis',
                                                         colorscale='Hot',
                                                         line=dict(color='rgb(50,50,50)', width=0.5)
                                                         ),
                                             ## adding labels
                                             # text=labels,
                                             text=l_nodes_str,
                                             hoverinfo='text'
                                             )

        traces_all =  traces_all + l_edge_traces + [l_node_trace]


    fo_log.write("\n\nLAYT TOTAL"  + str(layt_TOTAL))
    fo_log.flush()

    # INTER LAYER EDGES ################

    fo_x = open("sssssssssssssssssss.txt", "w")
    fo_x.write(str(len(list_of_edge_for_networkx_to_show_inter_sample_layer)) + "\n")
    fo_x.write(str(str(list_of_edge_for_networkx_to_show_inter_sample_layer)) + "\n")

    #################################################
    ### INTER SAMPLE-REF EDGES --------------
    list_of_edge_for_networkx_to_show_inter_layer = list_of_edge_for_networkx_to_show_inter_sample_ref_layer + list_of_edge_for_networkx_to_show_inter_sample_layer

    #list_of_edge_for_networkx_to_show_inter_layer = list_of_edge_for_networkx_to_show_inter_sample_ref_layer

    #l_inter_e = list_of_edge_for_networkx_to_show_inter_sample_ref_layer
    l_inter_sample_ref_e = list_of_edge_for_networkx_to_show_inter_sample_ref_layer
    l_INTER_SAMPLE_REF_LAYER_edge_Xe_coordinates = []
    l_INTER_SAMPLE_REF_LAYER_edge_Ye_coordinates = []
    l_INTER_SAMPLE_REF_LAYER_edge_Ze_coordinates = []

    l_INTER_SAMPLE_REF_LAYER_edge_color = []

    for n in range(len(l_inter_sample_ref_e)):
    # note   networkx "layout"  is a dictionary where key is node id and value is array of geometry of node.
        l_INTER_SAMPLE_REF_LAYER_edge_Xe_coordinates += [layt_TOTAL[l_inter_sample_ref_e[n][0]][0], layt_TOTAL[l_inter_sample_ref_e[n][1]][0], None]
        l_INTER_SAMPLE_REF_LAYER_edge_Ye_coordinates += [layt_TOTAL[l_inter_sample_ref_e[n][0]][1], layt_TOTAL[l_inter_sample_ref_e[n][1]][1], None]
        l_INTER_SAMPLE_REF_LAYER_edge_Ze_coordinates += [layt_TOTAL[l_inter_sample_ref_e[n][0]][2], layt_TOTAL[l_inter_sample_ref_e[n][1]][2], None]
        l_INTER_SAMPLE_REF_LAYER_edge_color.append(l_inter_sample_ref_e[n][2]["color"])
        l_INTER_SAMPLE_REF_LAYER_edge_color.append(l_inter_sample_ref_e[n][2]["color"])
        l_INTER_SAMPLE_REF_LAYER_edge_color.append(l_inter_sample_ref_e[n][2]["color"])

    fo_x.write("l_INTER_LAYER_edge_color:" + str(l_INTER_SAMPLE_REF_LAYER_edge_color) + "\n")
    fo_x.write("l_INTER_LAYER_edge_Xe_coordinates:" + str(len(l_INTER_SAMPLE_REF_LAYER_edge_Xe_coordinates)) + "\n")
    fo_x.write("l_INTER_LAYER_edge_color:" + str(len(l_INTER_SAMPLE_REF_LAYER_edge_color)) + "\n")

    # for n in range (len(list_edge_nx)):
    INTER_SAMPLE_REF_LAYER_edge_trace = plotly.graph_objs.Scatter3d(x=l_INTER_SAMPLE_REF_LAYER_edge_Xe_coordinates,
                                             y=l_INTER_SAMPLE_REF_LAYER_edge_Ye_coordinates,
                                             z=l_INTER_SAMPLE_REF_LAYER_edge_Ze_coordinates,
                                             mode='lines',
                                             name='INTER_SAMPLE_REF_LAYER_EDGES',
                                             line=dict(color=l_INTER_SAMPLE_REF_LAYER_edge_color, width=1),
                                             hoverinfo='none'
                                             )

    ########################
    # inter sample layers ------------------
    #l_inter_e = list_of_edge_for_networkx_to_show_inter_sample_ref_layer
    l_inter_sample_e = list_of_edge_for_networkx_to_show_inter_sample_layer
    l_INTER_SAMPLE_LAYER_edge_Xe_coordinates = []
    l_INTER_SAMPLE_LAYER_edge_Ye_coordinates = []
    l_INTER_SAMPLE_LAYER_edge_Ze_coordinates = []

    l_INTER_SAMPLE_LAYER_edge_color = []

    for n in range(len(l_inter_sample_e)):
    # note   networkx "layout"  is a dictionary where key is node id and value is array of geometry of node.
        l_INTER_SAMPLE_LAYER_edge_Xe_coordinates += [layt_TOTAL[l_inter_sample_e[n][0]][0], layt_TOTAL[l_inter_sample_e[n][1]][0], None]
        l_INTER_SAMPLE_LAYER_edge_Ye_coordinates += [layt_TOTAL[l_inter_sample_e[n][0]][1], layt_TOTAL[l_inter_sample_e[n][1]][1], None]
        l_INTER_SAMPLE_LAYER_edge_Ze_coordinates += [layt_TOTAL[l_inter_sample_e[n][0]][2], layt_TOTAL[l_inter_sample_e[n][1]][2], None]
        l_INTER_SAMPLE_LAYER_edge_color.append(l_inter_sample_e[n][2]["color"])
        l_INTER_SAMPLE_LAYER_edge_color.append(l_inter_sample_e[n][2]["color"])
        l_INTER_SAMPLE_LAYER_edge_color.append(l_inter_sample_e[n][2]["color"])

    fo_x.write("l_INTER_LAYER_edge_color:" + str(l_INTER_SAMPLE_LAYER_edge_color) + "\n")
    fo_x.write("l_INTER_LAYER_edge_Xe_coordinates:" + str(len(l_INTER_SAMPLE_LAYER_edge_Xe_coordinates)) + "\n")
    fo_x.write("l_INTER_LAYER_edge_color:" + str(len(l_INTER_SAMPLE_LAYER_edge_color)) + "\n")

    # for n in range (len(list_edge_nx)):
    INTER_SAMPLE_LAYER_edge_trace = plotly.graph_objs.Scatter3d(x=l_INTER_SAMPLE_LAYER_edge_Xe_coordinates,
                                             y=l_INTER_SAMPLE_LAYER_edge_Ye_coordinates,
                                             z=l_INTER_SAMPLE_LAYER_edge_Ze_coordinates,
                                             mode='lines',
                                             name='INTER_SAMPLE_LAYER_EDGES',
                                             line=dict(color=l_INTER_SAMPLE_LAYER_edge_color, width=3),
                                             hoverinfo='none'
                                             )


    #traces_all = traces_all + [INTER_LAYER_edge_trace]
    traces_all = traces_all + [INTER_SAMPLE_REF_LAYER_edge_trace] + [INTER_SAMPLE_LAYER_edge_trace]

    logger.debug("finishing   multilayer_3d_rescale_functions_a2,get_multi_traces_3d_network_x")

    return traces_all



def get_rescaled_geometry_for_3d_network_x(list_dic_edges_nodes_graph_by_layer , list_of_edge_for_networkx_to_show_inter_sample_ref_layer, list_of_edge_for_networkx_to_show_inter_sample_layer):
    # l_node_color is simply a list of value of the node.  [0.1, 1.0, 1.2]
    #  l_edges_color  is simply a list of int that indicate color  [ 0, 1, 0, 1, 1]

    list_dic_edges_inter_sample_ref_layer =[]
    list_dic_edges_inter_sample_layer =[]




    # note iterateing list of nodes, which consists of node id (could be name).
    # other function for igraph wont work same way...

    # !!!!
    # if there is input regarding changing node color based on quantification etc, takeover "node_color" arameter

    fo_log = open("log_multiple_traces.txt", "w")
    fo_log_layt = open("log_multiple_traces_layt.txt", "w")


    traces_all  = []
    dic_cluster_total_input_idx_MOD_vs_Xn_TOTAL = {}
    dic_cluster_total_input_idx_MOD_vs_Yn_TOTAL = {}
    dic_cluster_total_input_idx_MOD_vs_Zn_TOTAL = {}



    # this layout keeps all layout from all layers. This will be used for inter layer edges,
    layt_TOTAL = {}
    for dic in list_dic_edges_nodes_graph_by_layer:

        dic_cluster_total_input_idx_MOD_vs_node_info = dic["dic_cluster_total_input_idx_MOD_vs_node_info"]

        ##
        # !!!!!!!!!! checking if this is sample laers or not
        f_base_layer =0
        if dic["attribute_for_layer"] in ["sample", "sample1", "sample2", "sample3", "base", "base1", "base2", "base3",  "base4"]:
            f_base_layer = 1

        fo_log.write("\ndic[attribute_for_layer]:" + str(dic["attribute_for_layer"]) )
        fo_log.write("\ndic[attribute_for_layer]:" + str(dic ["layout_3d"]))

        l_nodes = list( dic["dic_cluster_total_input_idx_MOD_vs_node_info"].keys())
        layt = dic ["layout_3d"]


        fo_log_layt.write( "\n\n" +   str(layt))
        fo_log_layt.flush()


        # add layout of this layer to TOTAL layout
        for k,v in layt.items():
            layt_TOTAL[k]=v

        # get node color -------------

        l_node_color =  [ dic_cluster_total_input_idx_MOD_vs_node_info[ cluster_total_input_idx_MOD]   for  cluster_total_input_idx_MOD   in l_nodes         ]


        # get  edge color variation -------------
        l_edge_color = []
        l_node_color =  []
        for k, node_info in  dic["dic_cluster_total_input_idx_MOD_vs_node_info"].items():
            l_node_color.append( node_info["color"] )





        #MODIFIED20220714 ---------------------
        # change symbols if node match to supect mz list --------------------

        # first you have to figure out how many suspect files were used----
        #  for example:     "my_suspect_file_A" : "x" , "my_suspect_file_B" : "triangle-up"
        dic_filename_suspect_vs_symbol = {}
        for k, node_info in  dic["dic_cluster_total_input_idx_MOD_vs_node_info"].items():
            if len(node_info["suspect"]) > 0 :
                dic_filename_suspect_vs_symbol[node_info["suspect"]] = ""

        l_plotly_symbols_for_suspect_file = ["diamond","diamond-open", "square","square-open","circle-open"]

        # assign symbol to suspect filename ---------
        count = 0
        for k , v in dic_filename_suspect_vs_symbol.items():
            v = l_plotly_symbols_for_suspect_file[count]
            dic_filename_suspect_vs_symbol[k] = l_plotly_symbols_for_suspect_file[count]
            count = count + 1


        l_node_symbols = []

        for k, node_info in  dic["dic_cluster_total_input_idx_MOD_vs_node_info"].items():
            my_symbol = "circle"

            # if the current node has supect match to ONE suspect file.
            if len( node_info["suspect"] )> 0  and  len(node_info["l_source_suspect"]) == 1 :
                #my_symbol = "x"
                my_symbol = dic_filename_suspect_vs_symbol[node_info["suspect"]]

            # if the current node has suspect matchesa to MULTIPLE suspect file
            if len( node_info["suspect"] )> 0  and  len(node_info["l_source_suspect"]) > 1 :
                #my_symbol = "x"
                my_symbol = "x"

            l_node_symbols.append( my_symbol)



        # ------------------MODIFIED20220714

        list_edge_nx = dic["list_of_edge_for_networkx"]

        Xn = [layt[k][0] for k in l_nodes]
        Yn = [layt[k][1] for k in l_nodes]
        Zn = [layt[k][2] for k in l_nodes]


        # record laytout of all nodes. This will be used for INTER LAYER edges.
        for cluster_total_input_idx_MOD in l_nodes:
            dic_cluster_total_input_idx_MOD_vs_Xn_TOTAL[cluster_total_input_idx_MOD] = layt[cluster_total_input_idx_MOD][0]
            dic_cluster_total_input_idx_MOD_vs_Yn_TOTAL[cluster_total_input_idx_MOD] = layt[cluster_total_input_idx_MOD][1]
            dic_cluster_total_input_idx_MOD_vs_Zn_TOTAL[cluster_total_input_idx_MOD] = layt[cluster_total_input_idx_MOD][2]

        Xe = []
        Ye = []
        Ze = []

        dic_edge_tarce_idx_vs_l_Xe_coordinates = {}
        dic_edge_tarce_idx_vs_l_Ye_coordinates = {}
        dic_edge_tarce_idx_vs_l_Ze_coordinates = {}

        # make dictionary for each edge colors.   like dic of red edges, dic of black edge
        # if "l_edges_in_use_color"   is  [0,1,0,1]  ,  the color of edges with idx = 0,2 is  0 (black), and   edge with index = 1 and 3 is 1 (red)

        l_edge_Xe_coordinates = []
        l_edge_Ye_coordinates = []
        l_edge_Ze_coordinates = []
        l_edge_color = []

        for n in range(len(list_edge_nx)):
            # note   networkx "layout"  is a dictionary where key is node id and value is array of geometry of node.
            #dic_edge_tarce_idx_vs_l_Xe_coordinates[n] += [   layt[list_edge_nx[n][0]][0], layt[list_edge_nx[n][1]][0], None   ]
            #dic_edge_tarce_idx_vs_l_Ye_coordinates[n] += [   layt[list_edge_nx[n][0]][1], layt[list_edge_nx[n][1]][1], None   ]
            #dic_edge_tarce_idx_vs_l_Ze_coordinates[n] += [   layt[list_edge_nx[n][0]][2], layt[list_edge_nx[n][1]][2], None   ]
            l_edge_Xe_coordinates += [   layt[list_edge_nx[n][0]][0], layt[list_edge_nx[n][1]][0], None   ]
            l_edge_Ye_coordinates += [   layt[list_edge_nx[n][0]][1], layt[list_edge_nx[n][1]][1], None   ]
            l_edge_Ze_coordinates += [   layt[list_edge_nx[n][0]][2], layt[list_edge_nx[n][1]][2], None   ]

            l_edge_color.append(list_edge_nx[n][2]["color"])
            l_edge_color.append(list_edge_nx[n][2]["color"])
            l_edge_color.append(list_edge_nx[n][2]["color"])
        l_nodes_str = [str(node) for node in l_nodes]



        l_nodes_str =  []

        for node in l_nodes :
            node_info =  dic_cluster_total_input_idx_MOD_vs_node_info[node]
            #  2.7   node_info_str = "node:" + node.decode('utf-8') + ";\n" + node_info.spec_cluster.global_accession.decode('utf-8') + ";\n" + node_info.spec_cluster.compound_name.decode('utf-8')
            node_info_str = "{" +   r"'NODE':" + node +  \
                            r",'PREC_MZ':" + str(node_info["spec_cluster"]["represen_spec_uni"]["precursor_mz"]) + \
                            r",'RT_SEC':" + str(node_info["spec_cluster"]["represen_spec_uni"]["retention_time_in_sec"]) + \
                            r",'TOTAL_INPUT_IDX_MOD':"  +  node_info["total_input_idx_mod"] + \
                            ";'SUSPECT_NAME':" + node_info["name"] + \
                            ";'COMPOUND_NAME':" + node_info["spec_cluster"]["compound_name"] + \
                            r",'GLOBAL_ACCESSION':" + node_info["spec_cluster"]["global_accession"] +\
                            "}"

            l_nodes_str.append(node_info_str)

        # making traces of edges..
        # making as trace (trace 1 is node, 2 is edges)

        l_edge_traces = []
        dic_id_vs_color = {  0 : "black" , 1: "red" , 2:"purple" , 2:"green"}

        dic["l_edge_Xe_coordinates"] =  l_edge_Xe_coordinates
        dic["l_edge_Ye_coordinates"] =  l_edge_Ye_coordinates
        dic["l_edge_Ze_coordinates"] =  l_edge_Ze_coordinates
        dic["l_edge_color"] =  l_edge_color
        #for n in range (len(list_edge_nx)):
        edge_trace = plotly.graph_objs.Scatter3d(x=l_edge_Xe_coordinates,
                                             y=l_edge_Ye_coordinates,
                                             z=l_edge_Ze_coordinates,
                                             mode='lines',
                                             line=dict(color=  l_edge_color, width=2),
                                             hoverinfo='none'
                                             )
        l_edge_traces.append( edge_trace)


        # https://stackoverflow.com/questions/58335021/plotly-how-to-set-width-to-specific-line

        # trace for nodes ==============================


        dic["l_node_Xn_coordinates"] = Xn
        dic["l_node_Yn_coordinates"] = Yn
        dic["l_node_Zn_coordinates"] = Zn
        dic["l_node_symbols"] = l_node_symbols
        dic["l_node_color"] = l_node_color
        dic["l_nodes_str"] = l_nodes_str

        # setting node size
        if f_base_layer == 0 :
            node_size = 3
        if f_base_layer == 1 :
            node_size = 6

        l_node_trace = plotly.graph_objs.Scatter3d(x=Xn,
                                             y=Yn,
                                             z=Zn,
                                             mode='markers',
                                             name='COMPOUNDS',
                                             marker=dict(
                                                        # MODIFIED20220714-----
                                                        #symbol='circle',
                                                        symbol=l_node_symbols,
                                                        # ----- MODIFIED20220714
                                                         #size= node_size,
                                                         size=6,
                                                         # color=["red", "blue", "green", "gold", "purple"],

                                                        ## if color coded, node_color is actually a list.
                                                         color= l_node_color,
                                                         # see color map  https://awesomeopensource.com/project/bhaskarvk/colormap
                                                         #  https://community.plot.ly/t/what-colorscales-are-available-in-plotly-and-which-are-the-default/2079
                                                         # useful colorscale   'Viridis'  'YlGnBu', 'Portland'  'Rainbow','Hot'
                                                         # colorscale='Viridis',
                                                         colorscale='Hot',
                                                         line=dict(color='rgb(50,50,50)', width=0.5)
                                                         ),
                                             ## adding labels
                                             # text=labels,
                                             text=l_nodes_str,
                                             hoverinfo='text'
                                             )

        traces_all =  traces_all + l_edge_traces + [l_node_trace]

    # end itearate graph by layer

    fo_log.write("\n\nLAYT TOTAL"  + str(layt_TOTAL))
    fo_log.flush()

    # INTER LAYER EDGES ################

    fo_x = open("sssssssssssssssssss.txt", "w")
    fo_x.write(str(len(list_of_edge_for_networkx_to_show_inter_sample_layer)) + "\n")
    fo_x.write(str(str(list_of_edge_for_networkx_to_show_inter_sample_layer)) + "\n")

    #################################################
    ### INTER SAMPLE-REF EDGES --------------
    list_of_edge_for_networkx_to_show_inter_layer = list_of_edge_for_networkx_to_show_inter_sample_ref_layer + list_of_edge_for_networkx_to_show_inter_sample_layer

    #list_of_edge_for_networkx_to_show_inter_layer = list_of_edge_for_networkx_to_show_inter_sample_ref_layer

    #l_inter_e = list_of_edge_for_networkx_to_show_inter_sample_ref_layer
    l_inter_sample_ref_e = list_of_edge_for_networkx_to_show_inter_sample_ref_layer
    l_INTER_SAMPLE_REF_LAYER_edge_Xe_coordinates = []
    l_INTER_SAMPLE_REF_LAYER_edge_Ye_coordinates = []
    l_INTER_SAMPLE_REF_LAYER_edge_Ze_coordinates = []

    l_INTER_SAMPLE_REF_LAYER_edge_color = []

    for n in range(len(l_inter_sample_ref_e)):
    # note   networkx "layout"  is a dictionary where key is node id and value is array of geometry of node.
        l_INTER_SAMPLE_REF_LAYER_edge_Xe_coordinates += [layt_TOTAL[l_inter_sample_ref_e[n][0]][0], layt_TOTAL[l_inter_sample_ref_e[n][1]][0], None]
        l_INTER_SAMPLE_REF_LAYER_edge_Ye_coordinates += [layt_TOTAL[l_inter_sample_ref_e[n][0]][1], layt_TOTAL[l_inter_sample_ref_e[n][1]][1], None]
        l_INTER_SAMPLE_REF_LAYER_edge_Ze_coordinates += [layt_TOTAL[l_inter_sample_ref_e[n][0]][2], layt_TOTAL[l_inter_sample_ref_e[n][1]][2], None]
        l_INTER_SAMPLE_REF_LAYER_edge_color.append(l_inter_sample_ref_e[n][2]["color"])
        l_INTER_SAMPLE_REF_LAYER_edge_color.append(l_inter_sample_ref_e[n][2]["color"])
        l_INTER_SAMPLE_REF_LAYER_edge_color.append(l_inter_sample_ref_e[n][2]["color"])

    fo_x.write("l_INTER_LAYER_edge_color:" + str(l_INTER_SAMPLE_REF_LAYER_edge_color) + "\n")
    fo_x.write("l_INTER_LAYER_edge_Xe_coordinates:" + str(len(l_INTER_SAMPLE_REF_LAYER_edge_Xe_coordinates)) + "\n")
    fo_x.write("l_INTER_LAYER_edge_color:" + str(len(l_INTER_SAMPLE_REF_LAYER_edge_color)) + "\n")

    dic_edges_inter_sample_ref_layer = {}
    dic_edges_inter_sample_ref_layer["l_INTER_SAMPLE_REF_LAYER_edge_Xe_coordinates"]= l_INTER_SAMPLE_REF_LAYER_edge_Xe_coordinates
    dic_edges_inter_sample_ref_layer["l_INTER_SAMPLE_REF_LAYER_edge_Ye_coordinates"] = l_INTER_SAMPLE_REF_LAYER_edge_Ye_coordinates
    dic_edges_inter_sample_ref_layer["l_INTER_SAMPLE_REF_LAYER_edge_Ze_coordinates"] =l_INTER_SAMPLE_REF_LAYER_edge_Ze_coordinates
    dic_edges_inter_sample_ref_layer["l_INTER_SAMPLE_REF_LAYER_edge_color"] = l_INTER_SAMPLE_REF_LAYER_edge_color

    # for n in range (len(list_edge_nx)):
    INTER_SAMPLE_REF_LAYER_edge_trace = plotly.graph_objs.Scatter3d(x=l_INTER_SAMPLE_REF_LAYER_edge_Xe_coordinates,
                                             y=l_INTER_SAMPLE_REF_LAYER_edge_Ye_coordinates,
                                             z=l_INTER_SAMPLE_REF_LAYER_edge_Ze_coordinates,
                                             mode='lines',
                                             name='INTER_SAMPLE_REF_LAYER_EDGES',
                                             line=dict(color=l_INTER_SAMPLE_REF_LAYER_edge_color, width=1),
                                             hoverinfo='none'
                                             )

    ########################
    # inter sample layers ------------------
    #l_inter_e = list_of_edge_for_networkx_to_show_inter_sample_ref_layer

    l_inter_sample_e = list_of_edge_for_networkx_to_show_inter_sample_layer
    l_INTER_SAMPLE_LAYER_edge_Xe_coordinates = []
    l_INTER_SAMPLE_LAYER_edge_Ye_coordinates = []
    l_INTER_SAMPLE_LAYER_edge_Ze_coordinates = []

    l_INTER_SAMPLE_LAYER_edge_color = []

    for n in range(len(l_inter_sample_e)):
    # note   networkx "layout"  is a dictionary where key is node id and value is array of geometry of node.
        l_INTER_SAMPLE_LAYER_edge_Xe_coordinates += [layt_TOTAL[l_inter_sample_e[n][0]][0], layt_TOTAL[l_inter_sample_e[n][1]][0], None]
        l_INTER_SAMPLE_LAYER_edge_Ye_coordinates += [layt_TOTAL[l_inter_sample_e[n][0]][1], layt_TOTAL[l_inter_sample_e[n][1]][1], None]
        l_INTER_SAMPLE_LAYER_edge_Ze_coordinates += [layt_TOTAL[l_inter_sample_e[n][0]][2], layt_TOTAL[l_inter_sample_e[n][1]][2], None]
        l_INTER_SAMPLE_LAYER_edge_color.append(l_inter_sample_e[n][2]["color"])
        l_INTER_SAMPLE_LAYER_edge_color.append(l_inter_sample_e[n][2]["color"])
        l_INTER_SAMPLE_LAYER_edge_color.append(l_inter_sample_e[n][2]["color"])

    fo_x.write("l_INTER_LAYER_edge_color:" + str(l_INTER_SAMPLE_LAYER_edge_color) + "\n")
    fo_x.write("l_INTER_LAYER_edge_Xe_coordinates:" + str(len(l_INTER_SAMPLE_LAYER_edge_Xe_coordinates)) + "\n")
    fo_x.write("l_INTER_LAYER_edge_color:" + str(len(l_INTER_SAMPLE_LAYER_edge_color)) + "\n")

    dic_edges_inter_sample_layer = {}
    dic_edges_inter_sample_layer["l_INTER_SAMPLE_LAYER_edge_Xe_coordinates"] = l_INTER_SAMPLE_LAYER_edge_Xe_coordinates
    dic_edges_inter_sample_layer["l_INTER_SAMPLE_LAYER_edge_Ye_coordinates"] = l_INTER_SAMPLE_LAYER_edge_Ye_coordinates
    dic_edges_inter_sample_layer["l_INTER_SAMPLE_LAYER_edge_Ze_coordinates"] = l_INTER_SAMPLE_LAYER_edge_Ze_coordinates
    dic_edges_inter_sample_layer["l_INTER_SAMPLE_LAYER_edge_color"] = l_INTER_SAMPLE_LAYER_edge_color

    # for n in range (len(list_edge_nx)):
    INTER_SAMPLE_LAYER_edge_trace = plotly.graph_objs.Scatter3d(x=l_INTER_SAMPLE_LAYER_edge_Xe_coordinates,
                                             y=l_INTER_SAMPLE_LAYER_edge_Ye_coordinates,
                                             z=l_INTER_SAMPLE_LAYER_edge_Ze_coordinates,
                                             mode='lines',
                                             name='INTER_SAMPLE_LAYER_EDGES',
                                             line=dict(color=l_INTER_SAMPLE_LAYER_edge_color, width=3),
                                             hoverinfo='none'
                                             )


    dic_rescaled_3d_network_data = {}

    dic_rescaled_3d_network_data["list_dic_edges_nodes_graph_by_layer"] = list_dic_edges_nodes_graph_by_layer
    dic_rescaled_3d_network_data["dic_edges_inter_sample_ref_layer"] = dic_edges_inter_sample_ref_layer
    dic_rescaled_3d_network_data["dic_edges_inter_sample_layer"] = dic_edges_inter_sample_layer

    return dic_rescaled_3d_network_data








def get_traces_3d_network_x(layt, l_nodes, l_edges, l_node_color =[] , l_edge_color =[]):
    # note iterateing list of nodes, which consists of node id (could be name).
    # other function for igraph wont work same way...


    # !!!!
    # if there is input regarding changing node color based on quantification etc, takeover "node_color" arameter
    node_color = "blue"
    if len(l_node_color) != 0 :
        node_color = l_node_color

    """
    for n in l_node_color:
        print type(n), n
    """

    Xn = [layt[k][0] for k in l_nodes]
    Yn = [layt[k][1] for k in l_nodes]
    Zn = [layt[k][2] for k in l_nodes]
    Xe = []
    Ye = []
    Ze = []

    for edge in l_edges:
        Xe += [layt[edge[0]][0], layt[edge[1]][0], None]
        Ye += [layt[edge[0]][1], layt[edge[1]][1], None]
        Ze += [layt[edge[0]][2], layt[edge[1]][2], None]

    l_nodes_str = [str(node) for node in l_nodes]
    # making as trace (trace 1 is node, 2 is edges)



    trace1 = plotly.graph_objs.Scatter3d(x=Xe,
                          y=Ye,
                          z=Ze,
                          mode='lines',

                          line=dict(color = 'rgb(125,125,125)' ,  width=8),
                          hoverinfo='none'
                          )


    #https://stackoverflow.com/questions/58335021/plotly-how-to-set-width-to-specific-line



    #color='rgb(125,125,125)',
    # trace for nodes
    trace2 = plotly.graph_objs.Scatter3d(x=Xn,
                          y=Yn,
                          z=Zn,
                          mode='markers',
                          name='COMPOUNDS',
                          marker=dict(symbol='circle',
                                      size=6,
                                      #color=["red", "blue", "green", "gold", "purple"],


                                      color = node_color,
                                      # see color map  https://awesomeopensource.com/project/bhaskarvk/colormap
                                      #  https://community.plot.ly/t/what-colorscales-are-available-in-plotly-and-which-are-the-default/2079
                                      # useful colorscale   'Viridis'  'YlGnBu', 'Portland'  'Rainbow','Hot'
                                      #colorscale='Viridis',
                                      colorscale='Hot',
                                      line=dict(color='rgb(50,50,50)', width=0.5)
                                      ),
                          ## adding labels
                          # text=labels,
                          text=l_nodes_str,
                          hoverinfo='text'
                          )
    return [trace1, trace2]


########################################
#  !!!!!!!!!!   This is for igpraph
#   making set data for network visualization
############################################

# Set data for the Plotly plot of the graph:
# using the coordinates made so far.

def get_traces_3d_network_igraph(layt):
    N = len(layt)
    # coordinates for nodes
    #print layt
    Xn = [layt[k][0] for k in range(N)]  # x-coordinates of nodes
    Yn = [layt[k][1] for k in range(N)]  # y-coordinates
    Zn = [layt[k][2] for k in range(N)]  # z-coordinates

    # coordinates for edges
    Xe = []
    Ye = []
    Ze = []
    for e in Edges:
        Xe += [layt[e[0]][0], layt[e[1]][0], None]  # x-coordinates of edge ends
        Ye += [layt[e[0]][1], layt[e[1]][1], None]
        Ze += [layt[e[0]][2], layt[e[1]][2], None]

    import plotly.plotly as py
    import plotly.graph_objs as go

    # trace for edges
    trace1 = go.Scatter3d(x=Xe,
                          y=Ye,
                          z=Ze,
                          mode='lines',
                          line=dict(color='rgb(125,125,125)', width=2),
                          hoverinfo='none'
                          )

    # trace for nodes
    trace2 = go.Scatter3d(x=Xn,
                          y=Yn,
                          z=Zn,
                          mode='markers',
                          name='COMPOUNDS',
                          marker=dict(symbol='circle',
                                      size=6,
                                      color=group,
                                      colorscale='Viridis',
                                      line=dict(color='rgb(50,50,50)', width=0.5)
                                      ),
                          ## adding labels
                          text=labels,
                          hoverinfo='text'
                          )


    return [trace1, trace2]