#  E. Hayakawa 2019 Feb
#  These functions are for 3d multilayer network visualization ( mostly for using plotly )

# this functiona is to fractionate range
# practical use:
#
# l_my_ranges_x = fractionate_range (  [ -10, 10] , 3 , ratio_gap = 0.02  )
# l_my_ranges_y = fractionate_range (  [ -10, 10] , 2  , ratio_gap = 0.02  )
# then you will have list of ranges that covers   3x2 layers within x: -10 to 10, y: -10 to 10
def fractionate_range(list_range, num_frac, ratio_gap=0):
    length_range = list_range[1] - list_range[0]
    print(length_range)
    length_frac = float(length_range) / float(num_frac)
    list_list_range_fractionated = []

    for i in range(0, num_frac):
        val_start = list_range[0] + length_frac * i
        val_end = list_range[0] + length_frac * (i + 1)

        list_list_range_fractionated.append(
            [val_start + (length_range / num_frac * ratio_gap), val_end - (length_range / num_frac * ratio_gap)])

    list_list_range_fractionated[len(list_list_range_fractionated) - 1][1] = list_range[1] - (
    length_range / num_frac * ratio_gap)

    return list_list_range_fractionated


# this function is make list of depth in z axis.
# say you want to 3 levels in Z axis within z : -10 to 10, you can use like
#   list_z = fractionate_z_flat( [  -10 , 10  ], 3)
# list_range is the range in which you want to split z .  for instance if you want to split within -1 to 1, list_range = [-1,1]
# num_frac is number of fractions.  if you want to split [-1,1] into 3,(you get 3 meth layers between -1, 1), num_frac = 3

def fractionate_z_flat(list_range, num_frac):
    print ("starting  multilayer_3d_mesh_functions.fractionate_z_flat")
    length_range = list_range[1] - list_range[0]
    print("length_range", length_range)
    print("num_frac" , num_frac)

    #MODIFIED20220714-------------
    if num_frac > 1 :
        length_frac = float(length_range) / float(num_frac - 1)
        #length_frac = float(length_range) / float(num_frac)

        list_z = []

        for i in range(0, num_frac):
            list_z.append(list_range[0] + length_frac * i)
    if num_frac == 1 :
        list_z = [ list_range[0]]
    # ----------------------MODIFIED20220714

    print("returning  list_Z at the end of multilayer_3d_mesh_functions.fractionate_z_flat")
    return list_z



# this function is to create list of ranges for x, y and depth z
# for instance



#
# receive layer id and assign 3d mesh.
def make_dic_of_zflat_3d_mesh_coordinates(l_layer_id, l_ranges_x, l_ranges_y, l_z):
    dic_layer_id_vs_zflat_3d_mesh_coordinates = {}
    l_3d_coordinates = []
    for zi in range(0, len(l_z)):
        for yi in range(0, len(l_ranges_y)):
            for xi in range(0, len(l_ranges_x)):

                mesh_coord_x = [l_ranges_x[xi][0], l_ranges_x[xi][1], l_ranges_x[xi][0], l_ranges_x[xi][1]]
                mesh_coord_y = [l_ranges_y[yi][0], l_ranges_y[yi][0], l_ranges_y[yi][1], l_ranges_y[yi][1]]
                mesh_coord_z = [l_z[zi], l_z[zi], l_z[zi], l_z[zi]]

                # make coordiates for 3d (fub z-flat) mesh

                l_3d_coordinates.append([mesh_coord_x, mesh_coord_y, mesh_coord_z])

    count = 0

    for layer_id in l_layer_id:
        dic_layer_id_vs_zflat_3d_mesh_coordinates[layer_id] = l_3d_coordinates[count]
        count = count + 1

    return dic_layer_id_vs_zflat_3d_mesh_coordinates