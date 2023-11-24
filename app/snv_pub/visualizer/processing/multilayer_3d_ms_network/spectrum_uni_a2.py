# class spectrun_uni_a1
# this class is intended for the universal use of spectrum information in my projects.
# This class may be inherited to read/handle data from database (massbank) or experimental data

# last modified 201801219
#
#  et fragment structure annotation implemented
#
#



# compound_uni now has some features from metacycc compound


# How to use ########################################################################

# [Create spectrum uni object from single file]
# you want to probably use some external module to read spectrum file
# this case "read_massbank_3_a1" is reading "UT002689.txt" spectrum file and creating spectrum uni object
#  spec = read_massbank_3_a1.make_spec_uni_obj_from_massbank_file("UT002689.txt",1)

# # [Create list of spectrum uni object from folder contatining spectra data set ]
#  list_spec = read_massbank_3_a1.generate_list_of_spectum_uni_object_from_folder('E:\\Dropbox\\Dropbox\\Programming_drpbx\\python\\mass_spec_related\\read_massbank_3\\test_files\\*.txt')


# make relative intensity peak list
#
# for most spec data in massbank has relative intensity data, but some does not.
# such case you have to make peak list in relative intensity
#
# This example is making relative intensity spectrum with max value of 999
# spectrum_uni_a1.set_peak_list_rel_int( my_spec_uni_obj , 999)



# update/complement precursor mz
# often massbank data does not contain precursor m/z, although there is data of exact mass which is not ionized form.
# this function add calculated precursor mz to the spec uni object, based on ion types
#   update_precursor_mz_based_on_exact_mass(spec_uni)



# [clean potential contamination peaks]
# its often the case that spectra from a machine/institute contain common contaminant peaks.
# let get rid of those
# list_cleaned = spectrum_uni_a1.clean_common_contami_peaks_author_inst(list_spec  , 0.2,0.3, 2)
#
#  1st arg : list_spec is the list of spec uni object
# 2nd arg:




# [deisotope]
#  remove obvious isotopes (2nd or 3rd)
#  use "deisotope_simple_for_peaklist_a2_for_uni" function by providing with peak list [[mz,int],[]] and register as new peak list
#  Note you have to do it for peak list in absolute intensity and relative intensity modes.
#  You should update the number of peaks info in the object as well

# for spec_uni in list_spec :
#    spec_uni.peak_list_mz_int_abs =      spectrum_uni_a1.deisotope_simple_for_peaklist_a2_for_uni(  spec_uni.peak_list_mz_int_abs  , 0.1 , 3)
#    spec_uni.peak_list_mz_int_rel =      spectrum_uni_a1.deisotope_simple_for_peaklist_a2_for_uni(  spec_uni.peak_list_mz_int_rel ,0.1 , 3 )
#    sp.num_peaks = len(sp.peak_list_mz_int_rel)

# 1st arg: peaklist   2nd arg: mz tolerance to find isotope peak    3rd arg: relative intensity to judge potential isotope. If you set 3, 2nd isotope-like peak whose intensity is less than 1/3 of monoiso peak will be removed



# remove less intense peak within ranges
#   say you have intense peak 120 mz and weak peaks around it, you want to get rid of these within a range....
#
#   If you want to keep 2 most intense peaks with bins of 100 mz,
#   spectrum_uni_a1.set_topN_most_intense_peak_list_in_binned_rangesW( my_spectrum_uni_obj , 2 ,100)



# convert intensity to log(1+x).  Not like simple log, this prevent having negative valye
# USAGE
#spectrum_uni_a1.set_peaklist_intensity_in_log1p(spec,999)
# it does not return anything. change the content of input spectrum_uni object
# 1st arg : spectrum uni object   2nd arg: maximum value of relative intensity. Normally 999



# [convert intensity to log10 form]   "set_peaklist_intensity_in_log10"
# This function convert intensity of peak list in aboslute to log10 format.
# This consequently change peak list in relative intensity format, so you need tell what will be the maximum value for ralative intensity
# USAGE
#spectrum_uni_a1.set_peaklist_intensity_in_log10(spec,999)
# it does not return anything. change the content of input spectrum_uni object
# 1st arg : spectrum uni object   2nd arg: maximum value of relative intensity


###   introduce mass shift (delta)
#       adding varying delta mz value (range of 0 - 50.) in int to each peaks.
#       Both absolute intensity list and relative intensity list will affected in orchestrated manner.
#    spectrum_uni_a1.introduce_random_delta_to_spec_uni(o1,0,50)


import math
import sys


class spectrum_uni_class:
    def __init__(self):
            # source is where the data comes from
            self.source = ""

            self.accession_number = ""
            self.unique_number = 0
            #scan number can be string such as 26-611
            self.scan_number = ""

            self.authors =""

            self.name = ""
            self.cas_no = ""
            self.inchi_key =""
            self.smiles =""
            self.inchi =""

            self.pubchem_cid    = ""
            self.KEGG_id =""
            self.id    =  ""

            self.author = ""
            self.institution = ""
            self.ms_level = ""
            self.spectrum_type = ""
            self.ionization_mode =""

            self.instrument = ""
            self.instrument_type = ""
            self.fragmentation_type = ""
            self.fragmentation_energy = 0

            self.mass_resolution_class = 0

            self.exact_mass    =     0
            self.precursor_mz    =    0

            self.retention_time_in_sec = None
            self.retention_time_in_min = None
            #self.peak_list_mz =[]
            #self.peak_list_int =[]
            self.peak_list_mz_int_abs =[]
            self.peak_list_mz_int_rel = []


            # for fragment annotation
            self.peak_list_mz_annoid= []
            #  [    [100,[1,2]] , [200, [3,4]]      ]   # means frag 100mz , potential structure id is 1 and 2.
            self.dic_annoid_struct = {}
            #   { 1 : C=CCC#[CH2+] }   { 2:C=C(C)CC[OH2+]     }   # dictionary  key is annoation id, value is smiles, inchi or inchikey
            self.type_anno_struct = ""

            self.num_peaks = []


            self.charge_type    =    ""
            self.precursor_type = ""

            self.flag_rdkit_compatible = 0
            self.source_txt_data = ""
            self.source_txt_data_till_peak_info = ""



            self.list_cmpd_classification_kingdom = []
            self.list_cmpd_classification_class = []
            self.list_cmpd_classification_superclass = []
            self.list_cmpd_classification_alternative_parent = []

            self.list_cmpd_pathway = []

            self.source_filename =""
            self.source_filepath =""


# this class express compound in MS friendly manner.
# for instance, this class holds multiple spectra for single chemical compound entity.
# Say you have 3 MS2 spec for serotonin, then object for serotonin holds 3 MS2 spectra


class compound_uni_class:
    def __init__(self):

            self.name_represen = ""
            self.cas_no = ""
            self.inchi_key =""
            self.smiles =""
            self.inchi =""


            self.cas_no = ""
            self.pubchem_cid    = 0
            self.KEGG_id =""

            self.id    =  ""
            self.exact_mass    =     ""

            self.list_spectrum_uni  = []

            self.Pathway_metacyc =""

            self.Reaction_equation_metacyc = ""
            self.EC_metacyc =""

            self.list_classification_kingdom = []
            self.list_classification_class = []
            self.list_classification_superclass = []
            self.list_classification_alternative_parent = []



import random


############
#  Edit mass spectrum
############


# this function intentionally make peak list whose intensity is relative value
# You may need this function when reading ms2 spec file that only has absolute intensity.

def set_peak_list_rel_int(spec_uni,max_value=100):

    # empty to make sure
    spec_uni.peak_list_mz_int_rel =[]

    spec_uni.peak_list_mz_int_abs.sort(key=lambda e:e[1],reverse = True)

    if ( len(spec_uni.peak_list_mz_int_abs) >0 ):

        highest_int = spec_uni.peak_list_mz_int_abs[0][1]

        spec_uni.peak_list_mz_int_abs.sort(key=lambda e:e[0])

        for pk in spec_uni.peak_list_mz_int_abs:
            normalized_int =   (    float(pk[1]) / float (highest_int)  )   *  max_value
            spec_uni.peak_list_mz_int_rel.append(   [  pk[0]  ,  normalized_int  ]  )





def set_peaklist_intensity_in_log10(spec_uni,max_value_of_rel):


    list_peaklist_abs_log10 =[]
    for pk in spec_uni.peak_list_mz_int_abs:
        list_peaklist_abs_log10.append(  [ pk[0] ,   math.log10(pk[1]) ]   )

    spec_uni.peak_list_mz_int_abs = list_peaklist_abs_log10

    # make peak list in relative instensity++++++++++++++
    new_peaklist_rel =[]

    # sort by intensity
    spec_uni.peak_list_mz_int_abs.sort(key=lambda e:e[1],reverse = True)

    if ( len(spec_uni.peak_list_mz_int_abs) >0 ):
        highest_int = spec_uni.peak_list_mz_int_abs[0][1]

        spec_uni.peak_list_mz_int_abs.sort(key=lambda e:e[0])

        for pk in spec_uni.peak_list_mz_int_abs:
            normalized_int =   (    float(pk[1]) / float (highest_int)  )   *  max_value_of_rel
            new_peaklist_rel.append(   [  pk[0]  ,  normalized_int  ]  )

        spec_uni.peak_list_mz_int_rel = new_peaklist_rel


def set_peaklist_intensity_in_sqrt(spec_uni,max_value_of_rel):


    list_peaklist_abs_sqrt =[]
    for pk in spec_uni.peak_list_mz_int_abs:
        list_peaklist_abs_sqrt.append([pk[0],   math.sqrt(pk[1])])

    spec_uni.peak_list_mz_int_abs = list_peaklist_abs_sqrt

    # make peak list in relative instensity++++++++++++++
    new_peaklist_rel =[]


    # now you make relative intensity peak list
    # sort by intensity to get the highest intensity value
    spec_uni.peak_list_mz_int_abs.sort(key=lambda e:e[1],reverse = True)

    if ( len(spec_uni.peak_list_mz_int_abs) >0 ):
        highest_int = spec_uni.peak_list_mz_int_abs[0][1]

        # sort again.
        spec_uni.peak_list_mz_int_abs.sort(key=lambda e:e[0])

        for pk in spec_uni.peak_list_mz_int_abs:
            normalized_int =   (    float(pk[1]) / float (highest_int)  )   *  max_value_of_rel
            new_peaklist_rel.append(   [  pk[0]  ,  normalized_int  ]  )

        spec_uni.peak_list_mz_int_rel = new_peaklist_rel



def set_peaklist_intensity_in_log1p(spec_uni,max_value_of_rel):


    list_peaklist_abs_log1p =[]
    for pk in spec_uni.peak_list_mz_int_abs:
        list_peaklist_abs_log1p.append(  [ pk[0] ,   math.log1p(pk[1]) ]   )

    spec_uni.peak_list_mz_int_abs = list_peaklist_abs_log1p

    # make peak list in relative instensity++++++++++++++
    new_peaklist_rel =[]

    # sort by intensity
    spec_uni.peak_list_mz_int_abs.sort(key=lambda e:e[1],reverse = True)

    if ( len(spec_uni.peak_list_mz_int_abs) >0 ):
        highest_int = spec_uni.peak_list_mz_int_abs[0][1]

        spec_uni.peak_list_mz_int_abs.sort(key=lambda e:e[0])

        for pk in spec_uni.peak_list_mz_int_abs:
            normalized_int =   (    float(pk[1]) / float (highest_int)  )   *  max_value_of_rel
            new_peaklist_rel.append(   [  pk[0]  ,  normalized_int  ]  )

        spec_uni.peak_list_mz_int_rel = new_peaklist_rel








def get_peak_list_rel_int(peak_list,max_value):

    peak_list.sort(key=lambda e:e[1],reverse = True)
    new_peak_list=[]
    if ( len(peak_list) >0 ):

        highest_int = peak_list[0][1]

        peak_list.sort(key=lambda e:e[0])

        for pk in peak_list:
            #print pk
            normalized_int =   (    float(pk[1]) / float (highest_int)  )   *  max_value
            #print pk[0]  ,  normalized_int

            new_peak_list.append(   [  pk[0]  ,  normalized_int  ]  )

    return new_peak_list




#############
# remove low intense peaks
################
"""
if you want to remove peaks with lower than 5% of the most intense peak in the spectrum,
remove_low_int_peaks( my_spec_uni_obj, cutoff = 0.05 )

both peaklist with absolute intensity and relateive intensity will be processed

!!ATTENTION!!!
You need to have peaklist BOTH in abolute and  relative intensity ("peak_list_mz_int_rel" and "peak_list_mz_int_abs")


"""
def remove_low_int_peaks(  spec_uni , cutoff = 0.05        ):

    new_peak_list_mz_int_abs = []
    new_peak_list_mz_int_rel = []

    if len(spec_uni.peak_list_mz_int_abs ) == 0 :
        print("ERROR!!!   at spectrum)uni_a2, remove_low_int_peaks function. it seems your peak_list_mz_int_abs has 0 length. This function uses peak_list_mz_int_abs to edit peak list")
        sys.exit()

    spec_uni.peak_list_mz_int_abs.sort(key=lambda e:e[1], reverse = True)

    highest_int = spec_uni.peak_list_mz_int_abs[0][1]

    # sort again
    spec_uni.peak_list_mz_int_abs.sort(key=lambda e: e[0])

    for pk in spec_uni.peak_list_mz_int_abs:

        if float(pk[1]   >=  ( float(highest_int)) *  cutoff ) :
            # print pk[0] , pk[1] , "is retained"  ,  "highest is " , highest_int, "cutoff is" ,  float(( float(highest_int)) *  cutoff )
            new_peak_list_mz_int_abs.append(pk)


    if len(spec_uni.peak_list_mz_int_rel ) == 0 :
        print("ERROR!!!   at spectrum)uni_a2, remove_low_int_peaks function. it seems your peak_list_mz_int_abs has 0 length. This function uses peak_list_mz_int_abs to edit peak list")
        sys.exit()

    spec_uni.peak_list_mz_int_rel.sort(key=lambda e:e[1], reverse = True)

    highest_int = spec_uni.peak_list_mz_int_rel[0][1]

    # sort again
    spec_uni.peak_list_mz_int_rel.sort(key=lambda e: e[0])

    for pk in spec_uni.peak_list_mz_int_rel:

        if float(pk[1]   >=  ( float(highest_int)) *  cutoff ) :
            # print pk[0] , pk[1] , "is retained"  ,  "highest is " , highest_int , "cutoff is" ,  float(( float(highest_int)) *  cutoff )
            new_peak_list_mz_int_rel.append(pk)

    # update spectrum data
    spec_uni.peak_list_mz_int_abs = new_peak_list_mz_int_abs
    spec_uni.peak_list_mz_int_rel = new_peak_list_mz_int_rel




# usage   (it does not return object)
#    remove_precursor_peak_from_spectrum_uni(spec_uni , 2)

def remove_precursor_peak_from_spectrum_uni(spec_uni , mz_range , max_value_rel_int = 100):
    #

    prec_mz = spec_uni.precursor_mz

    new_peaklist_mz_int_abs =[]
    new_peaklist_mz_int_rel =[]

    # if this object has peaklist of absolute intensity
    if(  len (spec_uni.peak_list_mz_int_abs ) > 0):


        for pk in spec_uni.peak_list_mz_int_abs:
            # if  the current mz is far enough from  precursor mz     or  precursor mz +1  or precursor mz +2
            if (   abs(pk[0] - prec_mz) > mz_range     ) :
                new_peaklist_mz_int_abs.append(pk)

        spec_uni.peak_list_mz_int_abs = new_peaklist_mz_int_abs

        # then create relative intensity peaklist accordingly
        spec_uni.peak_list_mz_int_rel = get_peak_list_rel_int (new_peaklist_mz_int_abs, max_value_rel_int)

    '''
    for pk in spec_uni.peak_list_mz_int_rel:
        # if  the current mz is far enough from  precursor mz     or  precursor mz +1  or precursor mz +2
        if (   abs(pk[0] - prec_mz) > mz_range    ) :
            new_peaklist_mz_int_rel.append(pk)
    '''

    # if this object DOES NOT have peaklist in abosolute, but have peaklist in relative intensity
    if(  len (spec_uni.peak_list_mz_int_abs ) ==  0   and len( spec_uni.peak_list_mz_int_rel  ) > 0):


        for pk in spec_uni.peak_list_mz_int_rel:
            # if  the current mz is far enough from  precursor mz     or  precursor mz +1  or precursor mz +2
            if (   abs(pk[0] - prec_mz) > mz_range     ) :
                new_peaklist_mz_int_abs.append(pk)

        spec_uni.peak_list_mz_int_rel = get_peak_list_rel_int (new_peaklist_mz_int_abs, max_value_rel_int)


    #refresh number of peaks
    spec_uni.num_peaks = len(spec_uni.peak_list_mz_int_abs )





def update_precursor_mz_based_on_exact_mass(spec_uni):

    prec_mz = 0
    if(spec_uni.precursor_type == "[M+H]+"):
        prec_mz = spec_uni.exact_mass + 1.00728
    if(spec_uni.precursor_type == "[M+Na]+"):
        prec_mz = spec_uni.exact_mass + 22.989218
    if(spec_uni.precursor_type == "[M+NH4]+"):
        prec_mz = spec_uni.exact_mass + 18.033823

    if(spec_uni.precursor_type == "[M-H]-"):
        prec_mz = spec_uni.exact_mass - 1.00728
    if(spec_uni.precursor_type == "M+FA-H"):
        prec_mz = spec_uni.exact_mass - 44.998201
    if(spec_uni.precursor_type == "M+Hac-H"):
        prec_mz = spec_uni.exact_mass - 59.013851

    spec_uni.precursor_mz = prec_mz




def introduce_mass_error_to_peak_list_by_ppm(peak_list, error_ppm):



    peak_list_mod = []

    # for each peak
    for pk in peak_list:

        # error value to introduce
        mz_error_val =  (pk[0] * error_ppm ) / 1000000

        mz_error_val_to_introducde =  random.uniform( mz_error_val * -1.00   ,mz_error_val)


        # mz value modded
        mod_mz   = pk[0] + mz_error_val_to_introducde
        # modded peak
        peak_mod = [mod_mz, pk[1]]
        # add to new list
        peak_list_mod.append( peak_mod  )

    # return value
    return peak_list_mod




def introduce_mass_error_to_peak_list_by_mz(peak_list, error_mz):

    peak_list_mod = []

    # for each peak
    for pk in peak_list:

        # error value to introduce
        mz_error_val_to_introducde =  random.uniform( error_mz * -1.00   ,error_mz)
        #print mz_error_val_to_introducde
        # mz value modded
        mod_mz   = pk[0] + mz_error_val_to_introducde
        # modded peak
        peak_mod = [mod_mz, pk[1]]
        # add to new list
        peak_list_mod.append( peak_mod  )

    # return value
    return peak_list_mod



#########
# deisotope
#########

def deisotope_simple_for_peaklist_a2_for_uni(peak_list,mz_tol,int_ratio)    :


    list_mz_original = [ ]
    list_int_original = [ ]

    for pk in peak_list:
        list_mz_original.append(pk[0])
        list_int_original.append(pk[1])

    dic_mz_int = {}

    for n in range ( 0 , len(list_mz_original)):

        flag_to_be_removed = 0

        idx = len(list_mz_original) - n - 1

        curr_mz = list_mz_original[idx]
        curr_int = list_int_original[idx]

        for extend in range(0,10):

            if (  idx - extend < 0):
                break
            if (  list_mz_original[idx - extend]  <  curr_mz -2 ):
                break

            fabs_mz_curr_target =   math.fabs(   list_mz_original[idx - extend]    - curr_mz + 1 )


            # if the mz looks "isotope"ish enough, and the intensity of current peak is smaller enough...
            if(   fabs_mz_curr_target < mz_tol  and    list_int_original[idx - extend]  >  curr_int * int_ratio            ):

                flag_to_be_removed =  1
                #print "deisotoping",  curr_mz , "(", curr_int   , ")  since its likely tb the isotope of" , list_mz_original[idx - extend] , "(" , list_int_original[idx - extend] , ")"

        #if the current peak does not look like for isotope of any peaks in search...
        if (flag_to_be_removed ==  0 ):

            dic_mz_int[curr_mz] = curr_int


    list_mz_res =[]
    list_int_res =[]

    peak_list_deisotoped = []
    import collections
    od = collections.OrderedDict(sorted(dic_mz_int.items()))
    for key, value in (iter(od.items())):
        #print key, value
        list_mz_res.append(key)
        list_int_res.append(value)
        peak_list_deisotoped.append([key,value])

    '''
    print "res"
    for n in range(0,len(list_mz_res)):
        print list_mz_res[n] , list_int_res[n]
    '''

    return peak_list_deisotoped



import copy

# this func take peak list [ [100,10] , [200,5]] and return top N peak list
def make_topN_most_intense_peak_list(peak_list, topN):

    peak_list_mod = copy.deepcopy(peak_list)
    peak_list_mod.sort(key=lambda e:e[1], reverse=True)

    list_temp =[]
    count_n = 0

    #iterate list sorted with 3rd value
    # til reaches top_n
    for obj in peak_list_mod:
        list_temp.append(obj)
        count_n = count_n + 1
        if( count_n == topN):break

    list_temp.sort(key=lambda e:e[0])

    return list_temp


#--------------------------------------------------------
'''
spec_obj = read_mgf_python.ReadMGFfile("GGFMTSEKSQTP.mgf")[0]
peak_list_mod = introduce_mass_error_to_peak_list_by_ppm(spec_obj.peak_list , 500)
print spec_obj.peak_list

print "------------"
print spec_obj.peak_list
peak_list_mod_mz = introduce_mass_error_to_peak_list_by_mz(spec_obj.peak_list , 0.1)
print peak_list_mod_mz
spec_obj.peak_list = peak_list_mod_mz
print spec_obj.peak_list

print make_topN_most_intense_peak_list(spec_obj.peak_list , 5)
'''


import os


import glob



### Select/Extract spectrum ++++++++++++++++++++

# by spectrum_type
def select_list_of_spectrum_uni_object_by_sepctrum_type (list_spectrum_uni_object , key ):

    list_return = []
    for obj in list_spectrum_uni_object:
        if(obj.spectrum_type == key):
            list_return.append(obj)

    return list_return


# by ionization_mode
def select_list_of_spectrum_uni_object_by_ionization_mode (list_spectrum_uni_object , key ):

    list_return = []
    for obj in list_spectrum_uni_object:
        if(obj.ionization_mode == key):
            list_return.append(obj)
    return list_return

# by instrument_type
def select_list_of_spectrum_uni_object_by_instrument_type (list_spectrum_uni_object , key ):

    list_return = []
    for obj in list_spectrum_uni_object:
        if(obj.instrument_type == key):
            list_return.append(obj)

    return list_return

# by fragmentation_type
def select_list_of_spectrum_uni_object_by_fragmentation_type (list_spectrum_uni_object , key ):

    list_return = []
    for obj in list_spectrum_uni_object:
        if(obj.fragmentation_type == key):
            list_return.append(obj)

    return list_return


def select_list_of_spectrum_uni_object_by_precursor_type (list_spectrum_uni_object , key ):

    list_return = []
    for obj in list_spectrum_uni_object:
        if(obj.precursor_type == key):
            list_return.append(obj)

    return list_return

###########################
#  make non-redundant spectral dataset.
############################


# usage
# if you want to remove redandant spectra that has 50 % peak overe lap with 100 ppm
# list_sp_nr = spectrum_uni_a1.make_list_of_spec_uni_non_redundant(list_sp, 100, 0.5)

def make_list_of_spec_uni_non_redundant(     list_spectrum_uni_object    , mass_tol_ppm ,ratio_pk_b_explained  ) :


    # iterate list of spec uni
    list_idx_tb_deleted = []

    for idx_a in range( 0,  len(list_spectrum_uni_object)) :
        flag_idx_a_is_redundant = 0
        #print "\n\nidx_a", idx_a
        #for spec obj A
        sp_obj_a = list_spectrum_uni_object[idx_a]
        #print "prec mz", sp_obj_a.precursor_mz

        # iterate same lis of spec uni
        for idx_b in range( 0,  len(list_spectrum_uni_object)) :
            #print "---idx_b", idx_b

            if(idx_a != idx_b):
                # for spec obj B
                sp_obj_b = list_spectrum_uni_object[idx_b]

                # match prec m/z
                # if prec m/z matched,
                if(      math.fabs( (  sp_obj_a.precursor_mz - sp_obj_b.precursor_mz  ) /  sp_obj_a.precursor_mz )  * 1000000  <  mass_tol_ppm ) :
                        #print "------------prec mz matched" ,  sp_obj_a.precursor_mz , sp_obj_b.precursor_mz
                        count_peak_match = 0
                        # iterate exp frag peaks , but use n
                        for n_a in range( 0, len( sp_obj_a.peak_list_mz_int_abs )):

                            #print "pk_a" ,  sp_obj_a.peak_list_mz_int_abs[n_a][0]
                            # iterate theo frag peaks
                            for n_b in range( 0, len( sp_obj_b.peak_list_mz_int_abs)):
                                # if match, count
                                #print "---pk_b" ,  sp_obj_b.peak_list_mz_int_abs[n_b][0]
                                flag_pk_mz_matched = 0

                                if (   math.fabs(   (    (  sp_obj_a.peak_list_mz_int_abs[n_a][0] - sp_obj_b.peak_list_mz_int_abs[n_b][0]  )/ sp_obj_b.peak_list_mz_int_abs[n_b][0]  )  * 1000000  )   < mass_tol_ppm ):
                                    flag_pk_mz_matched = 1
                                    #print "ppm match"
                                if  (sp_obj_a.peak_list_mz_int_abs[n_a][0] - sp_obj_b.peak_list_mz_int_abs[n_b][0] == 0 ):
                                    #print "zero match"
                                    flag_pk_mz_matched = 1

                                if(flag_pk_mz_matched):
                                    count_peak_match = count_peak_match +1


                                # if current theo frag is larger than exp frag, break.

                                if (    sp_obj_a.peak_list_mz_int_abs[n_a][0]   <   sp_obj_b.peak_list_mz_int_abs[n_b][0]     + 2   )  :
                                    break



                        #print count_peak_match
                        # if more than X % of peaks in A match match to B   and number of peak of B is smaller than A,   B is not necessary (redundant)....

                        if (  float( count_peak_match ) /  float( len(sp_obj_b.peak_list_mz_int_abs))  >   ratio_pk_b_explained    and   len(sp_obj_a.peak_list_mz_int_abs)  >  len(sp_obj_b.peak_list_mz_int_abs)) :
                            #print idx_b , "added to kill"
                            list_idx_tb_deleted.append(idx_b )
                        # if number of peak is same, just prefer earlier one
                        if (  float( count_peak_match ) /  float( len(sp_obj_b.peak_list_mz_int_abs))  >   ratio_pk_b_explained    and   len(sp_obj_a.peak_list_mz_int_abs) ==  len(sp_obj_b.peak_list_mz_int_abs)) :
                            #print idx_b , "added to kill"
                            if(idx_a < idx_b):
                                list_idx_tb_deleted.append(idx_b )

    # make  "list_idx_tb_deleted" non redundant
    print("list_idx_tb_deleted" , list_idx_tb_deleted)
    list_idx_tb_deleted = list(set(list_idx_tb_deleted))


    list_spectrum_uni_object_nr = []

    for idx in range( 0,  len(list_spectrum_uni_object)) :

        flag_delete = 0
        if idx in list_idx_tb_deleted:
            flag_delete = 1

        if flag_delete == 0:

            list_spectrum_uni_object_nr.append(list_spectrum_uni_object[idx])

    return list_spectrum_uni_object_nr


#######################################
import sys
def get_top_x_fragment_rich_spectra( list_spectrum_uni_object , topx ):

    list_topx_spectra = []
    count = 0

    # making list of "number of peaks"
    list_len_peaklist = []
    for sp in list_spectrum_uni_object:
        list_len_peaklist.append( len(sp.peak_list_mz_int_abs) )

    #print "before" , list_len_peaklist


    list_len_peaklist.sort()
    list_len_peaklist.reverse()

    #print "after" , list_len_peaklist

    #print "len(list_len_peaklist)", len(list_len_peaklist)

    #print "topx"
    num_of_frag_of_topx = 0


    if(  len(list_len_peaklist) >=  topx   ) :
        num_of_frag_of_topx  = list_len_peaklist[ topx - 1 ]
        """
        try:
            num_of_frag_of_topx  = list_len_peaklist[ topx - 1 ]
        except:
            print "ERROR from spectrum_uni get_top_x_fragment_rich_spectra , top x is smaller than the length of list of spectra ?"
            sys.exit()
        """

        #print "num_of_frag_of_topx" , num_of_frag_of_topx
        for sp in list_spectrum_uni_object:

            #print mgf_spec.title   , len(mgf_spec.peak_list)
            if (len(sp.peak_list_mz_int_abs)  >= num_of_frag_of_topx ):
                print(len(sp.peak_list_mz_int_abs)  , num_of_frag_of_topx)
                list_topx_spectra.append( sp  )
                #print "to append"
    else :
        list_topx_spectra = list_spectrum_uni_object

    return list_topx_spectra

#
def get_spectra_with_precursor_type( list_spectrum_uni_object , precursor_type    ):


    list_spectra_prec_x = []


    for sp in list_spectrum_uni_object:

        if(sp.precursor_type == precursor_type):
            list_spectra_prec_x.append( sp)

    return list_spectra_prec_x
#####################
#  compound_uni related functions
#########################
# this function read list of spectrum uni objects and return a list of compound uni objects which may have multiple spectra for one compound.
# In other words, you can make list of compound together with bundled ms2 spectra,

def generate_list_of_compound_uni_with_spec_uni_set( list_spectrum_uni_object):

    list_inchikey = []

    for spec_o in list_spectrum_uni_object :
        #print spec_o.accession_number
        #print "spec_o.inchi_key:->" , spec_o.inchi_key , "+++"
        if( len(spec_o.inchi_key) > 2):
            list_inchikey.append(spec_o.inchi_key)

    # make non-redundant--------------
    list_inchikey_nr = list(set(list_inchikey))

    # make list of empty compound_uni object just having inchikey, which is now non-redundant
    list_compound_uni_obj = []

    for inchikey_nr in list_inchikey_nr :
        #compound_uni_obj = spectrum_uni_a1.compound_uni_class()
        compound_uni_obj = compound_uni_class()

        compound_uni_obj.inchi_key = inchikey_nr
        list_compound_uni_obj.append( compound_uni_obj   )


    # iterate spectrum_uni object list and register to compound object based on inchikeyz
    count_spectra_match = 0
    for compound_obj in list_compound_uni_obj:

        for spec_o in list_spectrum_uni_object :

            if( len(spec_o.inchi_key) > 2    and  len(compound_obj.inchi_key) > 2 ) :

                if(    spec_o.inchi_key == compound_obj.inchi_key          ):

                    #print  spec_o.inchi_key , spec_o.inchi_key
                    compound_obj.list_spectrum_uni.append(spec_o )
                    count_spectra_match = count_spectra_match +1
                    compound_obj.name_represen = spec_o.name
                    compound_obj.inchi = spec_o.inchi
                    compound_obj.smiles = spec_o.smiles

    #print count_spectra_match
    #print "len(list_compound_uni_obj)" ,len(list_compound_uni_obj)

    return list_compound_uni_obj

### this function returns list of spectrum, but 1 spectrum for one compound (inchikey)
# return list_represen_spectrum_uni
# mode 1 :  highest number    2:  lowest number of   3:  middle number   of peaks as representative
def generate_list_of_represen_spectrum_for_unique_compound(list_compound_uni_object , mode):

    list_represen_spectrum_uni = []
    import operator

    for cmpd_o in list_compound_uni_object:

        #for spec_o in cmpd_o.list_spectrum_uni :
        cmpd_o.list_spectrum_uni.sort(key=operator.attrgetter('num_peaks'))

        #for o in cmpd_o.list_spectrum_uni:
        #    print o.num_peaks


        #   choose spectrum with highest number of peaks as representative
        if(mode == 1) :
            spec_to_choose = cmpd_o.list_spectrum_uni[    len(cmpd_o.list_spectrum_uni) -1   ]
            list_represen_spectrum_uni.append(spec_to_choose)
            #print "chosen:" , spec_to_choose.num_peaks

        #   choose spectrum with lowest number of peaks as representative
        if(mode == 2) :
            spec_to_choose = cmpd_o.list_spectrum_uni[    0  ]
            list_represen_spectrum_uni.append(spec_to_choose)
            #print "chosen:" , spec_to_choose.num_peaks


        #   choose spectrum with middle number of peaks as representative
        if(mode == 3) :
            middle_pos =   int ( len(cmpd_o.list_spectrum_uni) / 2 )
            #print "middle pos is " ,  middle_pos

            spec_to_choose = cmpd_o.list_spectrum_uni[   middle_pos  ]
            list_represen_spectrum_uni.append(spec_to_choose)
            #print "chosen:" , spec_to_choose.num_peaks


    return list_represen_spectrum_uni


### modify


def introduce_mass_error_to_list_spec_uni_obj(list_spec_uni_obj , error_ppm):

    for obj_sp in list_spec_uni_obj:

        pkl_abs_mod =  edit_mass_spectrum_a1.introduce_mass_error_to_peak_list_by_ppm( obj_sp.peak_list_mz_int_abs , error_ppm)

        obj_sp.peak_list_mz_int_abs = pkl_abs_mod
        pkl_rel_mod =  edit_mass_spectrum_a1.introduce_mass_error_to_peak_list_by_ppm( obj_sp.peak_list_mz_int_rel , error_ppm)

        obj_sp.peak_list_mz_int_rel = pkl_rel_mod

    return list_spec_uni_obj




#############################
#  introduce random mass shift (delta) in int
#   you can specigy range of delta value  "start" , "end".   if start = 0  and end is 50, rondom delta is integer between 0 to 50
###############################################

import random

def introduce_random_delta_to_spec_uni(spec_uni, start , end):

    for n in range(0, len(spec_uni.peak_list_mz_int_abs)):

        rand_delta = random.randint(start , end)
        spec_uni.peak_list_mz_int_abs[n][0] =  spec_uni.peak_list_mz_int_abs[n][0] + rand_delta
        spec_uni.peak_list_mz_int_rel[n][0] =  spec_uni.peak_list_mz_int_rel[n][0] + rand_delta


    spec_uni.peak_list_mz_int_abs = sorted(spec_uni.peak_list_mz_int_abs)
    spec_uni.peak_list_mz_int_rel = sorted(spec_uni.peak_list_mz_int_rel)

def normalize_peak_list_rel_to_100(spec_uni):

    spec_uni.peak_list_mz_int_rel.sort(key=lambda e:e[1],reverse = True)

    highest_int = spec_uni.peak_list_mz_int_rel[0][1]
    # sort again by mz
    spec_uni.peak_list_mz_int_rel.sort(key=lambda e:e[0])

    for pk in spec_uni.peak_list_mz_int_rel:

        pk[1] =   float(pk[1]) / float (highest_int)


def get_topN_most_intense_peak_list_from_peaklist( list_peak, top_n):

    list_peak.sort(key=lambda e:e[1], reverse=True)
    list_temp =[]
    count_n = 0

    #iterate list sorted with 3rd value
    # til reaches top_n
    for obj in list_peak:
        list_temp.append(obj)
        count_n = count_n + 1
        if( count_n == top_n):break
    list_temp.sort(key=lambda e:e[0])

    list_peak.sort(key=lambda e:e[0])

    return list_temp


# this function returns peaklist with absolute intensity of TopN most intense peak from spec uni object
def get_topN_most_intense_peaklist_abs_from_spec_uni( spec_uni , top_n):

    return get_topN_most_intense_peak_list_from_peaklist( spec_uni.peak_list_mz_int_abs , top_n    )

# this function returns peaklist with relative intensity of TopN most intense peak  from spec uni object

def get_topN_most_intense_peaklist_rel_from_spec_uni( spec_uni , top_n):

    return get_topN_most_intense_peak_list_from_peaklist( spec_uni.peak_list_mz_int_rel , top_n    )




# Eaually pickup intesne peaks over segmented ranges

def get_topN_most_intense_peak_list_in_binned_rangesW( list_peak, top_n , bin_range_w):

    list_topN_peaklist_in_binned_rangesW = []

    if( len( list_peak) > 1):

        highest_mz = list_peak[len(list_peak) -1] [0]

        n = 0

        list_mz_segment = []
        while n * bin_range_w <= highest_mz :
            #print n , n + 1
            #print  n * float(bin_range_w)  ,"to", (n+1) * float(bin_range_w)
            list_mz_segment.append( n * float(bin_range_w) )
            n = n + 1

        list_mz_segment.append(highest_mz)

        #print list_mz_segment



        # iterate range (  like   0 to 50mz, 50 to 100mz )
        for n in range(  1, len(list_mz_segment)):

            #print "range" ,  list_mz_segment[n-1] , list_mz_segment[n]

            list_pk_in_bin = []

            # scan input peak list
            for pk in list_peak:

                # pk[0] is m/z value
                if (   list_mz_segment[n-1]  <   pk[0]       ) and (    pk[0] <=  list_mz_segment[n]     ) :

                    #print "hit", pk
                    list_pk_in_bin.append( pk )

            # pickup top N in this range-----------

            # reverse sort by intensity
            list_pk_in_bin.sort(key=lambda e:e[1], reverse=True)

            list_temp =[]
            count_n = 0
            #iterate list sorted with 3rd value
            # til reaches top_n
            for obj in list_pk_in_bin:
                list_temp.append(obj)
                count_n = count_n + 1
                if( count_n == top_n):break

            # sort by m/z
            list_temp.sort(key=lambda e:e[0])
            list_topN_peaklist_in_binned_rangesW  = list_topN_peaklist_in_binned_rangesW  + list_temp

    if( len( list_peak) <= 1):
        list_topN_peaklist_in_binned_rangesW = list_peak

    return list_topN_peaklist_in_binned_rangesW



def set_topN_most_intense_peak_list_in_binned_rangesW( spec_uni_obj, top_n , bin_range_w, max_value_rel_int = 100 ):

    spec_uni_obj.peak_list_mz_int_abs = get_topN_most_intense_peak_list_in_binned_rangesW( spec_uni_obj.peak_list_mz_int_abs, top_n , bin_range_w)
    spec_uni_obj.peak_list_mz_int_rel = get_topN_most_intense_peak_list_in_binned_rangesW( spec_uni_obj.peak_list_mz_int_rel, top_n , bin_range_w)
    set_peak_list_rel_int(spec_uni_obj,max_value_rel_int)


# Note this function is using deepcopy !!

def clean_common_contami_peaks_author_inst(list_spec_uni , mz_tol ,rate,min_hit_to_remove):

    list_new_cleaned_spec_uni =[]
    # iterate spec
    for spec_uni in list_spec_uni :
        # set evert\ything

        # this variable count how many spectra are from a combination of "author" x "instrument"
        #  if there are 100 spectra from the combination of "Hayakawa"  "Q-exactive", in this entire dataset, count_same_author_inst =  100
        count_same_author_inst =  0
        list_count_hit = [0] * len(spec_uni.peak_list_mz_int_abs )

        # do comparison for all the other
        for spec_uni_b in list_spec_uni :

            # if the authors and instrument tag is same
            if(  spec_uni.authors ==   spec_uni_b.authors      and        spec_uni.instrument ==   spec_uni_b.instrument   ):

                count_same_author_inst  = count_same_author_inst  +1

                # do comparison
                for n in range( 0,   len(spec_uni.peak_list_mz_int_abs)):

                    for pk_b in spec_uni_b.peak_list_mz_int_abs:

                        if ( abs(  spec_uni.peak_list_mz_int_abs[n][0] - pk_b[0])  < mz_tol):

                            list_count_hit[n] = list_count_hit[n]+1
                        if  ( pk_b[0]    >    spec_uni.peak_list_mz_int_abs[n][0] +  mz_tol ):
                            break

        # comparison for this spectrum (spec_uni) to all the other spec is done.
        # now  "count_same_author_inst" indicates the number of spectra from same "author", "instrument" setting as "spec_uni"
        #   "list_count_hit[n]" now holds which peak match how many times to spectra from same "author", "instrument" setting as "spec_uni"
        # now you decide which peaks to remove


        spec_uni_cleaned = copy.deepcopy(spec_uni)

        new_peak_list_mz_int_abs =[]
        new_peak_list_mz_int_rel =[]


        print(spec_uni_cleaned.peak_list_mz_int_abs)
        for n in range (0, len(list_count_hit)):
            print(n)
            flag_to_remove = 0
            # if the current peak is found more than twice in spectra( from same author, same instrument)
            if (   list_count_hit[n] > 0   ):
                if (     float(list_count_hit[n]) / float(count_same_author_inst) > rate      and  list_count_hit[n]  > min_hit_to_remove   ):
                    flag_to_remove = 1
                    #del spec_uni_cleaned.peak_list_mz[n]
                    #del spec_uni_cleaned.peak_list_int[n]
                    #del spec_uni_cleaned.peak_list_mz_int_abs[n]
                    #del spec_uni_cleaned.peak_list_mz_int_rel[n]

            if (flag_to_remove == 0 ) :
                new_peak_list_mz_int_abs.append(spec_uni.peak_list_mz_int_abs[n])
                new_peak_list_mz_int_rel.append(spec_uni.peak_list_mz_int_rel[n])


        # cleaned peak list is ready. Replace with old ones
        spec_uni_cleaned.peak_list_mz_int_abs =  new_peak_list_mz_int_abs
        spec_uni_cleaned.peak_list_mz_int_rel =  new_peak_list_mz_int_rel

        list_new_cleaned_spec_uni.append( spec_uni_cleaned )

    return list_new_cleaned_spec_uni






##############
#  Export spectrum
################

# this function is to write and export soectrum data to massbank-like format
def write_spectrum_uni_to_msbk_file(spec_uni_obj, folder_to_output ):

    #folder_to_output = "msbk_file2"
    if not os.path.exists(folder_to_output):
        os.makedirs(folder_to_output)


    # new filename
    ft2_filename =   (spec_uni_obj.source_filename).split(".")[0] + ".txt"
    # file stream


    fo_msbk = open( folder_to_output + "/"+ ft2_filename, 'w')

    # content to export as msbk
    str_msbk = "ACCESSION: " + str(spec_uni_obj.accession_number) + "\n"
    str_msbk = str_msbk  + "CH$NAME: " + str(spec_uni_obj.name) + "\n"
    str_msbk = str_msbk  + "CH$EXACT_MASS: " + str(spec_uni_obj.exact_mass)  + "\n"
    if  ( len(spec_uni_obj.smiles) ) :
        str_msbk = str_msbk  + "CH$SMILES: " + str(spec_uni_obj.smiles)  + "\n"

    if  ( len(spec_uni_obj.inchi) > 2 ) :
        str_msbk = str_msbk  + "CH$IUPAC: " + str(spec_uni_obj.inchi)  + "\n"
    if  ( len(spec_uni_obj.inchi_key) > 2 ) :
        str_msbk = str_msbk  + "CH$LINK: INCHIKEY " + str(spec_uni_obj.inchi_key)  + "\n"

    if  ( len(spec_uni_obj.spectrum_type) > 2 ) :
        str_msbk = str_msbk  + "AC$INSTRUMENT_TYPE: " + str(spec_uni_obj.instrument_type)  + "\n"

    str_msbk = str_msbk  + "AC$MASS_SPECTROMETRY: MS_TYPE " + str(spec_uni_obj.spectrum_type)  + "\n"
    str_msbk = str_msbk  + "AC$MASS_SPECTROMETRY: IONIZATION " + str(spec_uni_obj.spectrum_type)  + "\n"
    str_msbk = str_msbk  + "AC$MASS_SPECTROMETRY: FRAGMENTATION_MODE " + str(spec_uni_obj.fragmentation_type)  + "\n"
    str_msbk = str_msbk  + "AC$MASS_SPECTROMETRY: ION_MODE " + str(spec_uni_obj.ionization_mode)  + "\n"

    if ( spec_uni_obj.precursor_mz > 0):
        str_msbk = str_msbk  + "MS$FOCUSED_ION: PRECURSOR_M/Z " + str(spec_uni_obj.precursor_mz)  + "\n"

    if ( len(spec_uni_obj.precursor_type ) > 0):
        str_msbk = str_msbk  + "MS$FOCUSED_ION: PRECURSOR_TYPE " + str(spec_uni_obj.precursor_type )  + "\n"


    str_msbk = str_msbk  + "PK$NUM_PEAK: " + str( len(spec_uni_obj.peak_list_mz_int_abs) )  + "\n"
    str_msbk = str_msbk  + "PK$PEAK: m/z int. rel.int." +"\n"


    print("len(spec_uni_obj.peak_list_mz_int_abs)" , len(spec_uni_obj.peak_list_mz_int_abs))
    for n in range(0 , len( spec_uni_obj.peak_list_mz_int_abs )):

        print(spec_uni_obj.peak_list_mz_int_abs[n][0] , spec_uni_obj.peak_list_mz_int_abs[n][1] , spec_uni_obj.peak_list_mz_int_rel[n][1])

        str_msbk = str_msbk  +  "  "  + str(spec_uni_obj.peak_list_mz_int_abs[n][0])   +     " "  + str(spec_uni_obj.peak_list_mz_int_abs[n][1] ) +  " "  + str( spec_uni_obj.peak_list_mz_int_rel[n][1] )+ "\n"
    str_msbk = str_msbk  + "//" + "\n" + "\n"

    fo_msbk.write( str_msbk)


# this function write spectrum to massbank format, !!!! but after editing,,,
# Since you might have edited peak info original read from mass bank, now that you can't just simply write down info to
# for instance, peak list and number of peaks are no longer consistent (unless you managed to coordinate number of peaks and peak list, but it might be difficult for all of them)

# TO DO :  implement peak annotation
def write_spectrum_uni_to_msbk_file_after_edit(spec_uni_obj, folder_to_output ):

    #folder_to_output = "msbk_file2"
    if not os.path.exists(folder_to_output):
        os.makedirs(folder_to_output)

    # new filename
    output_filename =   (spec_uni_obj.source_filename).split(".")[0] + ".txt"

    fo_msbk = open( folder_to_output + "/"+ output_filename, 'w')

    str_msbk = spec_uni_obj.source_txt_data_till_peak_info

    str_msbk = str_msbk + "PK$NUM_PEAK: " +  str( len(spec_uni_obj.peak_list_mz_int_abs)) + "\n"
    str_msbk = str_msbk + "PK$PEAK: m/z int. rel.int.\n"

    if(  len(spec_uni_obj.peak_list_mz_int_abs)  < 1):
        print("ALERT: your peak list has no peak !!!!")

    # iterate peak list (potentially modified externally !!) and add to text
    for n in range(0 , len(spec_uni_obj.peak_list_mz_int_abs)):
        str_msbk = str_msbk + "  " +   str(spec_uni_obj.peak_list_mz_int_abs[n][0])          + " " +    str(spec_uni_obj.peak_list_mz_int_abs[n][1])      + " " +    str(spec_uni_obj.peak_list_mz_int_rel[n][1])                         + "\n"

    # add end_of_file mark
    str_msbk = str_msbk + "//\n"

    # write to file
    fo_msbk.write( str_msbk)





def convert_spectrum_uni_to_ft2(spec_uni_obj):

    # new filename
    ft2_filename =   (spec_uni_obj.source_filename).split(".")[0] + ".FT2"
    print("ft2",ft2_filename)

    fo_ft2 = open(ft2_filename, 'w')
    #   S   line
    fo_ft2.write( "S\t" +   str(spec_uni_obj.scan_number)  + "\t" +   str(spec_uni_obj.scan_number)  +  "\t" + str(spec_uni_obj.precursor_mz)     + "\n")
    #   Z   line
    fo_ft2.write( "I	RetentionTime	0\n")

    for n in range( 0 , len(spec_uni_obj.peak_list_mz_int_abs)):
        fo_ft2.write( " "+  str(spec_uni_obj.peak_list_mz_int_abs[0] ) + " " + str( spec_uni_obj.peak_list_mz_int_abs[1] ) + "\t0\t0\t0\t0\n"  )

    '''
    fo_ft2 = open('converted.ft2', 'w')
    for obj in list_msbk_cmpd_obj:
        #   S   line
        fo_ft2.write( "S\t" +   obj.accession  + "\t" +   obj.accession  +  "\t" + str(obj.prec_mz)     + "\n")
        #   Z   line
        fo_ft2.write( "Z\t" +   obj.charge_str  + "\t" +    str(obj.prec_mz)     + "\n")

        for pk in  obj.peak_list:
            fo_ft2.write( " "+  str(pk[0] ) + " " + str( pk[1] ) + "\t0\t0\t0\t0\n"  )

        #print obj.accession
        #print obj.peak_list
        continue
    '''
    '''
        current_FT2_file.write("S\t"+s_scanId+"\t"+s_scanId+"\t"+str(d_precursor_mz)+"\n")
    if (i_precursor_z == 0) :
        z_mz = 0
    else :
        z_mz = i_precursor_z*d_precursor_mz
    current_FT2_file.write("Z\t"+"-"+str(i_precursor_z)+"\t"+str(z_mz)+"\n")
    current_FT2_file.write("I\tRetentionTime\t"+str(d_retention_time)+"\n")
    for each_peak in current_spectrum_list :
        current_FT2_file.write(each_peak[0]+"\t"+each_peak[1]+"\t0\t0\t0\t0\n")
    '''

# Note this function is using deepcopy !!

def write_summary_table(list_spec_uni , table_file_name ):

    fo = open(table_file_name , 'w')

    str_out =  "ACCESSION" + "\t" + "NAME" + "\t" +  "InChi" + "\t" +   "SMILES" + "\n"

    # iterate spec
    for spec_uni in list_spec_uni :

        str_out =  spec_uni.accession_number + "\t" + spec_uni.name + "\t" +  spec_uni.inchi + "\t" +   spec_uni.smiles + "\n"
        fo.write(str_out)


    fo.close()

