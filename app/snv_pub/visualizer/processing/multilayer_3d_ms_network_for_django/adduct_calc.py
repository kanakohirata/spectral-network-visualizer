dic_adduct_vs_dic_calc = { \
    "M+3H": {"multi": 0.33, "mass": 1.007276}, \
    "M+2H+Na": {"multi": 0.33, "mass": 8.33459}	, \
    "M+H+2Na" :	{	"multi":	0.33	,	"mass":	15.76619  }		, \
    "M+3Na" :	{	"multi":	0.33	,	"mass":	22.989218  }		, \
    "M+2H" :	{	"multi":	0.5	,	"mass":	1.007276  }		, \
    "M+H+NH4" :	{	"multi":	0.5	,	"mass":	9.52055  }		, \
    "M+H+Na" :	{	"multi":	0.5	,	"mass":	11.998247  }		, \
    "M+H+K" :	{	"multi":	0.5	,	"mass":	19.985217  }		, \
    "M+ACN+2H" :	{	"multi":	0.5	,	"mass":	21.52055  }		, \
    "M+2Na" :	{	"multi":	0.5	,	"mass":	22.989218  }		, \
    "M+2ACN+2H" :	{	"multi":	0.5	,	"mass":	42.033823  }		, \
    "M+3ACN+2H" :	{	"multi":	0.5	,	"mass":	62.547097  }		, \
    "M+H" :	{	"multi":	1	,	"mass":	1.007276  }		, \
    "M+NH4" :	{	"multi":	1	,	"mass":	18.033823  }		, \
    "M+Na" :	{	"multi":	1	,	"mass":	22.989218  }		, \
    "M+CH3OH+H" :	{	"multi":	1	,	"mass":	33.033489	  }	, \
    "M+K" :	{	"multi":	1	,	"mass":	38.963158	  }	, \
    "M+ACN+H" :	{	"multi":	1	,	"mass":	42.033823	  }	, \
    "M+2Na-H" :	{	"multi":	1	,	"mass":	44.97116	  }	, \
    "M+IsoProp+H" :	{	"multi":	1	,	"mass":	61.06534	  }	, \
    "M+ACN+Na" :	{	"multi":	1	,	"mass":	64.015765	  }	, \
    "M+2K-H" :	{	"multi":	1	,	"mass":	76.91904	  }	, \
    "M+DMSO+H" :	{	"multi":	1	,	"mass":	79.02122	  }	, \
    "M+2ACN+H" :	{	"multi":	1	,	"mass":	83.06037	  }	, \
    "M+IsoProp+Na+H" :	{	"multi":	1	,	"mass":	84.05511	  }	, \
    "2M+H" :	{	"multi":	2	,	"mass":	1.007276	  }	, \
    "2M+NH4" :	{	"multi":	2	,	"mass":	18.033823	  }	, \
    "2M+Na" :	{	"multi":	2	,	"mass":	22.989218	  }	, \
    "2M+K" :	{	"multi":	2	,	"mass":	38.963158	  }	, \
    "2M+ACN+H" :	{	"multi":	2	,	"mass":	42.033823	  }	, \
    "2M+ACN+Na" :	{	"multi":	2	,	"mass":	64.015765	  }	, \
    "M-3H" :	{	"multi":	0.33	,	"mass":	-1.007276	  }	, \
    "M-2H" :	{	"multi":	0.5	,	"mass":	-1.007276	  }	, \
    "M-H2O-H" :	{	"multi":	1	,	"mass":	-19.01839	  }	, \
    "M-H" :	{	"multi":	1	,	"mass":	-1.007276	  }	, \
    "M+Na-2H" :	{	"multi":	1	,	"mass":	20.974666	  }	, \
    "M+Cl" :	{	"multi":	1	,	"mass":	34.969402	  }	, \
    "M+K-2H" :	{	"multi":	1	,	"mass":	36.948606	  }	, \
    "M+FA-H" :	{	"multi":	1	,	"mass":	44.998201	  }	, \
    "M+Hac-H" :	{	"multi":	1	,	"mass":	59.013851	  }	, \
    "M+Br" :	{	"multi":	1	,	"mass":	78.918885	  }	, \
    "M+TFA-H" :	{	"multi":	1	,	"mass":	112.985586	  }	, \
    "2M-H" :	{	"multi":	2	,	"mass":	-1.007276	  }	, \
    "2M+FA-H" :	{	"multi":	2	,	"mass":	44.998201	  }	, \
    "2M+Hac-H" :	{	"multi":	2	,	"mass":	59.013851	  }	, \
    "3M-H" :	{	"multi":	3	,	"mass":	-1.007276	  } \
    }


def get_mz_based_on_exactmw_and_adduct( exact_mw , str_adduct   ):


    if str_adduct in dic_adduct_vs_dic_calc:
        dic_calc = dic_adduct_vs_dic_calc[str_adduct]
        # do calc
        mz = exact_mw * dic_calc["multi"] + dic_calc["mass"]
    else :
        mz = 0
        print(" [adduct_calc] input aduct string "   + str_adduct  +  " is not available.")
    return mz



mz =  get_mz_based_on_exactmw_and_adduct( 123.01, "M-H"  )
print ("mz:" , mz)