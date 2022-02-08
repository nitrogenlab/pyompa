from __future__ import division, print_function
from collections import OrderedDict
import numpy as np


def assert_compatible_keys(the_dict, allowed, errorprefix):
    for key in the_dict:
        assert_in(value=key, allowed=allowed, errorprefix=errorprefix)


def assert_in(value, allowed, errorprefix):
    assert value in allowed,\
        errorprefix+str(value)+" not in allowed values: "+str(allowed)


def assert_has_keys(the_dict, required, errorprefix):
    for key in required:
        assert key in the_dict, (
            errorprefix+key+" must be specified. Provided keys are: "
            +str(the_dict.keys()))


def get_endmember_idx_mapping(endmember_names):
    """
    Get a mapping from endmember name to indexes with all the subtypes
     corresponding to that endmember. The subtypes should be denoted with
     _subtype.
    """
    endmembername_to_indices = OrderedDict()

    for idx,endmember_name in enumerate(endmember_names):
        if "_" in endmember_name:
            endmembername_core = endmember_name.split("_")[0]
        else:
            endmembername_core = endmember_name
        if endmembername_core not in endmembername_to_indices:
            endmembername_to_indices[endmembername_core] = [] 
        endmembername_to_indices[endmembername_core].append(idx) 

    return endmembername_to_indices 


def collapse_endmembers_by_idxmapping(endmember_fractions,
                                      endmembername_to_indices):
    endmember_names = endmembername_to_indices.keys()
    remapped_endmember_fractions = np.zeros(
        (len(endmember_fractions), len(endmember_names))) 
    for remapped_idx,endmember_name in enumerate(endmember_names):
        idxs = endmembername_to_indices[endmember_name]
        for idx in idxs:
            remapped_endmember_fractions[:,remapped_idx] +=\
                endmember_fractions[:,idx] 
    return remapped_endmember_fractions


def organize_converted_vars_by_groupname(converted_variables,
                                         convertedparam_groups):
    """
    converted_variables: a matrix of dimensions
        num_converted_variables X num_parameters
    """
    groupname_to_totalconvertedvariable = OrderedDict()
    groupname_to_effectiveconversionratios = OrderedDict()
    convar_idx = 0
    for convertedparam_group in convertedparam_groups:
        groupname = convertedparam_group.groupname
        #each type of converted variable may have multiple instances of
        # itself corresponding to different conversion ratios
        convar_vals = converted_variables[:,
            convar_idx:
            convar_idx+len(convertedparam_group.conversion_ratios)]
        convar_idx += len(convertedparam_group.conversion_ratios)
        #sanity check the signs; for each entry they
        # should either be all positive or all negative, within numerical
        # precision
        for convar_vals_row in convar_vals:
            if ((all(convar_vals_row >= 0) or
                 all(convar_vals_row <= 0))==False):
                print("WARNING: sign inconsistency in "
                      +groupname+":", convar_vals_row)
        total_convar = np.sum(convar_vals, axis=-1)
        groupname_to_totalconvertedvariable[groupname] =\
            total_convar 
        #proportions of the converted variable used at differnet ratios
        with np.errstate(divide='ignore', invalid='ignore'):
            convarusage_proportions = (
                convar_vals/total_convar[:,None])
            effective_conversion_ratios = OrderedDict()
            conversion_ratios_dict =\
                convertedparam_group.get_conversion_ratios_dict() 
            effective_conversion_ratios = OrderedDict([
                (relevant_param_name,
                 convarusage_proportions@conversion_ratio)
                for relevant_param_name,conversion_ratio in
                conversion_ratios_dict.items()
             ]) 
            groupname_to_effectiveconversionratios[groupname] =\
                effective_conversion_ratios

    return (groupname_to_totalconvertedvariable,
            groupname_to_effectiveconversionratios) 

#import gsw
#
#
##Remineralization ratios for computing PO, NO and SiO values
#R_PO = 155
#R_SiO = 15
#R_NO = 9.68
#
#
#def augment_df_with_PO_NO_SiO(df): 
#    """
#    Augments a data frame with the values for PO, NO and SiO. Assumes there are
#    columns 'phosphate', 'nitrate', 'silicate' and 'oxygen' in the data frame.
#
#    PO = phosphate*155 + oxygen
#    NO = nitrate*9.68 + oxygen 
#    SiO = silicate*15 + oxygen
#
#    Args:
#        df: the pandas  data frame
#    """
#    df["PO"] = df["oxygen"] + df["phosphate"]*R_PO
#    df["NO"] = df["oxygen"] + df["nitrate"]*R_NO
#    df["SiO"] = df["oxygen"] + df["silicate"]*R_SiO
#
#
#def prepare_water_mass_df(water_mass_arr):
#    df = pd.DataFrame(data=water_mass_arr,
#                      columns=["watermassname",
#                               "potential_temp", "practical_salinity",
#                               "oxygen", "phosphate", "silicate",
#                               "nitrate", "spiciness", "PV"])
#    augment_df_with_PO_NO_SiO(df)
#    return df
#
#
#def prepare_observations(input_file):
#    """
#    This code is specific to Rian's data format
#    """
#    header = ["c"+str(i) for i in range(1,30)]
#    header[4] = "bottle flag"
#    header[16] = "bottle salinity flag"
#    header[20] = "bottle oxygen flag"
#    header[22] = "silicate flag"
#    header[24] = "nitrate flag"
#    header[28] = "phosphate flag"
#
#    header[11] = "CTD pressure"
#    header[12] = "CTD temperature"
#    header[15] = "practical_salinity"
#    header[8] = "latitude"
#    header[9] = "longitude"
#
#    header[0] = "stnnbr"
#    header[5] = "geotrc_ID"
#    header[10] = "bottom depth"
#    header[19] = "oxygen"
#    header[21] = "silicate"
#    header[23] = "nitrate"
#    header[27] = "phosphate"
#
#    gp15_df = pd.read_csv(
#        "names_added_GP15OMPA_33RR20180918_only_gs_rosette_clean1_hy1.csv",
#         names=header, na_values = -999)
#
#    #remove bad data
#    for flag_type in ["bottle flag", "bottle salinity flag", "bottle oxygen flag",
#                  "silicate flag", "nitrate flag", "phosphate flag"]:
#        gp15_df = gp15_df[gp15_df[flag_type] <= 2]
#    gp15_df = pd.DataFrame(gp15_df)
#
#    augment_df_with_PO_NO_SiO(gp15_df)
#
#    #absolute salinty
#    absolute_salinity = gsw.SA_from_SP(SP=gp15_df["practical_salinity"],
#                                       p=gp15_df["CTD pressure"],
#                                       lon=gp15_df["longitude"],
#                                       lat=gp15_df["latitude"]) 
#    gp15_df["absolute_salinity"] = absolute_salinity
#    
#    #conservative temp
#    conservative_temp = gsw.CT_from_t(SA=absolute_salinity,
#                                  t=gp15_df["CTD temperature"],
#                                  p=gp15_df["CTD pressure"])
#    gp15_df["conservative_temp"] = conservative_temp
#
#    #potential temp
#    potential_temp = gsw.pt_from_CT(SA=absolute_salinity,
#                                CT=conservative_temp)
#    gp15_df["potential_temp"] = potential_temp
#
#    #sig0
#    sig0 = gsw.rho(SA=absolute_salinity, CT=conservative_temp, p=0) - 1000
#    gp15_df["sig0"] = sig0
#
#    z = gsw.z_from_p(p=gp15_df["CTD pressure"], lat=gp15_df["latitude"])
#    depth = -z #https://github.com/TEOS-10/python-gsw/blob/7d6ebe8114c5d8b4a64268d36100a70e226afaf6/gsw/gibbs/conversions.py#L577
#    gp15_df["depth"] = depth
#
#    spic0 = gsw.spiciness0(SA=absolute_salinity, CT=conservative_temp)
#    gp15_df["spiciness"] = spic0
#
#    #calculation of planetary vorticity
#    rho_ref = 1000.0 #reference density
#    Omega = 2*np.pi/86400;
#    f = np.mean(2*Omega*np.sin(gp15_df["latitude"]*(np.pi/180.0)))
#    PV = -(f/rho_ref)*(np.gradient(sig0, depth))
#    gp15_df["PV"] = PV #potential vorticity
#
#    print("Rows:",len(gp15_df))
#    gp15_df = gp15_df.dropna()
#    print("Rows without NA values:",len(gp15_df))
