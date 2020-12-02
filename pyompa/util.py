from __future__ import division, print_function

#Remineralization ratios for computing PO, NO and SiO values
R_PO = 155
R_SiO = 15
R_NO = 9.68


def augment_df_with_PO_NO_SiO(df): 
    """
    Augments a data frame with the values for PO, NO and SiO. Assumes there are
    columns 'phosphate', 'nitrate', 'silicate' and 'oxygen' in the data frame.

    PO = phosphate*155 + oxygen
    NO = nitrate*9.68 + oxygen 
    SiO = silicate*15 + oxygen

    Args:
        df: the pandas  data frame
    """
    df["PO"] = df["oxygen"] + df["phosphate"]*R_PO
    df["NO"] = df["oxygen"] + df["nitrate"]*R_NO
    df["SiO"] = df["oxygen"] + df["silicate"]*R_SiO

