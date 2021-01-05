from __future__ import division, print_function
from .ompacore import OMPAProblem
import numpy as np
import pandas as pd


def get_endmember_df_for_sig0_range(endmemnames_to_use,
                                    endmemname_to_df,
                                    sig0_bin_start,
                                    sig0_bin_end):
    #Idea: for each entry in endmemname_to_df, filter out the row
    # that lies within sig0_bin_start and sig0_bin_end (there should
    # only be one such row). Then assemble all those rows into a data frame
    # of the endmembers to use for the OMPA analysis, for that sig0 range

    #We'll create a list of the correct rows from each pandas df
    correct_rows = []
    sig0_bin_start = np.round(sig0_bin_start, decimals=2)
    sig0_bin_end = np.round(sig0_bin_end, decimals=2)
    for endmemname in endmemnames_to_use:
        #apply a filtering to endmemname_to_df to get the right
        # row corresponding to the sig0 range
        df = endmemname_to_df[endmemname]
        correct_rows_for_endmemname =\
        df[(df["sig0"] >= sig0_bin_start) & (df["sig0"] < sig0_bin_end)]
        #correct_rows_for_endmemname should have a length of 1 (there should
        # be only one row for each sig0 bin), so let's verify that with
        # an 'assert' statement.
        assert len(correct_rows_for_endmemname)==1
        #Store the row from this end member dataframe in
        # the list of correct rows
        correct_rows.append(correct_rows_for_endmemname)  
    #Now we just have to concatenate 'correct_rows' into a single
    # pandas DataFrame; pandas should take care of matching up the
    # columns correctly.
    paired_up_endmember_df = pd.concat(correct_rows) 

    return paired_up_endmember_df


class ThermoclineArrayOMPAProblem(object):

    def __init__(self, tc_lower_bound, tc_upper_bound, tc_step,
                       endmemnames_to_use, endmemname_to_df, obs_df,
                       paramsandweighting_conserved,
                       paramsandweighting_converted,
                       conversionratios,
                       watermassname_to_usagepenaltyfunc={}):
        self.tc_lower_bound = tc_lower_bound
        self.tc_upper_bound = tc_upper_bound
        self.tc_step = tc_step
        self.endmemnames_to_use = endmemnames_to_use
        self.endmemname_to_df = endmemname_to_df
        self.obs_df
        self.paramsandweighting_conserved = paramsandweighting_conserved
        self.paramsandweighting_converted = paramsandweighting_converted
        self.conversionratios = conversionratios
        self.watermassname_to_usagepenaltyfunc =\
            watermassname_to_usagepenaltyfunc

    def solve(self):      
        thermocline_ompa_results = []
        for sig0_bin_start in np.arange(self.tc_lower_bound,
                                        self.tc_upper_bound, self.tc_step):
            sig0_bin_end = sig0_bin_start + self.tc_step
            #Get the endmember dataframe for OMPA analysis corresponding to the
            #sig0 range 
            endmember_df_for_sig0_range =\
              get_endmember_df_for_sig0_range(
                  endmemnames_to_use=self.endmemnames_to_use,
                  endmemname_to_df=self.endmemname_to_df,
                  sig0_bin_start=sig0_bin_start,
                  sig0_bin_end=sig0_bin_end)

            #filter gp15_thermocline using sig0_bin_start and sig0_bin_end
            obs_df_for_sig0_range = self.obs_df[
                                  (self.obs_df["sig0"] >= sig0_bin_start)
                                  & (self.obs_df["sig0"] <= sig0_bin_end)]
            if (len(obs_df_for_sig0_range)==0):
              print("No observations for range", sig0_bin_start, sig0_bin_end)
              continue #skip this iteration of the loop
            
            #Now that you have the data frames for the observations and
            # end members, you can define the ompa problem
            ompa_problem = OMPAProblem(
                watermass_df=endmember_df_for_sig0_range,
                obs_df=obs_df_for_sig0_range, 
                paramsandweighting_conserved=self.paramsandweighting_conserved,
                paramsandweighting_converted=self.paramsandweighting_converted,
                conversionratios=self.conversionratios,
                smoothness_lambda=None,
                watermassname_to_usagepenaltyfunc=
                  self.watermassname_to_usagepenaltyfunc)
            ompa_problem.solve()
            if (ompa_problem.prob.status is not "infeasible"):
                thermocline_ompa_results.append(ompa_problem)
            else:
                print("Warning! Infeasible for:")
                print("obs df:")
                cols_to_print = (
                 [x[0] for x in thermocline_paramsandweighting[0]]
                 +[x[0] for x in thermocline_paramsandweighting[1]]
                 +["sig0"])
                print(obs_df_for_sig0_range[cols_to_print])
                print("water mass df:")
                print(endmember_df_for_sig0_range[cols_to_print])
                print("Try lowering the parameter weights!")

        self.thermocline_ompa_results = thermocline_ompa_results

        return thermocline_ompa_results
    

    
