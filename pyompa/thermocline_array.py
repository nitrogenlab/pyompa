from __future__ import division, print_function
from .ompacore import OMPAProblem, ExportToCsvMixin
import numpy as np
import pandas as pd
from collections import OrderedDict
from cvxpy.error import SolverError


def get_endmember_df_for_range(endmemnames_to_use,
                               endmemname_to_df,
                               endmember_name_column,
                               stratification_col,
                               bin_start, bin_end):
    #Idea: for each entry in endmemname_to_df, filter out the row
    # that lies within bin_start and bin_end (there should
    # only be one such row). Then assemble all those rows into a data frame
    # of the endmembers to use for the OMPA analysis, for that range

    #We'll create a list of the correct rows from each pandas df
    correct_rows = []
    bin_start = np.round(bin_start, decimals=2)
    bin_end = np.round(bin_end, decimals=2)
    endmem_names_present = []
    for endmemname in endmemnames_to_use:
        #apply a filtering to endmemname_to_df to get the right
        # row corresponding to the range
        df = endmemname_to_df[endmemname]
        correct_rows_for_endmemname = pd.DataFrame(
            df[(df[stratification_col] >= bin_start) &
               (df[stratification_col] < bin_end)])
        correct_rows_for_endmemname[endmember_name_column] = endmemname
        if len(correct_rows_for_endmemname)==0:
            continue
        else:
            endmem_names_present.append(endmemname) 
        #correct_rows_for_endmemname should have a length of 1 (there should
        # be only one row for each bin), so let's verify that with
        # an 'assert' statement.
        assert len(correct_rows_for_endmemname)==1, (
         "Too many rows for bin "+str(bin_start)+" to "+str(bin_end)
         +" for column "+stratification_col
         +":\n"+str(correct_rows_for_endmemname))
        #Store the row from this end member dataframe in
        # the list of correct rows
        correct_rows.append(correct_rows_for_endmemname)  
    #Now we just have to concatenate 'correct_rows' into a single
    # pandas DataFrame; pandas should take care of matching up the
    # columns correctly.
    paired_up_endmember_df = pd.concat(correct_rows) 

    return paired_up_endmember_df, endmem_names_present


class ThermoclineArraySoln(ExportToCsvMixin):

    def __init__(self, endmemname_to_df,
                       endmember_name_column,
                       endmemnames_to_use,
                       thermocline_ompa_problem, thermocline_ompa_results):
        self.endmemname_to_df = endmemname_to_df
        self.endmember_name_column = endmember_name_column
        self.endmemnames_to_use = endmemnames_to_use
        self.thermocline_ompa_problem = thermocline_ompa_problem
        self.thermocline_ompa_results = thermocline_ompa_results

        #get attributes in a format that is compatible with the API
        # of OMPASoln, for plotting and csv export purposes
        self.obs_df = pd.concat([x.obs_df for x in self])
        self.endmember_names = self[0].endmember_names
        self.param_names = self[0].param_names
        self.endmember_fractions = np.concatenate([
            x.endmember_fractions for x in self], axis=0) 
        self.converted_variables = np.concatenate([
            x.converted_variables for x in self], axis=0)
        self.param_residuals = np.concatenate([
            x.param_residuals for x in self], axis=0)
        #getting the groupname_to_xxx attributes ready
        self.groupname_to_totalconvertedvariable = OrderedDict()
        self.groupname_to_effectiveconversionratios = OrderedDict()
        for groupname in (self[0].groupname_to_totalconvertedvariable):
            self.groupname_to_totalconvertedvariable[groupname] =\
                np.concatenate([
                    x.groupname_to_totalconvertedvariable[groupname]
                    for x in self
                ], axis=0)
            self.groupname_to_effectiveconversionratios[groupname] =\
                    OrderedDict()
            for paramname in (self[0].
                            groupname_to_effectiveconversionratios[groupname]): 
                self.groupname_to_effectiveconversionratios[
                    groupname][paramname] = (
                 np.concatenate([
                  x.groupname_to_effectiveconversionratios[
                        groupname][paramname]
                  for x in self
                 ], axis=0))
        #getting endmembername to usage penalty
        self.endmembername_to_usagepenalty = OrderedDict()
        for endmembername in self.endmember_names:
            if endmembername in self[0].endmembername_to_usagepenalty: 
                self.endmembername_to_usagepenalty[endmembername] =\
                np.concatenate([x.endmembername_to_usagepenalty[endmembername]
                                for x in self], axis=0) 

    def core_quantify_ambiguity_via_residual_limits(self, *args, **kwargs):

        solns = [x.core_quantify_ambiguity_via_residual_limits(*args, **kwargs)
                 for x in self] 
        to_return =  ThermoclineArraySoln(
                    endmemname_to_df=self.endmemname_to_df,
                    endmember_name_column=self.endmember_name_column,
                    endmemnames_to_use=self.endmemnames_to_use,
                    thermocline_ompa_problem=None,
                    thermocline_ompa_results=solns) 
        to_return.perobs_obj = np.concatenate(
            [x.perobs_obj for x in to_return], axis=0) 
        return to_return

    def __len__(self):
        return len(self.thermocline_ompa_results)

    def __getitem__(self, i):
        return self.thermocline_ompa_results[i]

    def __iter__(self):
        return self.thermocline_ompa_results.__iter__()


class ThermoclineArrayOMPAProblem(object):

    def __init__(self, stratification_col,
                       tc_lower_bound, tc_upper_bound, tc_step,
                       obs_df,
                       **ompa_core_params):
        self.stratification_col = stratification_col
        self.tc_lower_bound = tc_lower_bound
        self.tc_upper_bound = tc_upper_bound
        self.tc_step = tc_step
        self.obs_df = obs_df
        self.ompa_core_params = ompa_core_params
        if (np.min(self.obs_df[self.stratification_col])
            < self.tc_lower_bound):
            print("==============================")
            print("Heads up! You specified a tc lower bound of",
                  tc_lower_bound,"but the observations df contains samples"
                  +" with a",self.stratification_col,
                  "as low as",np.min(self.obs_df[self.stratification_col]))
            print("==============================")
        if (np.max(self.obs_df[self.stratification_col])
            > self.tc_upper_bound):
            print("==============================")
            print("Heads up! You specified a tc upper bound of",
                  tc_upper_bound,"but the observations df contains samples"
                  +" with a",self.stratification_col,
                  "as high as",np.max(self.obs_df[self.stratification_col]))
            print("==============================")

    def solve(self, endmemname_to_df, endmember_name_column="endmember_name",
                    endmemnames_to_use=None,
                    **ompa_core_solve_params): 

        if (endmemnames_to_use is None):
            endmemnames_to_use = sorted(endmemname_to_df.keys())

        thermocline_ompa_results = []
        for bin_start in np.arange(self.tc_lower_bound,
                                   self.tc_upper_bound, self.tc_step):
            bin_end = bin_start + self.tc_step
            #Get the endmember dataframe for OMPA analysis corresponding to the
            #range 
            endmember_df_for_range, endmem_names_present =\
              get_endmember_df_for_range(
                  stratification_col=self.stratification_col,
                  endmemnames_to_use=endmemnames_to_use,
                  endmemname_to_df=endmemname_to_df,
                  endmember_name_column=endmember_name_column,
                  bin_start=bin_start,
                  bin_end=bin_end)

            #filter gp15_thermocline using bin_start and bin_end
            obs_df_for_range = self.obs_df[
                          (self.obs_df[self.stratification_col] >= bin_start)
                          & (self.obs_df[self.stratification_col] <= bin_end)]
            if (len(obs_df_for_range)==0):
              print("No observations for range", bin_start, bin_end)
              continue #skip this iteration of the loop
            
            #Now that you have the data frames for the observations and
            # end members, you can define the ompa problem
            #Do one observation at a time so that if the problem is infeasible
            # we don't lose more observations than needed
            for rowidx,one_obs in obs_df_for_range.iterrows():
                one_obs = one_obs.to_frame().transpose()
                print(one_obs)
                try:
                    ompa_soln = OMPAProblem(
                                 obs_df=one_obs,
                                 **self.ompa_core_params).solve(
                                   endmember_df=endmember_df_for_range,
                                   endmember_name_column=endmember_name_column,
                                   **ompa_core_solve_params)

                    if (ompa_soln.status != "infeasible"):
                        thermocline_ompa_results.append(
                            ompa_soln.insert_blank_endmembers_as_needed(
                                       new_endmember_names=endmemnames_to_use))
                    else:
                        print("Warning! Infeasible for:")
                        print("obs df:")
                        cols_to_print = self.ompa_core_params['param_names']
                        print(one_obs[cols_to_print])
                        print("endmember df:")
                        print(endmember_df_for_range[cols_to_print])
                        print("Try lowering the parameter weights!")
                except SolverError as e:
                    print("Encountered SolverError "+str(e))
                    print("obs df:")
                    cols_to_print = self.ompa_core_params['param_names']  
                    print(one_obs[cols_to_print])
                    print("endmember df:")
                    print(endmember_df_for_range[cols_to_print])

        self.thermocline_ompa_results = thermocline_ompa_results

        return ThermoclineArraySoln(
                 endmemname_to_df=endmemname_to_df,
                 endmember_name_column=endmember_name_column,
                 endmemnames_to_use=endmemnames_to_use,
                 thermocline_ompa_problem=self,
                 thermocline_ompa_results=thermocline_ompa_results)

