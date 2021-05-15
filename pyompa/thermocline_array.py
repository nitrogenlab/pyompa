from __future__ import division, print_function
from .ompacore import OMPAProblem
import numpy as np
import pandas as pd
from collections import OrderedDict


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
    for endmemname in endmemnames_to_use:
        #apply a filtering to endmemname_to_df to get the right
        # row corresponding to the range
        df = endmemname_to_df[endmemname]
        correct_rows_for_endmemname = pd.DataFrame(
            df[(df[stratification_col] >= bin_start) &
               (df[stratification_col] < bin_end)])
        correct_rows_for_endmemname[endmember_name_column] = endmemname
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

    return paired_up_endmember_df


class ThermoclineArraySoln(object):

    def __init__(self, endmemname_to_df,
                       endmember_name_column,
                       endmemnames_to_use,
                       thermocline_ompa_problem, thermocline_ompa_results):
        self.endmemname_to_df = endmemname_to_df
        self.endmember_name_column = endmember_name_column
        self.endmemnames_to_use = endmemnames_to_use
        self.thermocline_ompa_problem = thermocline_ompa_problem
        self.thermocline_ompa_results = thermocline_ompa_results

    def export_to_csv(self, csv_output_name,
                              orig_cols_to_include=[],
                              export_orig_param_vals=True,
                              export_residuals=True,
                              export_endmember_fracs=True,
                              export_oxygen_deficit=True,
                              export_conversion_ratios=True):
  
        print("writing to",csv_output_name)	
      
        toexport_df_dict = OrderedDict()

        conserved_params_to_use = (
            self.thermocline_ompa_results[0].conserved_params_to_use)
        converted_params_to_use = (
            self.thermocline_ompa_results[0].converted_params_to_use)
        param_names = conserved_params_to_use + converted_params_to_use 

        if (export_orig_param_vals):
            orig_cols_to_include += param_names

        for orig_col in orig_cols_to_include:
            assert orig_col in self.thermocline_ompa_results[0].obs_df, (
             "Export settings asked that "+orig_col+" be copied from the"
             +" original observations data frame into the output, but"
             +" no such column header was present in the observations data"
             +" frame; the column headers are: "
             +str(self.thermocline_ompa_results[0].columns)) 
            col_vals = np.concatenate([np.array(x.obs_df[orig_col])
                        for x in self.thermocline_ompa_results])
            toexport_df_dict[orig_col] = col_vals

        if (export_residuals):
             param_residuals = np.concatenate([
                x.param_residuals
                for x in self.thermocline_ompa_results], axis=0)
             for param_idx,param_name in enumerate(
                conserved_params_to_use+converted_params_to_use):
                toexport_df_dict[param_name+"_resid"] = (
                    param_residuals[:,param_idx])

        endmembernames=list(self.
         thermocline_ompa_results[0].endmember_df[
         self.thermocline_ompa_results[0].endmember_name_column])
        if (export_endmember_fracs):
            endmember_fractions = np.concatenate([
                x.endmember_fractions
                for x in self.thermocline_ompa_results], axis=0)
            for endmember_idx in range(len(endmembernames)):
                toexport_df_dict[endmembernames[endmember_idx]+"_frac"] =\
                    endmember_fractions[:,endmember_idx]

        if (export_oxygen_deficit and
            ((self.thermocline_ompa_results[0]
                  .total_oxygen_deficit) is not None)):
            total_oxygen_deficit = np.concatenate([x.total_oxygen_deficit
                               for x in self.thermocline_ompa_results], axis=0)
            toexport_df_dict["total_oxygen_deficit"] = total_oxygen_deficit

        if (export_conversion_ratios and
            ((self.thermocline_ompa_results[0]
                  .total_oxygen_deficit) is not None)): 
            for converted_param_idx in range(
                                        len(converted_params_to_use)):
                param_name = converted_params_to_use[converted_param_idx]
                effective_conversion_ratios = np.concatenate([
                    x.effective_conversion_ratios
                    for x in self.thermocline_ompa_results], axis=0)
                toexport_df_dict["oxygen_to_"+param_name+"_ratio"] =\
                    1.0/effective_conversion_ratios[:,converted_param_idx]

        toexport_df = pd.DataFrame(toexport_df_dict)
        toexport_df.to_csv(csv_output_name, index=False)

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
                       paramsandweighting_conserved,
                       paramsandweighting_converted,
                       conversionratios,
                       endmembername_to_usagepenaltyfunc={}):
        self.stratification_col = stratification_col
        self.tc_lower_bound = tc_lower_bound
        self.tc_upper_bound = tc_upper_bound
        self.tc_step = tc_step
        self.obs_df = obs_df
        self.paramsandweighting_conserved = paramsandweighting_conserved
        self.paramsandweighting_converted = paramsandweighting_converted
        self.conversionratios = conversionratios
        self.endmembername_to_usagepenaltyfunc =\
            endmembername_to_usagepenaltyfunc

    def solve(self, endmemname_to_df, endmember_name_column="endmember_name",
                    endmemnames_to_use=None): 
        if (endmemnames_to_use is None):
            endmemnames_to_use = sorted(endmemname_to_df.keys())
        thermocline_ompa_results = []
        for bin_start in np.arange(self.tc_lower_bound,
                                        self.tc_upper_bound, self.tc_step):
            bin_end = bin_start + self.tc_step
            #Get the endmember dataframe for OMPA analysis corresponding to the
            #range 
            endmember_df_for_range =\
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
            ompa_soln = OMPAProblem(
                obs_df=obs_df_for_range, 
                paramsandweighting_conserved=self.paramsandweighting_conserved,
                paramsandweighting_converted=self.paramsandweighting_converted,
                conversionratios=self.conversionratios,
                smoothness_lambda=None,
                endmembername_to_usagepenaltyfunc=
                  self.endmembername_to_usagepenaltyfunc).solve(
                    endmember_df=endmember_df_for_range,
                    endmember_name_column=endmember_name_column)
            if (ompa_soln.status != "infeasible"):
                thermocline_ompa_results.append(ompa_soln)
            else:
                print("Warning! Infeasible for:")
                print("obs df:")
                cols_to_print = (
                 [x[0] for x in thermocline_paramsandweighting[0]]
                 +[x[0] for x in thermocline_paramsandweighting[1]]
                 +[self.stratification_col])
                print(obs_df_for_range[cols_to_print])
                print("endmember df:")
                print(endmember_df_for_range[cols_to_print])
                print("Try lowering the parameter weights!")

        self.thermocline_ompa_results = thermocline_ompa_results

        return ThermoclineArraySoln(
                 endmemname_to_df=endmemname_to_df,
                 endmember_name_column=endmember_name_column,
                 endmemnames_to_use=endmemnames_to_use,
                 thermocline_ompa_problem=self,
                 thermocline_ompa_results=thermocline_ompa_results)
