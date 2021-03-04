from __future__ import division, print_function
import cvxpy as cp
import numpy as np
import pandas as pd
import scipy
import scipy.spatial
import scipy.linalg
from collections import OrderedDict


class OMPASoln(object):

    def __init__(self, endmember_df, endmember_name_column,
                       ompa_problem,
                       endmember_fractions,
                       oxygen_deficits,
                       total_oxygen_deficit,
                       param_residuals,
                       effective_conversion_ratios,
                       nullspace_A,
                       **kwargs):
        self.endmember_df = endmember_df
        self.endmember_name_column = endmember_name_column
        self.ompa_problem = ompa_problem
        self.endmember_fractions = endmember_fractions
        self.oxygen_deficits = oxygen_deficits
        self.total_oxygen_deficit = total_oxygen_deficit
        self.param_residuals = param_residuals
        self.effective_conversion_ratios = effective_conversion_ratios
        self.obs_df = ompa_problem.obs_df
        self.conserved_params_to_use = ompa_problem.conserved_params_to_use
        self.converted_params_to_use = ompa_problem.converted_params_to_use   
        self.endmembername_to_usagepenalty =\
            ompa_problem.endmembername_to_usagepenalty
        self.nullspace_A = nullspace_A
        self.__dict__.update(kwargs)

    def core_quantify_ambiguity_via_nullspace(self, obj_weights):
        #obj_weights should be an array of weights that define the objective
        # of the linear program, in the form "o @ (s + N(A) @ v)" (where
        # o is the objective, s is the solution, N(A) is the null space of A)

        #For each soln s (length "num_endmem + num_conv_ratios"), the
        # relevant linear programming constraints are
        #For rows pertaining to end member fractions:
        # "s + N(A) @ v >= 0" and "s + N(A) @ v <= 1" for endmember usgae
        # "s + N(A) @ v >= 0" OR "s + N(A) v <= 0" for oxygen usage 
        # (the 'or' means our solution space will basically contain two
        #  polytopes that may be disconnected (if there are multiple
        #  remin ratios; one polytope for positive and one for
        #  negative oxygen usage)
        #If there are any points with nonzero endmember usage constraints,
        # we can further add the constraint that P @ (N(A) @ v) = 0

        #Notes on how to determine the full feasible set (time complexity
        # seems to increase exponentially:)
        #If we are in |v| dimensional space, then the intersection of
        # |v| nonredundant equality constraints will define the corners. 
        #The available equality constraints are:
        # num_endmem constraints of the form endmem == 0
        # num_endmem constraints of the form endmem == 1
        # num_convratio constraints of the form endmem == 0
        # optionally one more constraint if there are nonzero usage penalties
        #Full solution:
        #From the (2*num_endmem + num_convratio + 1?) equations, we have
        # to check all (2*num_endmem + num_convratio + 1?)-choose-|v| subsets
        #Then for each subset, we need to solve the equation, and check
        # whether the solution satisfies the remaining constraints.
        #If the solution has negative O usage (has to be neg for all
        # conv ratios), we store the point in the negative O usage polytope;
        # otherwise, we store the point in the positive O usage polytope
        #We would then return the two polytopes
   
        #In practice, however, we are probably interested in quantifying the
        # limits of something (e.g. max and min values of an end member, or
        # max/min amounts of a tracer that could be explained by mixing),
        # in which case it is not necessary to solve for the full convex
        # polytope; the problem can just be solved as a linear program

        endmember_names = list(self.endmember_df[self.endmember_name_column])
        endmember_usagepenalty =\
            self.ompa_problem.prep_endmember_usagepenalty_mat(endmember_names)

        #For each observation, we can solve the linear program
        new_perobs_endmember_fractions = []
        new_perobs_oxygen_deficits = []
        perobs_obj = []
        for obs_idx in range(len(self.endmember_fractions)):
            endmem_fracs = self.endmember_fractions[obs_idx] 
            assert len(endmem_fracs)==len(endmember_names)
            oxy_def = self.oxygen_deficits[obs_idx] 
            usagepenalty = endmember_usagepenalty[obs_idx]

            assert len(usagepenalty) == len(endmem_fracs)
            assert ((len(endmem_fracs) + len(oxy_def))
                    == self.nullspace_A.shape[0])
            assert len(obj_weights) == self.nullspace_A.shape[0]

            def compute_soln(oxy_def_sign):
                v = cp.Variable(shape=(self.nullspace_A.shape[1]))
                obj = cp.Maximize(cp.sum(obj_weights @ (self.nullspace_A@v)))
                endmem_frac_deltas = self.nullspace_A[:len(endmem_fracs)] @ v
                new_endmem_fracs = (endmem_fracs + endmem_frac_deltas)
                new_oxy_def = (oxy_def
                  + self.nullspace_A[len(endmem_fracs):] @ v)
                constraints = [
                   new_endmem_fracs >= 0,
                   new_endmem_fracs <= 1,
                   new_oxy_def*oxy_def_sign >= 0,
                   usagepenalty@endmem_frac_deltas == 0 
                ] 
                prob = cp.Problem(obj, constraints)
                prob.solve(verbose=False, max_iter=50000)
                if (prob.status == "infeasible"):
                    return (None, None), -np.inf
                else:
                    return ((new_endmem_fracs.value, new_oxy_def.value),
                            prob.value) #soln and optimal value

            if (self.nullspace_A.shape[1] > 0):
                posodsign_soln, posodsign_obj = compute_soln(1) 
                negodsign_soln, negodsign_obj = compute_soln(-1) 
                pos_wins = (posodsign_obj > negodsign_obj)
                if (pos_wins):
                    (new_endmem_fracs, new_oxy_def) = (
                     posodsign_soln, posodsign_obj)
                else:
                    (new_endmem_fracs, new_oxy_def) = (
                     negodsign_soln, negodsign_obj)

                assert new_endmem_fracs is not None
                assert np.abs(np.sum(new_endmem_fracs) - 1) < 1e-7
                #fix numerical issues with soln
                new_endmem_fracs = np.maximum(new_endmem_fracs, 0)
                new_endmem_fracs = new_endmem_fracs/(np.sum(new_endmem_fracs))
                if (pos_wins):
                    new_oxy_def = np.maximum(new_oxy_def, 0)
                else:
                    new_oxy_def = np.minimum(new_oxy_def, 0)
            else:
                new_endmem_fracs = endmem_fracs
                new_oxy_def = oxy_def 

            obj = obj_weights@np.concatenate(
                   [new_endmem_fracs, new_oxy_def], axis=0) 

            new_perobs_endmember_fractions.append(new_endmem_fracs)
            new_perobs_oxy_defs.append(new_oxy_def)
            perobs_obj.append(obj) 

        return (np.array(new_perobs_endmember_fractions),
                np.array(new_perobs_oxy_defs),
                np.array(perobs_obj))

    def export_to_csv(self, csv_output_name,
                            orig_cols_to_include=[],
                            export_orig_param_vals=True,
                            export_residuals=True,
                            export_endmember_fracs=True,
                            export_oxygen_deficit=True,
                            export_conversion_ratios=True,
                            export_endmember_usage_penalties=False):

        print("writing to",csv_output_name)

        toexport_df_dict = OrderedDict()

        if (export_orig_param_vals):
            orig_cols_to_include += (
             self.conserved_params_to_use+self.converted_params_to_use)

        for orig_col in orig_cols_to_include:
            assert orig_col in self.obs_df, (
             "Export settings asked that "+orig_col+" be copied from the"
             +" original observations data frame into the output, but"
             +" no such column header was present in the observations data"
             +" frame; the column headers are: "+str(self.obs_df.columns)) 
            toexport_df_dict[orig_col] = self.obs_df[orig_col]

        if (export_residuals):
             for param_idx,param_name in enumerate(
                self.conserved_params_to_use+self.converted_params_to_use):
                toexport_df_dict[param_name+"_resid"] =\
                    self.param_residuals[:,param_idx] 

        endmembernames=list(self.endmember_df[self.endmember_name_column])
        if (export_endmember_fracs):
            for endmember_idx in range(len(endmembernames)):
                toexport_df_dict[endmembernames[endmember_idx]+"_frac"] =\
                    self.endmember_fractions[:,endmember_idx]

        if (export_oxygen_deficit and
            (self.total_oxygen_deficit is not None)):
            toexport_df_dict["total_oxygen_deficit"] = self.total_oxygen_deficit

        if (export_conversion_ratios and
            (self.total_oxygen_deficit is not None)): 
            for converted_param_idx in range(len(self.converted_params_to_use)):
                param_name = self.converted_params_to_use[converted_param_idx]
                toexport_df_dict["oxygen_to_"+param_name+"_ratio"] =\
                    1.0/self.effective_conversion_ratios[:,converted_param_idx]

        if (export_endmember_usage_penalties):
            for endmembername in endmembernames:
                if (endmembername in\
                    self.ompa_problem.endmembername_to_usagepenalty): 
                    endmember_usagepenalty = (self.ompa_problem.
                                  endmembername_to_usagepenalty[endmembername])
                    toexport_df_dict[endmembername+"_penalty"] =\
                        endmember_usagepenalty
        
        toexport_df = pd.DataFrame(toexport_df_dict)
        toexport_df.to_csv(csv_output_name, index=False)

    def iteratively_refine_ompa_soln(self, num_iterations):
        init_endmember_df = self.ompa_problem.construct_ideal_endmembers(
            ompa_soln=self)
        ompa_solns = self.ompa_problem.iteratively_refine_ompa_solns(
            init_endmember_df=init_endmember_df,
            endmember_name_column=self.endmember_name_column,
            num_iterations=num_iterations)
        return [self]+ompa_solns


class OMPAProblem(object):
    """
        Core class for conducting OMPA analysis using cvxpy
    """     
    def __init__(self, obs_df,
                       paramsandweighting_conserved,
                       paramsandweighting_converted,
                       conversionratios,
                       smoothness_lambda,
                       endmembername_to_usagepenaltyfunc):
        self.obs_df = obs_df
        self.paramsandweighting_conserved = paramsandweighting_conserved
        self.paramsandweighting_converted = paramsandweighting_converted
        self.conversionratios = conversionratios
        self.smoothness_lambda = smoothness_lambda
        self.endmembername_to_usagepenaltyfunc =\
          endmembername_to_usagepenaltyfunc
        self.process_params()
        self.prep_endmember_usagepenalties() 


    def process_params(self):
        #check that every param in self.paramsandweighting_converted is
        # specified in convertedparams_ratios
        
        #paramsandweighting_conserved is a list of tuples; split them up
        self.conserved_params_to_use, self.conserved_weighting = [
          list(x) for x in zip(*self.paramsandweighting_conserved)]
        if (len(self.paramsandweighting_converted) > 0):
            self.converted_params_to_use, self.converted_weighting = [
              list(x) for x in zip(*self.paramsandweighting_converted)]
        else:
            self.converted_params_to_use, self.converted_weighting = [], []

        param_names = self.converted_params_to_use+self.conserved_params_to_use 
        for param_name in (param_names):
            assert param_name in self.obs_df,\
                (param_name+" not specified in observations data frame;"
                 +" observations data frame columns are "
                 +str(list(self.obs_df.columns)))

        with_drop_na = self.obs_df.dropna(subset=param_names)
        if (len(with_drop_na) < len(self.obs_df)):
            print("Dropping "+str(len(self.obs_df)-len(with_drop_na))
                  +" rows that have NA values in the observations")
            self.obs_df = with_drop_na

        if (max(self.converted_weighting+self.conserved_weighting) > 100):
            print("Warning: having very large param weights can lead to"
                  +" instability in the optimizer! Consider scaling down"
                  +" your weights")
        
        #make sure every parameter in converted_params_to_use is specified in
        # convertedparams_ratios:
        assert all([(x in self.conversionratios)
                     for x in self.converted_params_to_use])
        if (len(self.converted_params_to_use) > 0):
            #also assert that every entry in convertedratios has the same length
            assert len(set([len(x) for x in 
                       self.conversionratios.values()])) == 1, conversionratios
            self.num_conversion_ratios = len(
                list(self.conversionratios.values())[0])
        else:
            self.num_conversion_ratios = 0

    def prep_endmember_usagepenalties(self):
        self.endmembername_to_usagepenalty = OrderedDict()
        for endmembername, penalty_func in\
         self.endmembername_to_usagepenaltyfunc.items():
            print("Adding penalty for",endmembername)
            penalty = penalty_func(self.obs_df)
            self.endmembername_to_usagepenalty[endmembername] = penalty

    def prep_endmember_usagepenalty_mat(self, endmember_names):
        endmember_usagepenalty = np.zeros((len(self.obs_df),
                                           len(endmember_names)))
        for endmemberidx,endmembername in enumerate(endmember_names):
            if endmembername in self.endmembername_to_usagepenalty:
                endmember_usagepenalty[:,endmemberidx] =\
                    self.endmembername_to_usagepenalty[endmembername]
        #print a warning if specified a usage penalty that was not used
        for endmembername in self.endmembername_to_usagepenalty:
            if endmembername not in endmember_names:
                raise RuntimeError("You specified a usage penalty for "
                 +endmembername+" but that endmember did not appear "
                 +"in the endmember data frame; endmembers are "
                 +str(endmember_names))
                 
        return endmember_usagepenalty

    def get_b(self):
        b = np.array(self.obs_df[self.conserved_params_to_use
                                 +self.converted_params_to_use])
        return b

    def get_param_weighting(self):
        return np.array(self.conserved_weighting+self.converted_weighting)

    def get_conversion_ratio_rows_of_A(self):
        #conversion_ratios has dimensions:
        # num_conversion_ratios x num_converted_params
        conversion_ratios = np.array([
            [self.conversionratios[param][i] for param
             in self.converted_params_to_use]
           for i in range(self.num_conversion_ratios)])
        return conversion_ratios, np.array([
            [0 for x in self.conserved_params_to_use]+list(converted_ratio_row)
            for converted_ratio_row in conversion_ratios])

    def get_endmem_mat(self, endmember_df):
        return np.array(endmember_df[self.conserved_params_to_use
                                     +self.converted_params_to_use])

    def get_nullspace(self, M, R):
        #Let M represent the end-member matrix - dimensions of end_mem x params
        # the R represent the remin ratios - dims of num_remin_ratios x params
        # |M.T R.T| x = 0
        # |1   0  | 
        # (where |1 0| represents the mass conservation constraint)
        #The solutions of x will give us the null space
        mat = np.concatenate(
               [(np.concatenate([np.transpose(M), np.transpose(R)], axis=1)
                 if len(R) > 0 else np.transpose(M)),
                np.array([1 for i in range(len(M))]
                         +[0 for i in range(len(R))])[None,:]
               ], axis=0)
        ns = scipy.linalg.null_space(mat) 
        if ns.shape[1] > 0:
            assert np.allclose(mat.dot(ns), 0)
        return ns

    def solve(self, endmember_df, endmember_name_column):

        for param_name in (self.conserved_params_to_use
                           +self.converted_params_to_use):
            assert param_name in endmember_df,\
                (param_name+" not specified in endmemebr_df where columns are "
                 +str(list(endmember_df.columns)))

        endmember_names = list(endmember_df[endmember_name_column])

        weighting = self.get_param_weighting() 
        smoothness_lambda = self.smoothness_lambda

        endmember_usagepenalty =\
            self.prep_endmember_usagepenalty_mat(endmember_names)
        self.endmember_usagepenalty = endmember_usagepenalty

        #Prepare A
        conversion_ratios, conversion_ratio_rows =\
            self.get_conversion_ratio_rows_of_A()
        #conversion_ratio_rows = 'R'
        endmem_mat = self.get_endmem_mat(endmember_df) #'M'
        #add a row to A for the ratios
        if (len(conversion_ratio_rows > 0)):
            A = np.concatenate([endmem_mat, conversion_ratio_rows], axis=0)
        else:
            A = endmem_mat

        #compute the nullspace of A - will be useful for disentangling
        # ambiguity in the solution
        nullspace_A = self.get_nullspace(M=endmem_mat,
                                         R=conversion_ratio_rows)

        #prepare b
        b = self.get_b()
        
        #Rescale by param weighting
        print("params to use:", self.conserved_params_to_use,
                                self.converted_params_to_use)        
        print("param weighting:", weighting)
        print("Conversion ratios:\n"+str(conversion_ratios))
        A = A*weighting[None,:]
        b = b*weighting[None,:]

        if (smoothness_lambda is not None):
            pairs_matrix = make_pairs_matrix(
              obs_df=self.obs_df,
              depth_metric="depth",
              depth_scale=1.0,
              nneighb=4)
        else:
            pairs_matrix = None

        #first run with only a positive conversion ratio allowed
        _, _, _, perobs_weighted_resid_sq_positiveconversionsign, _ =\
          self.core_solve(
            A=A, b=b,
            num_conversion_ratios=self.num_conversion_ratios,
            num_converted_params=len(self.converted_params_to_use),
            pairs_matrix=None,
            endmember_usagepenalty=endmember_usagepenalty,
            conversion_sign_constraints=1,
            smoothness_lambda=None)
        _, _, _, perobs_weighted_resid_sq_negativeconversionsign, _ =\
          self.core_solve(
            A=A, b=b,
            num_conversion_ratios=self.num_conversion_ratios,
            num_converted_params=len(self.converted_params_to_use),
            pairs_matrix=None,
            endmember_usagepenalty=endmember_usagepenalty,
            conversion_sign_constraints=-1,
            smoothness_lambda=None)
        
        #determine which conversion sign is better
        positive_conversionsign_isbetter = (
            perobs_weighted_resid_sq_positiveconversionsign <
            perobs_weighted_resid_sq_negativeconversionsign)
        final_conversion_signconstraints = (
            1.0*positive_conversionsign_isbetter
            + -1.0*(positive_conversionsign_isbetter==False))
        
        (x, endmember_fractions,
         oxygen_deficits,
         perobs_weighted_resid_sq, prob) = self.core_solve(
            A=A, b=b,
            num_conversion_ratios=self.num_conversion_ratios,
            num_converted_params=len(self.converted_params_to_use),
            pairs_matrix=pairs_matrix,
            endmember_usagepenalty=endmember_usagepenalty,
            conversion_sign_constraints=final_conversion_signconstraints,
            smoothness_lambda=smoothness_lambda)
        
        if (endmember_fractions is not None):
            print("objective:", np.sum(perobs_weighted_resid_sq))
            param_reconstruction = (x@A)/weighting[None,:]
            param_residuals = b/weighting[None,:] - param_reconstruction
        else:
            param_residuals = None

        if (oxygen_deficits is not None):
            #sanity check the signs of the oxygen deficits; for each entry they
            # should either be all positive or all negative, within numerical
            # precision
            for oxygen_deficit in oxygen_deficits:
                if (len(oxygen_deficit) > 0):
                    if ((all(oxygen_deficit >= 0) or
                         all(oxygen_deficit <= 0))==False):
                        print("WARNING: sign inconsistency in"
                              +" oxygen deficits:", oxygen_deficit)
            total_oxygen_deficit = np.sum(oxygen_deficits, axis=-1)
            #proportions of oxygen use at differnet ratios
            with np.errstate(divide='ignore', invalid='ignore'):
                oxygen_usage_proportions = (
                    oxygen_deficits/total_oxygen_deficit[:,None])
                #Reminder: conversion_ratios has dims of
                # num_conversion_ratios x num_converted_params
                #oxygen_usage_proportions has dims of
                # num_examples X num_conversion_ratios
                #Note: should disregard conversion ratios computed when oxygen
                # usage is very small.
                effective_conversion_ratios = (
                    oxygen_usage_proportions@conversion_ratios)         
        else:
            total_oxygen_deficit = None
            effective_conversion_ratios = None
           

        return OMPASoln(endmember_df=endmember_df,
                  ompa_problem=self,
                  endmember_names=endmember_names,
                  endmember_name_column=endmember_name_column,
                  status=prob.status,
                  endmember_fractions=endmember_fractions,
                  oxygen_deficits=oxygen_deficits,
                  resid_wsumsq=np.sum(perobs_weighted_resid_sq),
                  param_residuals=param_residuals,
                  total_oxygen_deficit=total_oxygen_deficit,
                  effective_conversion_ratios=effective_conversion_ratios,
                  nullspace_A=nullspace_A)

    def core_solve(self, A, b, num_conversion_ratios, num_converted_params,
                   pairs_matrix, endmember_usagepenalty,
                   conversion_sign_constraints, smoothness_lambda,
                   verbose=True):
  
        #We are going to solve the following problem:
        #P is the penalty matrix. It has dimensions of
        #  (observations X end_members)
        #Minimize (x@A - b)^2 + (x[:,:-num_conversion_ratios]*P)^2
        #Subject to x[:,:-num_conversion_ratios] >= 0,
        #           cp.sum(x[:,:-num_conversion_ratios], axis=1) == 1
        # x has dimensions of observations X (end_members+num_conversion_ratios)
        # the +1 row represents O2 deficit, for remineralization purposes
        # A has dimensions of (end_members+num_conversion_ratios) X parameteres
        # b has dimensions of observations X parameters 
        
        num_endmembers = len(A)-num_conversion_ratios
        x = cp.Variable(shape=(len(b), len(A)))
        obj = (cp.sum_squares(x@A - b) +
                cp.sum_squares(cp.atoms.affine.binary_operators.multiply(
                                    x[:,:num_endmembers],
                                    endmember_usagepenalty) ))
        if (smoothness_lambda is not None):
            #leave out O2 deficit column from the smoothness penality as it's
            # on a bit of a different scale.
            obj += smoothness_lambda*cp.sum_squares(
                    pairs_matrix@x[:,:num_endmembers])
        obj = cp.Minimize(obj)
        
        #leave out the last column as it's the conversion ratio
        constraints = [
           x[:,:num_endmembers] >= 0,
           cp.sum(x[:,:num_endmembers],axis=1)==1]
        if (num_conversion_ratios > 0):
            if (hasattr(conversion_sign_constraints, '__len__')==False):
                constraints.append(
                  cp.atoms.affine.binary_operators.multiply(
                      conversion_sign_constraints,
                      x[:,num_endmembers:]) >= 0)
            else:
                constraints.append(
                  cp.atoms.affine.binary_operators.multiply(
                      np.tile(A=conversion_sign_constraints[:,None],
                              reps=(1,num_conversion_ratios)),
                      x[:,num_endmembers:]) >= 0)
        prob = cp.Problem(obj, constraints)
        prob.solve(verbose=False, max_iter=50000)
        #settign verbose=True will generate more print statements and slow down the analysis
        
        print("status:", prob.status)
        print("optimal value", prob.value)

        if (prob.status=="infeasible"):
            raise RuntimeError("Optimization failed - "
                               +"try lowering the parameter weights?")
        else:
            #weighted sum of squared of the residuals
            original_resid_wsumsq = np.sum(np.square((x.value@A) - b))

            endmember_fractions = x.value[:,:num_endmembers]
            ##enforce the constraints (nonnegativity, sum to 1) on
            ## endmember_fractions
            endmember_fractions = np.maximum(endmember_fractions, 0) 
            endmember_fractions = (endmember_fractions/
                np.sum(endmember_fractions,axis=-1)[:,None])

            if (num_converted_params > 0):
                oxygen_deficits = x.value[:,num_endmembers:]
                if (hasattr(conversion_sign_constraints, '__len__')):
                    gt_zero_mask = conversion_sign_constraints > 0
                    oxygen_deficits[gt_zero_mask] = np.maximum(
                        oxygen_deficits[gt_zero_mask], 0.0) 
                    lt_zero_mask = conversion_sign_constraints < 0
                    oxygen_deficits[lt_zero_mask] = np.minimum(
                        oxygen_deficits[lt_zero_mask], 0.0) 
                else:
                    if (conversion_sign_constraints > 0):
                        oxygen_deficits = np.maximum(oxygen_deficits, 0.0)
                    elif (conversion_sign_constraints < 0):
                        oxygen_deficits = np.minimum(oxygen_deficits, 0.0)
                    else:
                        assert False, conversion_sign_constraints
                #fixed_x is x that is forced to satisfy the constraints
                fixed_x = np.concatenate(
                            [endmember_fractions, oxygen_deficits], axis=-1)
            else:
                oxygen_deficits = None
                fixed_x = np.concatenate([endmember_fractions,
                 np.zeros((len(endmember_fractions), num_conversion_ratios))],
                 axis=-1)

            afterfixing_resid_wsumsq = np.sum(np.square(fixed_x@A - b))
            print("Original weighted sum squares:",original_resid_wsumsq)
            print("Post fix weighted sum squared:",afterfixing_resid_wsumsq)

            perobs_weighted_resid_sq =\
                np.sum(np.square((fixed_x@A) - b), axis=-1)
            #TODO: enforce that the constraints are satisfied, and
            #recompute residuals accordingly.
        
        return (fixed_x, endmember_fractions, oxygen_deficits,
                perobs_weighted_resid_sq, prob)

    def construct_ideal_endmembers(self, ompa_soln):

        print("Constructing ideal end members")

        b = self.get_b() # dims of num_obs X params

        if (ompa_soln.oxygen_deficits is not None):
            #self.oxygen_deficits has dims of num_obs X num_conversion_ratios
            #conversion_ratio_rows has dims of num_conversion_ratios X params
            conversion_ratios, conversion_ratio_rows =\
                self.get_conversion_ratio_rows_of_A()
            deltas_due_to_oxygen_deficits = ( #<- num_obs X params
             ompa_soln.oxygen_deficits@conversion_ratio_rows) 
            b = b - deltas_due_to_oxygen_deficits

        #Do a sanity check to make sure that, with the existing end member
        # matrix, we end up recapitulating the residuals
        #existing_endmemmat has dims of num_endmembers X params
        existing_endmemmat = self.get_endmem_mat(ompa_soln.endmember_df)
        #self.endmember_fractions has dims of num_obs X num_endmembers
        #Note: b and existing_endmemmat haven't been scaled by weighting yet,
        # so there is no need
        # to un-scale it here. Also note deltas_due_to_oxygen_deficits has
        # already been subtracted
        old_param_residuals = (b
          - ompa_soln.endmember_fractions@existing_endmemmat)
        np.testing.assert_almost_equal(old_param_residuals,
                                       ompa_soln.param_residuals, decimal=5)

        #Now solve for a better end-member mat
        weighting = self.get_param_weighting()
        new_endmemmat_var =  cp.Variable(shape=existing_endmemmat.shape)  
        
        #keeping the endmember_fractions, what are the best end members?
        obj = cp.Minimize(
            cp.sum_squares(ompa_soln.endmember_fractions@new_endmemmat_var
                           - b*weighting[None,:]))
        constraints = [] #no constraints for now
        prob = cp.Problem(obj, constraints)
        prob.solve(verbose=False, max_iter=50000)
        
        print("status:", prob.status)
        print("optimal value", prob.value)

        if (prob.status=="infeasible"):
            raise RuntimeError("Something went wrong "
                               +"- the optimization failed")
        else:
            new_endmemmat = new_endmemmat_var.value/weighting[None,:]
            #Sanity check that the residuals got better
            new_param_residuals = (b
                - ompa_soln.endmember_fractions@new_endmemmat)
            new_param_resid_wsumsq =\
                np.sum(np.square(new_param_residuals*weighting[None,:]))
            old_param_resid_wsumsq =\
                np.sum(np.square(old_param_residuals*weighting[None,:]))
            print("Old weighted residuals sumsquared:",old_param_resid_wsumsq)
            print("New weighted residuals sumsquared:",new_param_resid_wsumsq)
            assert new_param_resid_wsumsq <= old_param_resid_wsumsq

            #make a endmember data frame
            new_endmemmat_df = pd.DataFrame(OrderedDict(
             [(ompa_soln.endmember_name_column,
               list(ompa_soln.endmember_df[ompa_soln.endmember_name_column]))]
             +[(paramname, values) for paramname,values in
             zip(self.conserved_params_to_use+self.converted_params_to_use,
                 new_endmemmat.T) ])) 
            return new_endmemmat_df

    def iteratively_refine_ompa_solns(self, init_endmember_df,
            endmember_name_column, num_iterations):
        assert num_iterations >= 1,\
            "num_iterations must be >= 1; is "+str(num_iterations)
        print("On iteration 1")
        ompa_solns = [self.solve(endmember_df=init_endmember_df,
                                 endmember_name_column=endmember_name_column)]
        for i in range(1,num_iterations):
            print("On iteration "+str(i+1))
            new_endmember_df =\
                self.construct_ideal_endmembers(ompa_soln=ompa_solns[-1]) 
            ompa_solns.append(self.solve(endmember_df=new_endmember_df,
                                  endmember_name_column=endmember_name_column))
        return ompa_solns 


def spherical_to_surface_cartesian(lat, lon):
    r = 6.371*(1E3) #earth radius
    theta = ((1-lat)/180.0)*np.pi
    phi = (lon/180.0)*np.pi
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    return (x,y)


def add_surface_cartesian_coordinates_to_df(df):
    latitudes = df["latitude"]
    longitudes = df["longitude"]
    xs,ys = list(zip(*[spherical_to_surface_cartesian(*x)
                       for x in zip(latitudes, longitudes)]))
    df["x"] = xs
    df["y"] = ys
    #plt.scatter(xs, ys)
    #plt.show()


def compute_pairwise_distances_depthmetric(df, depth_metric, depth_scale):
    xs = df["x"]
    ys = df["y"]
    
    depth_metric = np.array(df[depth_metric])
    depth_diffs = np.abs(depth_metric[:,None] -
                         depth_metric[None,:])*depth_scale

    #plt.hist(depth_diffs.ravel(), bins=20)
    #plt.show()

    coors = np.array([xs, ys]).transpose((1,0))
    euclidean_distances = scipy.spatial.distance.squareform(
        scipy.spatial.distance.pdist(coors))
    #plt.hist(euclidean_distances.ravel(), bins=100)
    #plt.show()

    weighted_distances = np.sqrt(np.square(euclidean_distances)
                                 + np.square(depth_diffs))
    #plt.hist(weighted_distances.ravel(), bins=20)
    #plt.show()
    return weighted_distances


def make_pairs_matrix(obs_df, depth_metric, depth_scale, nneighb):
    obs_df = pd.DataFrame(obs_df)
    add_surface_cartesian_coordinates_to_df(obs_df)
    pairwise_distances = compute_pairwise_distances_depthmetric(
        obs_df, depth_metric=depth_metric, depth_scale=depth_scale)
    #plt.hist(pairwise_distances.ravel(), bins=20)
    #plt.show()
    nneighb_thresh = np.sort(pairwise_distances, axis=-1)[:,nneighb]
    masked_pairwise_distances =\
      (pairwise_distances*(pairwise_distances <= nneighb_thresh[:,None])
                         *(pairwise_distances > 0))
    pairs_to_consider_indices = np.nonzero(masked_pairwise_distances)
    print("Constrained pairs:",len(pairs_to_consider_indices[0]))
    pairs_distances = pairwise_distances[
        pairs_to_consider_indices[0],
        pairs_to_consider_indices[1]]
    #plt.hist(pairs_distances.ravel(), bins=20)
    #plt.show()
    pairs_matrix = np.zeros((len(pairs_to_consider_indices[0]),
                              len(obs_df)))
    pairs_matrix[np.arange(len(pairs_distances)),
                  pairs_to_consider_indices[0]] = 1.0/nneighb#(
                      #1/pairs_distances)
    pairs_matrix[np.arange(len(pairs_distances)),
                  pairs_to_consider_indices[1]] = -1.0/nneighb#(
                      #1/pairs_distances)
    return pairs_matrix
    
