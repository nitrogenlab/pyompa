from __future__ import division, print_function
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from .util import collapse_endmembers_by_idxmapping
from collections import OrderedDict


def plot_endmember_usagepenalties(endmembername_to_usagepenalty,
                                  xaxis_vals, xaxis_label,
                                  yaxis_vals, yaxis_label,
                                  flip_y=True):
    endmembernames = sorted(endmembername_to_usagepenalty)
    num_figs = len(endmembernames)
    fig, ax = plt.subplots(nrows=1, ncols=num_figs, figsize=(5*num_figs,4))
    for i in range(len(endmembernames)):
        endmembername = endmembernames[i]
        plt.sca(ax[i])
        plt.scatter(xaxis_vals, yaxis_vals,
                    c=endmembername_to_usagepenalty[endmembername])
        plt.xlabel(xaxis_label)
        if (i==0):
            plt.ylabel(yaxis_label)
        if (flip_y):
            plt.ylim(plt.ylim()[1], plt.ylim()[0])
        plt.colorbar()
        plt.title(endmembername)
    plt.show()


def plot_ompasoln_endmember_usagepenalties(
        ompa_soln, xaxis_colname, yaxis_colname, flip_y=True):
    plot_endmember_usagepenalties(
        endmembername_to_usagepenalty=
            ompa_soln.endmembername_to_usagepenalty,
        xaxis_vals=ompa_soln.obs_df[xaxis_colname],
        xaxis_label=xaxis_colname,
        yaxis_vals=ompa_soln.obs_df[yaxis_colname],
        yaxis_label=yaxis_colname,
        flip_y=True)


def plot_endmember_fractions(xaxis_vals, xaxis_label, yaxis_vals, yaxis_label,
        endmember_fractions, endmembernames,
        groupname_to_totalconvertedvariable,
        groupname_to_effectiveconversionratios,
        flip_y=True):
    num_endmembers = endmember_fractions.shape[1]
    num_figs = (num_endmembers +
                len(groupname_to_totalconvertedvariable) +
                sum([len(x) for x in
                     groupname_to_effectiveconversionratios.values()]))
    fig, ax = plt.subplots(nrows=1, ncols=num_figs,
                           figsize=(5*num_figs,4))
    for i in range(num_endmembers):
        plt.sca(ax[i])
        plt.scatter(xaxis_vals, yaxis_vals, c=endmember_fractions[:,i])
        plt.xlabel(xaxis_label)
        if (i==0):
            plt.ylabel(yaxis_label)
        if (flip_y):
            plt.ylim(plt.ylim()[1], plt.ylim()[0])
        plt.colorbar()
        plt.clim(0.0, max(1.0,np.max(endmember_fractions[:,i])) )
        plt.title(endmembernames[i])
    plotidx = num_endmembers
    for groupname in groupname_to_totalconvertedvariable:
        plt.sca(ax[plotidx])
        plotidx += 1
        convar_vals = groupname_to_totalconvertedvariable[groupname]
        plt.scatter(xaxis_vals, yaxis_vals,
                    c=convar_vals,
                    cmap=("RdBu"))
        max_abs_convval = np.max(np.abs(convar_vals))
        if (flip_y):
            plt.ylim(plt.ylim()[1], plt.ylim()[0])
        plt.colorbar()
        plt.clim(-max_abs_convval, max_abs_convval)
        plt.xlabel(xaxis_label)
        plt.title(groupname)
        effective_conversion_ratios =\
            groupname_to_effectiveconversionratios[groupname]
        for converted_param in effective_conversion_ratios:
            plt.sca(ax[plotidx])
            plotidx += 1
            plt.scatter(xaxis_vals, yaxis_vals,
                        c=effective_conversion_ratios[converted_param])
            if (flip_y):
                plt.ylim(plt.ylim()[1], plt.ylim()[0])
            plt.colorbar()
            plt.xlabel(xaxis_label)
            plt.title(converted_param+":"+groupname+" ratio")
    plt.show()


def plot_ompasoln_endmember_fractions(ompa_soln, xaxis_colname,
                                      yaxis_colname, flip_y=True,
                                      group_endmembers=True):

    if (group_endmembers):
        endmembername_to_indices = ompa_soln.endmembername_to_indices
        endmember_names = list(endmembername_to_indices.keys())
        endmember_fractions = collapse_endmembers_by_idxmapping(
            endmember_fractions=ompa_soln.endmember_fractions,
            endmembername_to_indices=endmembername_to_indices) 
    else:
        endmember_names = ompa_soln.endmember_names 
        endmember_fractions = ompa_soln.endmember_fractions

    plot_endmember_fractions(
        xaxis_vals=ompa_soln.obs_df[xaxis_colname],
        xaxis_label=xaxis_colname,
        yaxis_vals=ompa_soln.obs_df[yaxis_colname],
        yaxis_label=yaxis_colname,
        endmember_fractions=endmember_fractions,
        endmembernames=endmember_names,
        groupname_to_totalconvertedvariable=
         ompa_soln.groupname_to_totalconvertedvariable,
        groupname_to_effectiveconversionratios=
         ompa_soln.groupname_to_effectiveconversionratios,
        flip_y=flip_y)


def plot_residuals(param_residuals, param_names, xaxis_vals, xaxis_label,
                   yaxis_vals, yaxis_label, flip_y=True,
                   param_residual_weights=None):
    num_params = param_residuals.shape[1]
    ncols = num_params + (1 if param_residual_weights is not None
                          else 0)
    fig, ax = plt.subplots(nrows=1,
                           ncols=num_cols,
                           figsize=(5*num_cols,4))
    for i in range(param_residuals.shape[1]):
        plt.sca(ax[i])
        param_resid_maxabs = np.max(np.abs(param_residuals[:,i]))
        plt.scatter(x=xaxis_vals,
                    y=yaxis_vals,
                    c=param_residuals[:,i],
                    cmap="RdBu")
        plt.colorbar()
        plt.clim(-param_resid_maxabs, param_resid_maxabs)
        plt.xlabel(xaxis_label)
        if (i==0):
            plt.ylabel(yaxis_label)
        if (flip_y):
            plt.ylim(plt.ylim()[1], plt.ylim()[0])
        plt.title(param_names[i])

    if (param_residual_weights is not None):
        plt.sca(ax[param_residuals.shape[1]]) 
        #sum of the squared weighted resids
        sum_squared_weighted_resid = np.sum(np.square(param_residuals
                                            *param_residual_weights[None,:]),
                                            axis=1) 
        plt.scatter(x=xaxis_vals,
                    y=yaxis_vals,
                    c=param_residuals[:,i],
                    cmap="viridis")
        plt.colorbar()
        plt.xlabel(xaxis_label)
        if (i==0):
            plt.ylabel(yaxis_label)
        if (flip_y):
            plt.ylim(plt.ylim()[1], plt.ylim()[0])
        plt.title("Sum of squared\nweighted residuals")

    plt.show()


def plot_ompasoln_residuals(ompa_soln, xaxis_colname,
                            yaxis_colname, flip_y=True):
    plot_residuals(
        param_residuals=ompa_soln.param_residuals,
        param_names=ompa_soln.param_names,
        xaxis_vals=ompa_soln.obs_df[xaxis_colname],
        xaxis_label=xaxis_colname,
        yaxis_vals=ompa_soln.obs_df[yaxis_colname],
        yaxis_label=yaxis_colname, flip_y=flip_y,
        param_residuals_weights=ompa_soln.effective_param_weighting)


#deprecated now; api of ThermoclineArraySoln was updated such that can
# just call plot_ompasoln_endmember_fractions
def plot_thermocline_endmember_fractions(*args, **kwargs):
    return plot_ompasoln_endmember_fractions(*args, **kwargs)


#deprecated now; api of ThermoclineArraySoln was updated such that can
# just call plot_ompasoln_residuals
def plot_thermocline_residuals(*args, **kwargs):
    return plot_ompasoln_residuals(*args, **kwargs)


def nozero_xaxis(field_name):
  import altair as alt
  return alt.X(field_name, scale=alt.Scale(zero=False))


def nozero_yaxis(field_name, domain=None):
  import altair as alt
  if (domain is None):
    return alt.Y(field_name, scale=alt.Scale(zero=False))
  else:
    return alt.Y(field_name, scale=alt.Scale(zero=False, domain=domain))


def transect_scatterplot(basechart, selection,
                         property_name, altairdf, xaxis_colname, yaxis_colname,
                         zerocenter=False, flip_y=True):
    import altair as alt
    additional_color_kwargs = {}
    if (zerocenter):
        max_abs_property = np.max(np.abs(altairdf[property_name]))
        scale = alt.Scale(scheme='blueorange',
                          domain=[-max_abs_property, max_abs_property])
        additional_color_kwargs['scale'] = scale
    color = alt.condition(selection, property_name, alt.value('lightgray'),
                          title="",
                          **additional_color_kwargs)
    max_depth = np.max(altairdf[yaxis_colname])*1.05
    to_return = basechart.encode(
              nozero_xaxis(xaxis_colname),
              nozero_yaxis(yaxis_colname, domain=((max_depth,0) if flip_y
                                                  else (0,max_depth)) )
           ).encode(color=color).properties(title=property_name)
    return to_return
           

def wrap_scatterplots(scatterplots, resolve_scale='shared', rowsize=7):
    import altair as alt
    hconcats = [
        alt.hconcat(*scatterplots[i:i+rowsize]).resolve_scale(
            color=resolve_scale)
        for i in range(0,len(scatterplots), rowsize)  
    ]
    return alt.vconcat(*hconcats)


def pp_scatterplot(obs_basechart, selection,
                   endmember_basechart,
                   property1, property2, opacity):
    import altair as alt
    color = alt.condition(selection, alt.value('lightblue'),
                          alt.value('lightgray'))
    return (obs_basechart.mark_point(opacity=opacity).encode(
                nozero_xaxis(property1),
                nozero_yaxis(property2),
                color=color)
            + endmember_basechart.encode(nozero_xaxis(property1),
                                         nozero_yaxis(property2)))
    

def build_altair_viz(ompa_problem, xaxis_colname, yaxis_colname,
                     flip_y=True, chart_width=200, chart_height=200,
                     extra_tooltip_cols=[]):
    import altair as alt

    endmember_name_column = ompa_problem.endmember_name_column
    altairdf = pd.DataFrame(ompa_problem.obs_df)
    endmember_names = []
    for endmember_idx in range(ompa_problem.endmember_fractions.shape[1]):
        endmember_name =\
          ompa_problem.endmember_df[endmember_name_column][endmember_idx]
        endmember_names.append(endmember_name)
        altairdf[endmember_name] =\
            ompa_problem.endmember_fractions[:,endmember_idx]
    
    if (ompa_problem.total_oxygen_deficit is not None):
        altairdf["total O2 deficit"] = ompa_problem.total_oxygen_deficit
        for i in range(len(ompa_problem.converted_params_to_use)):
            converted_param_name = ompa_problem.converted_params_to_use[i]
            altairdf[converted_param_name+" eff conv ratio"] =\
                1.0/ompa_problem.effective_conversion_ratios[:,i]

    param_names = (ompa_problem.converted_params_to_use
                     +ompa_problem.conserved_params_to_use)
    for param_idx in range(ompa_problem.param_residuals.shape[1]):
        param_name = param_names[param_idx]
        altairdf[param_name+"_resid"] =\
          ompa_problem.param_residuals[:,param_idx]

    interval_selection = alt.selection_interval()
    tooltip_columns = (param_names
      + [xaxis_colname, yaxis_colname] + extra_tooltip_cols
      + endmember_names
      + [x+"_resid" for x in param_names]
      + (["total O2 deficit"]+
         [x+" eff conv ratio" for x in
          ompa_problem.converted_params_to_use]
         if ompa_problem.total_oxygen_deficit is not None
         else []))
    #make the linked property-property plots
    obs_basechart = alt.Chart(altairdf).mark_point().encode(
      tooltip=tooltip_columns,
      color=alt.condition(interval_selection,
                          alt.value("lightblue"),#"NPIW",
                          alt.value('lightgray'))
      ).add_selection(interval_selection).properties(
          width=chart_width,
          height=chart_height)
    
    endmember_basechart =\
      alt.Chart(ompa_problem.endmember_df).mark_point(
          shape="diamond", size=50).encode(
              color=endmember_name_column).properties(
                width=chart_width,
                height=chart_height)

    #display a row that is the endmember fractions
    endmember_fraction_scatterplots = [
        transect_scatterplot(
            basechart=obs_basechart,
            selection=interval_selection,
            property_name=property_name,
            altairdf=altairdf,
            xaxis_colname=xaxis_colname, yaxis_colname=yaxis_colname,
            flip_y=flip_y)
        for property_name in endmember_names[:-1]]

    if (ompa_problem.total_oxygen_deficit is not None):
        oxygen_deficit_scatterplot = transect_scatterplot(
                basechart=obs_basechart,
                selection=interval_selection,
                property_name="total O2 deficit",
                altairdf=altairdf,
                xaxis_colname=xaxis_colname, yaxis_colname=yaxis_colname,
                flip_y=flip_y)

    the_pp_scatterplots = []
    for i in range(len(param_names)):
        for j in range(i+1,len(param_names)):
            the_pp_scatterplots.append(
                pp_scatterplot(
                    obs_basechart=obs_basechart,
                    selection=interval_selection,
                    endmember_basechart=endmember_basechart,
                    property1=param_names[i],
                    property2=param_names[j],
                    opacity=0.2)
            )
    
    resid_scatterplots = [
        transect_scatterplot(
            basechart=obs_basechart,
            selection=interval_selection,
            property_name=property_name+"_resid",
            altairdf=altairdf,
            zerocenter=True,
            xaxis_colname=xaxis_colname, yaxis_colname=yaxis_colname,
            flip_y=flip_y)
        for property_name in param_names]

    prop_scatterplots = [
        transect_scatterplot(
            basechart=obs_basechart,
            selection=interval_selection,
            property_name=property_name,
            altairdf=altairdf,
            zerocenter=False,
            xaxis_colname=xaxis_colname, yaxis_colname=yaxis_colname,
            flip_y=flip_y)
        for property_name in param_names]

    conversion_ratio_scatterplots = [
      transect_scatterplot(
          basechart=obs_basechart,
            selection=interval_selection,
            property_name=paramname+" eff conv ratio",
            altairdf=altairdf,
            zerocenter=False,
            xaxis_colname=xaxis_colname, yaxis_colname=yaxis_colname,
            flip_y=flip_y)
      for paramname in ompa_problem.converted_params_to_use
    ]
    #assert False

    return alt.vconcat(*(
        
              [wrap_scatterplots(endmember_fraction_scatterplots)]
              +([oxygen_deficit_scatterplot,
                 wrap_scatterplots(conversion_ratio_scatterplots,
                                   resolve_scale="independent")]
                 if ompa_problem.total_oxygen_deficit is not None else [])
              +[wrap_scatterplots(prop_scatterplots,
                                  resolve_scale="independent"),
                wrap_scatterplots(resid_scatterplots,
                                  resolve_scale='independent'),
                wrap_scatterplots(the_pp_scatterplots)]
              
              ))


def build_thermocline_altair_viz(ompa_problems_arr,
                                 xaxis_colname, yaxis_colname,
                                 chart_width=200, chart_height=200,
                                 flip_y=True, extra_tooltip_cols=[]):
    import altair as alt
    #verify endmember names are the same for all
    endmember_name_column = ompa_problems_arr[0].endmember_name_column
    endmember_names = tuple(ompa_problems_arr[0].endmember_df[
                             endmember_name_column])
    assert all(tuple(x.endmember_df[endmember_name_column])==endmember_names
               for x in ompa_problems_arr)
    #similar verification for param_names
    param_names = tuple(ompa_problems_arr[0].conserved_params_to_use
                        + ompa_problems_arr[0].converted_params_to_use)
    assert all(tuple(x.conserved_params_to_use
                     + x.converted_params_to_use)==param_names
               for x in ompa_problems_arr)

    #concatenate all the observations to get a new obs_df
    altairdf = pd.concat([x.obs_df for x in ompa_problems_arr])

    for endmember_idx in range(ompa_problems_arr[0]
                                .endmember_fractions.shape[1]):
      endmember_name = endmember_names[endmember_idx]
      altairdf[endmember_name] = np.concatenate([
        x.endmember_fractions[:,endmember_idx]
        for x in ompa_problems_arr])
      
    if (ompa_problems_arr[0].total_oxygen_deficit is not None):
        altairdf["total O2 deficit"] = np.concatenate([
          x.total_oxygen_deficit
          for x in ompa_problems_arr])
        for i in range(len(ompa_problems_arr[0].converted_params_to_use)):
            converted_param_name =\
              ompa_problems_arr[0].converted_params_to_use[i]
            effective_conversion_ratios = np.concatenate([
              x.effective_conversion_ratios for x in ompa_problems_arr])
            altairdf[converted_param_name+" eff conv ratio"] =\
                1.0/effective_conversion_ratios[:,i]

    for param_idx in range(ompa_problems_arr[0].param_residuals.shape[1]):
        param_name = param_names[param_idx]
        altairdf[param_name+"_resid"] = np.concatenate([
          x.param_residuals[:,param_idx] for x in ompa_problems_arr])          

    param_names = (ompa_problems_arr[0].conserved_params_to_use
                     +ompa_problems_arr[0].converted_params_to_use)
    interval_selection = alt.selection_interval()
    tooltip_columns = (list(param_names)
      + [xaxis_colname, yaxis_colname] + extra_tooltip_cols
      + list(endmember_names)
      + [x+"_resid" for x in param_names]
      + (["total O2 deficit"]+
         [x+" eff conv ratio" for x in
          ompa_problems_arr[0].converted_params_to_use]
         if ompa_problems_arr[0].total_oxygen_deficit is not None
         else []) )
    
    #make the linked property-property plots
    obs_basechart = alt.Chart(altairdf).mark_point().encode(
      tooltip=tooltip_columns,
      color=alt.condition(interval_selection,
                          alt.value("lightblue"),#"NPIW",
                          alt.value('lightgray'))
      ).add_selection(interval_selection).properties(
          width=chart_width,
          height=chart_height)
    
    #prepare a endmember df
    endmember_df = pd.concat([x.endmember_df for x in ompa_problems_arr])

    endmember_basechart =\
      alt.Chart(endmember_df).mark_point(
          shape="diamond", size=50).encode(
              tooltip=list(param_names),
              color=endmember_name_column).properties(
                width=chart_width,
                height=chart_height)

    #display a row that is the endmember fractions
    endmember_fraction_scatterplots = [
        transect_scatterplot(
            basechart=obs_basechart,
            selection=interval_selection,
            property_name=property_name,
            altairdf=altairdf,
            xaxis_colname=xaxis_colname, yaxis_colname=yaxis_colname,
            flip_y=flip_y)
        for property_name in endmember_names]

    if (ompa_problems_arr[0].total_oxygen_deficit is not None):
        oxygen_deficit_scatterplot = transect_scatterplot(
                basechart=obs_basechart,
                selection=interval_selection,
                property_name="total O2 deficit",
                altairdf=altairdf,
                xaxis_colname=xaxis_colname, yaxis_colname=yaxis_colname,
                flip_y=flip_y)

    the_pp_scatterplots = []
    for i in range(len(param_names)):
        for j in range(i+1,len(param_names)):
            the_pp_scatterplots.append(
                pp_scatterplot(
                    obs_basechart=obs_basechart,
                    selection=interval_selection,
                    endmember_basechart=endmember_basechart,
                    property1=param_names[i],
                    property2=param_names[j],
                    opacity=1)
            )
    
    resid_scatterplots = [
        transect_scatterplot(
            basechart=obs_basechart,
            selection=interval_selection,
            property_name=property_name+"_resid",
            altairdf=altairdf,
            zerocenter=True,
            xaxis_colname=xaxis_colname, yaxis_colname=yaxis_colname,
            flip_y=flip_y)
        for property_name in param_names]

    prop_scatterplots = [
        transect_scatterplot(
            basechart=obs_basechart,
            selection=interval_selection,
            property_name=property_name,
            altairdf=altairdf,
            zerocenter=True,
            xaxis_colname=xaxis_colname, yaxis_colname=yaxis_colname,
            flip_y=flip_y)
        for property_name in param_names]

    conversion_ratio_scatterplots = [
      transect_scatterplot(
          basechart=obs_basechart,
            selection=interval_selection,
            property_name=paramname+" eff conv ratio",
            altairdf=altairdf,
            zerocenter=False,
            xaxis_colname=xaxis_colname, yaxis_colname=yaxis_colname,
            flip_y=flip_y)
      for paramname in ompa_problems_arr[0].converted_params_to_use
    ]

    return alt.vconcat(*(
      [wrap_scatterplots(endmember_fraction_scatterplots)]
      +([oxygen_deficit_scatterplot,
        wrap_scatterplots(conversion_ratio_scatterplots,
                          resolve_scale="independent")]
        if ompa_problems_arr[0].total_oxygen_deficit is not None else [])
      +[wrap_scatterplots(prop_scatterplots,
                          resolve_scale="independent"),
        wrap_scatterplots(resid_scatterplots,
                          resolve_scale='independent'),
        wrap_scatterplots(the_pp_scatterplots)]      
      ))
