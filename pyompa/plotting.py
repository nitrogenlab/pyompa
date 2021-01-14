from __future__ import division, print_function
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import altair as alt


def plot_endmember_usagepenalties(endmembername_to_usagepenalty,
                                  xaxis_vals, xaxis_label,
                                  yaxis_vals, yaxis_label,
                                  flip_y=True):
    endmembernames = sorted(endmembername_to_usagepenalty)
    num_figs = len(endmembernames)
    fig, ax = plt.subplots(nrows=1, ncols=num_figs, figsize=(5*num_figs,4))
    for i in range(endmembernames):
        plt.sca(ax[i])
        plt.scatter(xaxis_vals, yaxis_vals, c=penalty)
        plt.xlabel(xaxis_label)
        if (i==0):
            plt.ylabel(yaxis_label)
        if (flip_y):
            plt.ylim(plt.ylim()[1], plt.ylim()[0])
        plt.colorbar()
        plt.title(endmembernames[i])
    plt.show()


def plot_ompasoln_endmember_usagepenalties(
        ompa_soln, endmembername_to_usagepenalty,
        xaxis_colname, yaxis_colname, flip_y):
    plot_endmember_usagepenalties(
        endmembername_to_usagepenalty=
            ompa_soln.ompa_problem.endmembername_to_usagepenalty,
        xaxis_vals=ompa_soln.obs_df[xaxis_colname],
        xaxis_label=xaxis_colname,
        yaxis_vals=ompa_soln.obs_df[yaxis_colname],
        yaxis_label=yaxis_colname,
        flip_y=True)


def plot_endmember_fractions(xaxis_vals, xaxis_label, yaxis_vals, yaxis_label,
        endmember_fractions, endmembernames, total_oxygen_deficit,
        effective_conversion_ratios, converted_param_names,
        flip_y=True):
    num_endmembers = endmember_fractions.shape[1]
    num_figs = (num_endmembers +
                (1+len(converted_param_names)
                 if total_oxygen_deficit is not None else 0))
    fig, ax = plt.subplots(nrows=1, ncols=num_figs, figsize=(5*num_figs,4))
    for i in range(num_endmembers):
        plt.sca(ax[i])
        plt.scatter(xaxis_vals, yaxis_vals, c=endmember_fractions[:,i])
        plt.xlabel(xaxis_label)
        if (i==0):
            plt.ylabel(yaxis_label)
        if (flip_y):
            plt.ylim(plt.ylim()[1], plt.ylim()[0])
        plt.colorbar()
        plt.clim(0.0, 1.0)
        plt.title(endmembernames[i])
    if (total_oxygen_deficit is not None):
        plt.sca(ax[num_endmembers])
        plt.scatter(xaxis_vals, yaxis_vals, c=total_oxygen_deficit)
        if (flip_y):
            plt.ylim(plt.ylim()[1], plt.ylim()[0])
        plt.colorbar()
        plt.xlabel(xaxis_label)
        plt.title("oxygen deficit")
    for i in range(len(converted_param_names)):
        plt.sca(ax[num_endmembers+1+i])
        plt.scatter(xaxis_vals, yaxis_vals,
                    c=1.0/effective_conversion_ratios[:,i])
        if (flip_y):
            plt.ylim(plt.ylim()[1], plt.ylim()[0])
        plt.colorbar()
        plt.xlabel(xaxis_label)
        plt.title(converted_param_names[i]
                  +" \nconversion ratio (relative to oxygen)")
    plt.show()


def plot_ompasoln_endmember_fractions(ompa_soln, xaxis_colname,
                                      yaxis_colname, flip_y=True):
    plot_endmember_fractions(
        xaxis_vals=ompa_soln.obs_df[xaxis_colname],
        xaxis_label=xaxis_colname,
        yaxis_vals=ompa_soln.obs_df[yaxis_colname],
        yaxis_label=yaxis_colname,
        endmember_fractions=ompa_soln.endmember_fractions,
        endmembernames=list(
            ompa_soln.endmember_df[ompa_soln.endmember_name_column]),
        total_oxygen_deficit=ompa_soln.total_oxygen_deficit,
        effective_conversion_ratios=ompa_soln.effective_conversion_ratios,
        converted_param_names=ompa_soln.converted_params_to_use,
        flip_y=flip_y)


def plot_residuals(param_residuals, param_names, xaxis_vals, xaxis_label,
                   yaxis_vals, yaxis_label, flip_y=True):
    num_params = param_residuals.shape[1]
    fig, ax = plt.subplots(nrows=1, ncols=num_params, figsize=(5*num_params,4))
    for i in range(param_residuals.shape[1]):
        plt.sca(ax[i])
        plt.scatter(x=xaxis_vals,
                    y=yaxis_vals,
                    c=np.abs(param_residuals[:,i]),
                    cmap="viridis")
        plt.colorbar()
        plt.xlabel(xaxis_label)
        if (i==0):
            plt.ylabel(yaxis_label)
        if (flip_y):
            plt.ylim(plt.ylim()[1], plt.ylim()[0])
        plt.title(param_names[i])
    plt.show()


def plot_ompasoln_residuals(ompa_soln, xaxis_colname,
                            yaxis_colname, flip_y=True):
    plot_residuals(
        param_residuals=ompa_soln.param_residuals,
        param_names=(ompa_soln.conserved_params_to_use
                      + ompa_soln.converted_params_to_use),
        xaxis_vals=ompa_soln.obs_df[xaxis_colname],
        xaxis_label=xaxis_colname,
        yaxis_vals=ompa_soln.obs_df[yaxis_colname],
        yaxis_label=yaxis_colname, flip_y=flip_y)


def plot_thermocline_endmember_fractions(ompa_problems_arr,
                                         xaxis_colname, yaxis_colname,
                                         flip_y=True):
    #first, check to make sure all entries in ompa_problems_arr have the same
    # number of endmembers; the assert statement will throw an error if
    # that is not the case

    endmember_name_column = ompa_problems_arr[0].endmember_name_column
    num_endmembers = ompa_problems_arr[0].endmember_fractions.shape[1]
    assert all([x.endmember_fractions.shape[1]==num_endmembers
                for x in ompa_problems_arr])

    xaxis_vals = np.concatenate([
            np.array(x.obs_df[xaxis_colname]) for x in ompa_problems_arr])
    yaxis_vals = np.concatenate([
        np.array(x.obs_df[yaxis_colname]) for x in ompa_problems_arr])
    endmember_fractions = np.concatenate([
            x.endmember_fractions for x in ompa_problems_arr], axis=0)
    endmembernames = list(ompa_problems_arr[0].endmember_df[
                            endmember_name_column])
    converted_param_names = ompa_problems_arr[0].converted_params_to_use
    if (ompa_problems_arr[0].total_oxygen_deficit is not None):
        total_oxygen_deficit = np.concatenate([x.total_oxygen_deficit
                                           for x in ompa_problems_arr], axis=0)
        effective_conversion_ratios = np.concatenate([
              x.effective_conversion_ratios for x in ompa_problems_arr], axis=0)
    else:
        total_oxygen_deficit = None
        effective_conversion_ratios = None

    plot_endmember_fractions(
        xaxis_vals=xaxis_vals,
        yaxis_vals=yaxis_vals,
        xaxis_label=xaxis_colname,
        yaxis_label=yaxis_colname,
        endmember_fractions=endmember_fractions,
        endmembernames=endmembernames,
        total_oxygen_deficit=total_oxygen_deficit,
        converted_param_names=converted_param_names,
        effective_conversion_ratios=effective_conversion_ratios,
        flip_y=flip_y)


def plot_thermocline_residuals(ompa_problems_arr, xaxis_colname, yaxis_colname,
                               flip_y=True):

    param_residuals = np.concatenate([
            x.param_residuals for x in ompa_problems_arr], axis=0)
    param_names = (ompa_problems_arr[0].conserved_params_to_use
                     +ompa_problems_arr[0].converted_params_to_use)
    xaxis_vals = np.concatenate([
            np.array(x.obs_df[xaxis_colname]) for x in ompa_problems_arr])
    yaxis_vals = np.concatenate([
        np.array(x.obs_df[yaxis_colname]) for x in ompa_problems_arr])

    plot_residuals(
        param_residuals=param_residuals,
        param_names=param_names,
        xaxis_vals=xaxis_vals,
        yaxis_vals=yaxis_vals,
        xaxis_label=xaxis_colname,
        yaxis_label=yaxis_colname,
        flip_y=flip_y)


def nozero_xaxis(field_name):
  return alt.X(field_name, scale=alt.Scale(zero=False))


def nozero_yaxis(field_name, domain=None):
  if (domain is None):
    return alt.Y(field_name, scale=alt.Scale(zero=False))
  else:
    return alt.Y(field_name, scale=alt.Scale(zero=False, domain=domain))


def transect_scatterplot(basechart, selection,
                         property_name, altairdf, xaxis_colname, yaxis_colname,
                         zerocenter=False, flip_y=True):
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
    hconcats = [
        alt.hconcat(*scatterplots[i:i+rowsize]).resolve_scale(
            color=resolve_scale)
        for i in range(0,len(scatterplots), rowsize)  
    ]
    return alt.vconcat(*hconcats)


def pp_scatterplot(obs_basechart, selection,
                   endmember_basechart,
                   property1, property2, opacity):
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
