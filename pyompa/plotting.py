from __future__ import division, print_function
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import altair as alt


def plot_water_mass_fractions(latitudes, depths,
    water_mass_fractions, watermassnames, total_oxygen_deficit,
    converted_params_to_use, effective_conversion_ratios):
    num_watermasses = water_mass_fractions.shape[1]
    num_figs = (num_watermasses +
                (1+len(converted_params_to_use)
                 if total_oxygen_deficit is not None else 0))
    print("numfigs:", num_figs)
    fig, ax = plt.subplots(nrows=1, ncols=num_figs, figsize=(5*num_figs,4))
    for i in range(num_watermasses):
        plt.sca(ax[i])
        plt.scatter(latitudes, depths, c=water_mass_fractions[:,i])
        plt.xlabel("latitude")
        if (i==0):
            plt.ylabel("depth")
        plt.ylim(plt.ylim()[1], plt.ylim()[0])
        plt.colorbar()
        plt.title(watermassnames[i])
    if (total_oxygen_deficit is not None):
        plt.sca(ax[num_watermasses])
        plt.scatter(latitudes, depths, c=total_oxygen_deficit)
        plt.ylim(plt.ylim()[1], plt.ylim()[0])
        plt.colorbar()
        plt.xlabel("latitude")
        plt.title("oxygen deficit")
    for i in range(len(converted_params_to_use)):
        plt.sca(ax[num_watermasses+1+i])
        plt.scatter(latitudes, depths,
                    c=1.0/effective_conversion_ratios[:,i])
        plt.ylim(plt.ylim()[1], plt.ylim()[0])
        plt.colorbar()
        plt.xlabel("latitude")
        plt.title(converted_params_to_use[i]
                  +" \nconversion ratio (relative to oxygen)")
    plt.show()


def plot_residuals(param_residuals, params_to_use, latitudes, depths):
    num_params = param_residuals.shape[1]
    fig, ax = plt.subplots(nrows=1, ncols=num_params, figsize=(5*num_params,4))
    for i in range(param_residuals.shape[1]):
        plt.sca(ax[i])
        plt.scatter(x=latitudes,
                    y=depths,
                    c=np.abs(param_residuals[:,i]),
                    cmap="viridis")
        plt.colorbar()
        plt.xlabel("latitude")
        if (i==0):
            plt.ylabel("depth")
        plt.ylim(plt.ylim()[1], plt.ylim()[0])
        plt.title(params_to_use[i])
    plt.show()


def plot_thermocline_water_mass_fractions(ompa_problems_arr):
    #first, check to make sure all entries in ompa_problems_arr have the same
    # number of water masses; the assert statement will throw an error if
    # that is not the case

    num_watermasses = ompa_problems_arr[0].water_mass_fractions.shape[1]
    assert all([x.water_mass_fractions.shape[1]==num_watermasses
                for x in ompa_problems_arr])

    latitudes = np.concatenate([
            np.array(x.obs_df["latitude"]) for x in ompa_problems_arr])
    depths = np.concatenate([
        np.array(x.obs_df["depth"]) for x in ompa_problems_arr])
    water_mass_fractions = np.concatenate([
            x.water_mass_fractions for x in ompa_problems_arr], axis=0)
    watermassnames = list(ompa_problems_arr[0].watermass_df["watermassname"])
    converted_params_to_use = ompa_problems_arr[0].converted_params_to_use
    if (ompa_problems_arr[0].total_oxygen_deficit is not None):
        total_oxygen_deficit = np.concatenate([x.total_oxygen_deficit
                                           for x in ompa_problems_arr], axis=0)
        effective_conversion_ratios = np.concatenate([
              x.effective_conversion_ratios for x in ompa_problems_arr], axis=0)
    else:
        total_oxygen_deficit = None
        effective_conversion_ratios = None

    plot_water_mass_fractions(
        latitudes=latitudes,
        depths=depths,
        water_mass_fractions=water_mass_fractions,
        watermassnames=watermassnames,
        total_oxygen_deficit=total_oxygen_deficit,
        converted_params_to_use=converted_params_to_use,
        effective_conversion_ratios=effective_conversion_ratios)


def plot_thermocline_residuals(ompa_problems_arr):

    param_residuals = np.concatenate([
            x.param_residuals for x in ompa_problems_arr], axis=0)
    params_to_use = (ompa_problems_arr[0].conserved_params_to_use
                     +ompa_problems_arr[0].converted_params_to_use)
    latitudes = np.concatenate([
            np.array(x.obs_df["latitude"]) for x in ompa_problems_arr])
    depths = np.concatenate([
        np.array(x.obs_df["depth"]) for x in ompa_problems_arr])

    plot_residuals(
        param_residuals=param_residuals,
        params_to_use=params_to_use,
        latitudes=latitudes,
        depths=depths)


def nozero_xaxis(field_name):
  return alt.X(field_name, scale=alt.Scale(zero=False))


def nozero_yaxis(field_name, domain=None):
  if (domain is None):
    return alt.Y(field_name, scale=alt.Scale(zero=False))
  else:
    return alt.Y(field_name, scale=alt.Scale(zero=False, domain=domain))


def latdepth_scatterplot(basechart, selection,
                         property_name, altairdf, zerocenter=False):
    additional_color_kwargs = {}
    if (zerocenter):
        max_abs_property = np.max(np.abs(altairdf[property_name]))
        scale = alt.Scale(scheme='blueorange',
                          domain=[-max_abs_property, max_abs_property])
        additional_color_kwargs['scale'] = scale
    color = alt.condition(selection, property_name, alt.value('lightgray'),
                          title="",
                          **additional_color_kwargs)
    max_depth = np.max(altairdf["depth"])*1.05
    to_return = basechart.encode(
              nozero_xaxis("latitude"),
              nozero_yaxis("depth", domain=(max_depth,0))
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
                   watermass_basechart,
                   property1, property2, opacity):
    color = alt.condition(selection, alt.value('lightblue'),
                          alt.value('lightgray'))
    return (obs_basechart.mark_point(opacity=opacity).encode(
                nozero_xaxis(property1),
                nozero_yaxis(property2),
                color=color)
            + watermass_basechart.encode(nozero_xaxis(property1),
                                         nozero_yaxis(property2)))
    

def build_altair_viz(ompa_problem, chart_width=200, chart_height=200):
    import altair as alt

    altairdf = pd.DataFrame(ompa_problem.obs_df)
    watermass_names = []
    for water_mass_idx in range(ompa_problem.water_mass_fractions.shape[1]):
        watermass_name =\
          ompa_problem.watermass_df["watermassname"][water_mass_idx]
        watermass_names.append(watermass_name)
        altairdf[watermass_name] =\
            ompa_problem.water_mass_fractions[:,water_mass_idx]
    
    if (ompa_problem.total_oxygen_deficit is not None):
        altairdf["total O2 deficit"] = ompa_problem.total_oxygen_deficit
        for i in range(len(ompa_problem.converted_params_to_use)):
            converted_param_name = ompa_problem.converted_params_to_use[i]
            altairdf[converted_param_name+" eff conv ratio"] =\
                1.0/ompa_problem.effective_conversion_ratios[:,i]

    params_to_use = (ompa_problem.converted_params_to_use
                     +ompa_problem.conserved_params_to_use)
    for param_idx in range(ompa_problem.param_residuals.shape[1]):
        param_name = params_to_use[param_idx]
        altairdf[param_name+"_resid"] =\
          ompa_problem.param_residuals[:,param_idx]

    interval_selection = alt.selection_interval()
    tooltip_columns = (params_to_use
      + ["latitude", "longitude"]
      + watermass_names
      + [x+"_resid" for x in params_to_use]
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
    
    watermass_basechart =\
      alt.Chart(ompa_problem.watermass_df).mark_point(
          shape="diamond", size=50).encode(
              color="watermassname").properties(
                width=chart_width,
                height=chart_height)

    #display a row that is the water mass fractions
    watermass_fraction_scatterplots = [
        latdepth_scatterplot(
            basechart=obs_basechart,
            selection=interval_selection,
            property_name=property_name,
            altairdf=altairdf)
        for property_name in watermass_names[:-1]]

    if (ompa_problem.total_oxygen_deficit is not None):
        oxygen_deficit_scatterplot = latdepth_scatterplot(
                basechart=obs_basechart,
                selection=interval_selection,
                property_name="total O2 deficit",
                altairdf=altairdf)

    the_pp_scatterplots = []
    for i in range(len(params_to_use)):
        for j in range(i+1,len(params_to_use)):
            the_pp_scatterplots.append(
                pp_scatterplot(
                    obs_basechart=obs_basechart,
                    selection=interval_selection,
                    watermass_basechart=watermass_basechart,
                    property1=params_to_use[i],
                    property2=params_to_use[j],
                    opacity=0.2)
            )
    
    resid_scatterplots = [
        latdepth_scatterplot(
            basechart=obs_basechart,
            selection=interval_selection,
            property_name=property_name+"_resid",
            altairdf=altairdf,
            zerocenter=True)
        for property_name in params_to_use]

    prop_scatterplots = [
        latdepth_scatterplot(
            basechart=obs_basechart,
            selection=interval_selection,
            property_name=property_name,
            altairdf=altairdf,
            zerocenter=False)
        for property_name in params_to_use]

    conversion_ratio_scatterplots = [
      latdepth_scatterplot(
          basechart=obs_basechart,
            selection=interval_selection,
            property_name=paramname+" eff conv ratio",
            altairdf=altairdf,
            zerocenter=False)
      for paramname in ompa_problem.converted_params_to_use
    ]
    #assert False

    return alt.vconcat(*(
        
              [wrap_scatterplots(watermass_fraction_scatterplots)]
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
                                 chart_width=200, chart_height=200):

    #verify watermass names are the same for all
    watermass_names = tuple(ompa_problems_arr[0].watermass_df["watermassname"])
    assert all(tuple(x.watermass_df["watermassname"])==watermass_names
               for x in ompa_problems_arr)
    #similar verification for param_names
    param_names = tuple(ompa_problems_arr[0].conserved_params_to_use
                        + ompa_problems_arr[0].converted_params_to_use)
    assert all(tuple(x.conserved_params_to_use
                     + x.converted_params_to_use)==param_names
               for x in ompa_problems_arr)

    #concatenate all the observations to get a new obs_df
    altairdf = pd.concat([x.obs_df for x in ompa_problems_arr])

    for water_mass_idx in range(ompa_problems_arr[0]
                                .water_mass_fractions.shape[1]):
      watermass_name = watermass_names[water_mass_idx]
      altairdf[watermass_name] = np.concatenate([
        x.water_mass_fractions[:,water_mass_idx]
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

    params_to_use = (ompa_problems_arr[0].conserved_params_to_use
                     +ompa_problems_arr[0].converted_params_to_use)
    interval_selection = alt.selection_interval()
    tooltip_columns = (list(param_names)
      + ["latitude", "longitude"]
      + list(watermass_names)
      + [x+"_resid" for x in params_to_use]
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
    
    #prepare a watermass df
    watermass_df = pd.concat([x.watermass_df for x in ompa_problems_arr])

    watermass_basechart =\
      alt.Chart(watermass_df).mark_point(
          shape="diamond", size=50).encode(
              tooltip=list(param_names),
              color="watermassname").properties(
                width=chart_width,
                height=chart_height)

    #display a row that is the water mass fractions
    watermass_fraction_scatterplots = [
        latdepth_scatterplot(
            basechart=obs_basechart,
            selection=interval_selection,
            property_name=property_name,
            altairdf=altairdf)
        for property_name in watermass_names]

    if (ompa_problems_arr[0].total_oxygen_deficit is not None):
        oxygen_deficit_scatterplot = latdepth_scatterplot(
                basechart=obs_basechart,
                selection=interval_selection,
                property_name="total O2 deficit",
                altairdf=altairdf)

    the_pp_scatterplots = []
    for i in range(len(params_to_use)):
        for j in range(i+1,len(params_to_use)):
            the_pp_scatterplots.append(
                pp_scatterplot(
                    obs_basechart=obs_basechart,
                    selection=interval_selection,
                    watermass_basechart=watermass_basechart,
                    property1=params_to_use[i],
                    property2=params_to_use[j],
                    opacity=1)
            )
    
    resid_scatterplots = [
        latdepth_scatterplot(
            basechart=obs_basechart,
            selection=interval_selection,
            property_name=property_name+"_resid",
            altairdf=altairdf,
            zerocenter=True)
        for property_name in params_to_use]

    prop_scatterplots = [
        latdepth_scatterplot(
            basechart=obs_basechart,
            selection=interval_selection,
            property_name=property_name,
            altairdf=altairdf,
            zerocenter=True)
        for property_name in params_to_use]

    conversion_ratio_scatterplots = [
      latdepth_scatterplot(
          basechart=obs_basechart,
            selection=interval_selection,
            property_name=paramname+" eff conv ratio",
            altairdf=altairdf,
            zerocenter=False)
      for paramname in ompa_problems_arr[0].converted_params_to_use
    ]

    return alt.vconcat(*(
      [wrap_scatterplots(watermass_fraction_scatterplots)]
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
