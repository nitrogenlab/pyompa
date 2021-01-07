from __future__ import division, print_function
from matplotlib import pyplot as plt
import numpy as np


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
        plt.colorbar()
        plt.xlabel("latitude")
        plt.title("oxygen deficit")
    for i in range(len(converted_params_to_use)):
        plt.sca(ax[num_watermasses+1+i])
        plt.scatter(latitudes, depths,
                    c=1.0/effective_conversion_ratios[:,i])
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
