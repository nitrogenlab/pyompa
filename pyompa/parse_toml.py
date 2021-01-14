from __future__ import division
import toml
import pandas as pd
from .ompacore import OMPAProblem
from .endmemberpenaltyfunc import 
from .util import assert_in, assert_compatible_keys, assert_has_keys
from collections import OrderedDict


PARSE_DF_ALLOWED_KEYS = ["csv_file", "na_values"]


def parse_df_from_config(config, config_file_type):
    assert "csv_file" in config,\
        "Need argument 'csv_file' when parsing "+config_file_type+" config"
    assert_has_keys(the_dict=config,
          allowed=["csv_file"],
          errorprefix="Issue when parsing "+config_file_type+" config: ")  
    kwargs = {}
    if ("na_values" in config):
        config["na_values"] = config["na_values"] 
    df = pd.read_csv(config["csv_file"], **kwargs)
    return df


def parse_observations_config(config):
    assert_compatible_keys(the_dict=config, allowed=PARSE_DF_ALLOWED_KEYS,
         errorprefix="Issue when parsing observations config: ")  
    return parse_df_from_config(config=config, config_file_type="obervations")


def parse_endmembers_config(config):
    assert_compatible_keys(the_dict=config,
          allowed=PARSE_DF_ALLOWED_KEYS+["endmember_name_column"],
          errorprefix="Issue when parsing endmembers config: ")  
    assert_has_keys(the_dict=config,
          allowed=["endmembers"],
          errorprefix="Issue when parsing endmembers config: ")  
    endmembers_df =\
        parse_df_from_config(config=config, config_file_type="endmembers")
    endmember_name_column = config["endmember_name_column"]

    return endmembers_df, endmember_name_column


def parse_endmember_penalty_from_config(config):
    return OrderedDict([(endmember_name, EndMemExpPenaltyFunc(subconfig))
                        for endmember_name, subconfig in config.items()]) 


def parse_params(config):
    PARSE_PARAMS_ALLOWED_KEYS = ["weight", "remineralized", "ratios"] 
    paramsandweighting_conserved = []
    paramsandweighting_converted = []
    conversionratios = {}
    for param_name in config:
        param_config = config[paramname]
        assert_compatible_keys(the_dict=param_config,
            allowed=PARSE_PARAMS_ALLOWED_KEYS,
            errorprefix="Issue when parsing param config for "+param_name+": ") 
        assert_has_keys(the_dict=param_config,
            allowed=["weight", "remineralized"],
            errorprefix="Issue when parsing param config for "+param_name+": ")
        weight = param_config["weight"]
        remineralized = param_config["remineralized"]
        if (remineralized == True):
            paramsandweighting_converted.append((param_name, weight))
            assert "ratios" in param_config, ("'ratios' must be specified for "
              +"param "+param_name+" if remineralized=true")
            ratios = param_config["ratios"]
            conversion_ratios[param_name] = ratios 
        else:
            paramsandweighting_conserved.append((param_name, weight))
            assert "ratios" not in param_config, ("'ratios' is only applicable "
              +"for param "+param_name+" if remineralized=true")
    return (paramsandweighting_conserved, paramsandweighting_converted,
            conversionratios)


def run_analysis_from_toml_file(toml_config_file):
    config = toml.loads(open(toml_config_file).read())
    assert_compatible_keys(the_dict=config,
          allowed=["observations", "params", "endmembers",
                   "endmember_penalties"],
          errorprefix="Issue when parsing toml config: ")  
    assert_has_keys(the_dict=config,
          allowed=["observations", "params", "endmembers"],
          errorprefix="Issue when parsing toml config: ")  

    obs_df = parse_observations_config(config=config["observations"]) 

    (endmember_df, endmember_name_column) = (
        parse_observations_config(config=config["endmembers"]))

    (paramsandweighting_conserved, paramsandweighting_converted,
     conversionratios) = parse_params(config=config["params"])

    if "endmember_penalties" in config:
        endmembername_to_usagepenaltyfunc = (
            parse_endmember_penalty_from_config(
                config=config["endmember_penalties"]))
    else:
        endmembername_to_usagepenaltyfunc = {}

    ompa_soln = OMPAProblem(
          obs_df = obs_df,
          paramsandweighting_conserved=paramsandweighting_conserved,
          paramsandweighting_converted=paramsandweighting_converted,
          conversionratios=conversionratios,
          smoothness_lambda=None,
          endmembername_to_usagepenaltyfunc=
            endmembername_to_usagepenaltyfunc).solve(
              endmember_df=endmember_df,
              endmember_name_column="endmember_name") 

    return ompa_soln
