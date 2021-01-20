# pyompa
This is a python package for OMPA analysis. It is currently under development and has not been formally published yet. If you are interested in using this package, please contact Professor Karen Casciotti (kcasciot [at] stanford [dot] edu).

## Installation

pyompa is on pypi, and can be installed using the following pip command:
```bash
pip install pyompa
```
pip comes installed with Python 2 >=2.7.9 or Python 3 >=3.4. If you find that pip is not installed, refer to [these instructions](https://pip.pypa.io/en/stable/installing/).

## Running pyompa from the command line

After installing pyompa, it can be run with a single command:
```bash
run_ompa_given_config path/to/the/config/file.txt
```
The configuration file specifies all the settings about where to find the observations .csv, the end-members .csv, which parameters to use, etc. An example configuration file is pasted below:

```bash
#Comments can be added to the configuration file using the hashtag symbol (#).
# Anything following the hashtag is ignored. FYI, the 'type' of configuration
# file used by pyompa is a TOML configuration file (don't worry if you don't
# know what that means). The main thing to note is that sections in this type
# of configuration file are denoted using [section name].

[observations] #this specifies that the next section is going to talk about
#               settings pertaining to the observations file.
# csv_file specifies the path to a .csv file containing the observations. For
# an example, refer to
# https://github.com/nitrogenlab/GP15_watermassanalysis/blob/988fb74f44b7c62b32c6a11752a9bcbb06471783/gp15_intermediateanddeep_obs.csv
#The csv_file can have any columns you want, but it must contain column headers
# corresponding to the parameters you wish to use in the OMPA analysis.
csv_file = "gp15_intermediateanddeep_obs.csv"
na_values = -999 #this specifies how "NA" values are denoted in the
#                 observations .csv file

[endmembers]
#Like the observations file, the endmembers file must have column headers
# corresponding to the parameters you wish to use in theOMPA analysis. See
#https://github.com/nitrogenlab/GP15_watermassanalysis/blob/988fb74f44b7c62b32c6a11752a9bcbb06471783/endmember_df_intermediateanddeep.csv
# for an example
csv_file = "endmember_df_intermediateanddeep.csv" #this is a path to a .csv
#                                                  file containing the
#                                                  end-members.
endmember_name_column = "watermass_name" #what is the header of the column
#                                         specifying the endmember names.

#If a parameter is specified with params.xxxx, then xxxx must correspond to a
# column header in both the observations and the end-members .csv file
[params.potential_temp] #This specifies that "potential_temp" will be a
#                        parameter used in the OMPA analysis
weight = 56.0
remineralized = false

[params.practical_salinity]
weight = 80.0
remineralized = false

[params.silicate]
weight = 3.0
remineralized = false

[params.nitrate]
weight = 5.0
remineralized = true #remineralized parameters must be specified using
#                     remineralized=true
#Remineralized parameters also need the remineralization ratios to be specified.
#TODO: figure out how to explain these ratios better...these are the ratios
# relative to oxygen, but with a sign-flip, and also the rows are 'coupled',
# meaning that you put the rows for nitrate/phosphate/oxygen one below the
# other and then look at each column to get the different redfield ratios
# specified by each column. That is why each remineralized parameter needs
# to have the same number of entries under 'ratios', even if the entries
# are identical. Also, the effective remineralization ratio that is used
# in the solution for each observation will be somwhere within the convex hull
# of the redfield ratios specified in this configuration file.
ratios = [0.10330578512, 0.10330578512]

[params.phosphate]
weight = 5.0
remineralized = true
ratios = [0.01036168132, 0.00327210989]

[params.oxygen]
weight = 1.0
remineralized = true
ratios = [-1, -1] #sign is negative because oxygen is consumed when
#                  nitrate/phosphate are produced

#If a penalty is specified with endmember_penalties.endmembername.fieldname,
# then fieldname must be a column header in the observations file and
# endmembername must be the name of an endmember in the endmembers file.
#Note: multiple types of penalties (on latitude, on sigma0, etc) can be
# specified for an endmember; each one would just need to given a section
# like endmember_penalties.endmembername.fieldname, and the penalties would
# be added up under-the-hood.
[endmember_penalties.LCDW.sigma0]
type="density_default"
lowerbound = 27.72

[endmember_penalties.AABW.sigma0]
type="density_default"
lowerbound = 27.72

[endmember_penalties.PSUW.lat]
type="latlon_default"
lowerbound = 10
#'upperbound' can also be specified here if desired.

#type="latlon_default" and type="density_default" specify default values for
# the shape of the exponential function used to impose the penalty. These
# are determined by 'alpha', which controls the rate of exponential increase,
# and 'beta' which is a linear scaling of the penalty. The default values
# of alpha and beta can be overridden. You can also set type="other", which
# specifies no defaults for alpha and beta and forces the user to
# provide both.
[endmember_penalties.PSUW.depth]
type="other"
upperbound = 3000
alpha=0.001 #controls rate of exponential increase
beta=50 #linearly scales penalty

#Settings specifiying how to export the results
[export]
csv_output_name="ompa_soln.csv"
orig_cols_to_include = ["lat", "lon", "depth", "stnnbr", "geotrc_ID"] #columns
#                              from the obs_df you want repeated in the output
export_orig_param_vals=true #whether to repeat the original parameter values.
export_residuals=true #residuals are observed-reconstructed. Columns will have
#                      the title "[parametername]_resid"
export_endmember_fracs=true #the endmember fractional compositions. Columns
#                            will have the title "[endmembername]_frac"
export_oxygen_deficit=true #the oxygen deficit used when doing the OMPA fit.
#                           Column will have the title "total_oxygen_deficit"
export_conversion_ratios=true #the effective remineralization ratio used when
#                              doing the OMPA fit. Columns will have the
#                              name "oxygen_to_[paramname]_ratio".
export_endmember_usage_penalties=true #whether to include the penalty function
#                                      that ended up being used. Columns will
#                                      have the title "[endmembername]_penalty"
```
As of version v0.3.0.3-alpha, you can supply multiple configuration files to pyompa, and the contents of the configuration files will be concatenated together under-the-hood and treated as though they were one single configuration file. This can be useful when you want to try out different settings for only one part of the configuration (e.g. if you want to explore different parameter weightings, you can separate out the parameters configuration into its own file).

## Running from a colab notebook/in python
See https://github.com/nitrogenlab/GP15_watermassanalysis/blob/main/config_files_demo.ipynb for an example.
