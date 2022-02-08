from __future__ import division, print_function


class BasicThermoclineSplitResult(object):

    def __init__(self, thermocline_obs_df, belowthermocline_obs_df):
        self.thermocline_obs_df = thermocline_obs_df
        self.belowthermocline_obs_df = belowthermocline_obs_df

    def viz_perstation_thermoclineboundary(self,
            station_number, xaxis_colname, yaxis_colname,
            colorbar_colname, flip_y=True):
        
                   
class AbstractThermoclineSplitter(object): #for defining the API

    #obs_df = observations data frame
    def __call__(self, obs_df): 
        #Should return thermocline_obs_df and belowthermocline_obs_df
        raise NotImplementedError() 


class FixedThresholdThermoclineSplitter(AbstractThermoclineSplitter):

    def __init__(self, splitby_colname, therm_begin_val, therm_end_val,
                       increasing_with_depth):

        #station_colname is optional
        if (increasing_with_depth):
            assert therm_begin_val <= therm_end_val,\
             ("increasing_with_depth is set to true, but therm_begin_val ("
              +str(therm_begin_val)+") is > therm_end_val ("
              +str(therm_end_val)+")")
        else:
            assert therm_begin_val >= therm_end_val,\
             ("increasing_with_depth is set to false, but therm_begin_val ("
              +str(therm_begin_val)+") is < therm_end_val ("
              +str(therm_end_val)+")")

        self.splitby_colname = splitby_colname
        self.therm_begin_val = therm_begin_val
        self.therm_end_val = therm_end_val
        self.increasing_with_depth = increasing_with_depth

    def __call__(self, obs_df):
        if (self.increasing_with_depth):
            thermocline_obs_df = obs_df[
                 (obs_df[splitby_colname] >= self.therm_begin_val) &&
                 (obs_df[splitby_colname] <= self.therm_end_val)]
            belowthermocline_obs_df = obs_df[
                 (obs_df[splitby_colname] > self.therm_end_val)]
        else:
            thermocline_obs_df = obs_df[
                 (obs_df[splitby_colname] <= self.therm_begin_val) &&
                 (obs_df[splitby_colname] >= self.therm_end_val)]
            belowthermocline_obs_df = obs_df[
                 (obs_df[splitby_colname] < self.therm_end_val)]

        return BasicThermoclineSplitResult(
                thermocline_obs_df=thermocline_obs_df,
                belowthermocline_obs_df=belowthermocline_obs_df)


class PerStationThermoclineSplitter(object):

    #station_colname is used for grouping the observations
    def __init__(self, station_colname):
        self.station_colname = station_colname



