from __future__ import division, print_function
import cvxpy as cp
import numpy as np


class OMPAProblem(object):
    """
    Core class for conducting the OMPA with cvxpy
    """

    def __init__(self, watermass_df, obs_df, paramsandweighting,
                       smoothness_lambda):
        self.watermass_df = watermass_df
        self.obs_df = obs_df
        self.paramsandweighting = paramsandweighting
        #split up the paramsandweighting into params and weighting
        self.params_to_use, self.weighting = [list(x) for x in 
                                              zip(*paramsandweighting)]
        self.smoothness_lambda = smoothness_lambda

    def core_solve(self, A, b, pairs_matrix, verbose=True):
        #We are going to solve the following problem:
        # Minimize (x A - b)^2
        # Subject to x >= 0, cp.sum(x, axis=1) == 1
        # x has dimensions of observations X end_members
        # A has dimensions of end_members X parameteres
        # b has dimensions of observations X parameters
       
        x = cp.Variable(shape=(len(b), len(A)))
        obj = cp.sum_squares(x@A - b)
        if (self.smoothness_lambda is not None):
            obj += self.smoothness_lambda*cp.sum_squares(pairs_matrix@x)
        obj = cp.Minimize(obj)
        
        constraints = [x >= 0, cp.sum(x,axis=1)==1]
        prob = cp.Problem(obj, constraints)
        prob.solve(verbose=False, max_iter=50000)
        #settign verbose=True will generate more print statements and slow down the analysis
        
        print("status:", prob.status)
        print("optimal value", prob.value)

        if (prob.status=="infeasible"):
            water_mass_fractions = None
            residuals_squared = None
        else:
          water_mass_fractions = x.value
          residuals_squared = np.sum(
              np.square((water_mass_fractions@A) - b), axis=-1)
        
        return water_mass_fractions, residuals_squared, prob

    def solve(self):

        watermass_df = self.watermass_df
        obs_df = pd.DataFrame(self.obs_df)
        params_to_use = self.params_to_use
        weighting = np.array(self.weighting)
            
        A = np.array(watermass_df[params_to_use])
        b = np.array(obs_df[params_to_use])
        
        print("params to use:", params_to_use)
        print("param weighting:", weighting)
        A = A*weighting[None,:]
        b = b*weighting[None,:]

        if (self.smoothness_lambda is not None):
            pairs_matrix = make_pairs_matrix(
              obs_df=obs_df,
              depth_metric="depth",
              depth_scale=1.0,
              nneighb=4)
        else:
            pairs_matrix = None

        water_mass_fractions, residuals_squared, prob = self.core_solve(
            A=A, b=b, pairs_matrix=pairs_matrix)      
        self.prob = prob    
        
        if (water_mass_fractions is not None):
            print("objective:", np.sum(residuals_squared))
            param_reconstruction = (water_mass_fractions@A)/weighting[None,:]
            param_residuals = b/weighting[None,:] - param_reconstruction

            self.water_mass_fractions = water_mass_fractions
            self.param_reconstruction = param_reconstruction
            self.param_residuals = param_residuals


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


