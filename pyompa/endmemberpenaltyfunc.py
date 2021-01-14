from __future__ import division, print_function
import numpy as np
from .util import assert_in, assert_compatible_keys 


def get_exponential_from_bounds_func(alpha, beta,
                                     lowerbound=-np.inf,
                                     upperbound=np.inf):
    assert upperbound > lowerbound, (upperbound, lowerbound)
    return lambda x: beta*(np.exp(
        alpha*np.maximum(0, np.maximum(lowerbound-x, x-upperbound)))-1)


#same as get_exponential_from_bounds_func but with defaults for alpha and beta
def get_default_density_exp_penalty_func(alpha=0.5, beta=1000, **kwargs):
    return get_exponential_from_bounds_func(
            alpha=alpha, beta=beta, **kwargs) 


#same as get_exponential_from_bounds_func but with defaults for alpha and beta
def get_default_latlon_exp_penalty_func(alpha=0.05, beta=100, **kwargs):
    return get_exponential_from_bounds_func(
            alpha=alpha, beta=beta, **kwargs) 


def get_combined_penalty_func(colname_to_penaltyfunc):
    def penalty_func(df):
        total_penalty = None
        for colname, penaltyfunc in colname_to_penaltyfunc.items():
            values = np.array(df[colname])
            penalty = penaltyfunc(values) 
            if (total_penalty is None):
                total_penalty = penalty
            else:
                total_penalty += penalty
        return total_penalty
    return penalty_func


class EndMemExpPenaltyFunc(object):
    
    #mapping from spectype to factory functions that manufacture the
    # penalty functions
    SPECTYPE_TO_FACTORYFUNC = {
        'density_default': get_default_density_exp_penalty_func,
        'latlon_default': get_default_latlon_exp_penalty_func,
        'other': get_exponential_from_bounds_func}

    def __init__(self, spec):
        self.spec = spec 
        self.validate_spec()
        self.process_spec()

    def validate_spec(self):
        for colname,spec_for_col in self.spec.items():
            assert_compatible_keys(
              the_dict=spec_for_col,
              allowed=['type', 'lowerbound', 'upperbound', 'alpha', 'beta'],
              errorprefix="Problem with "+str(colname)+" penalty config: ") 
            assert_has_keys(the_dict=spec_for_col, required_keys=["type"],
                            errorprefix="Problem with "+str(colname)
                                        +" penalty config: ")
            assert_in(
                value=spec_for_col['type'],
                allowed=list(self.SPECTYPE_TO_FACTORYFUNC.keys()),
                errorprefix="Problem with type for "
                            +str(colname)+" penalty config; ")

    def process_spec(self):
        colname_to_penaltyfunc = {}
        for colname, spec_for_col in self.spec.items(): 
            factory_func = self.SPECTYPE_TO_FACTORYFUNC[spec_for_col['type']] 
            #gather the arguments going into the factory function
            factory_func_keys = spec_for_col.copy() 
            del factory_func_keys['type']
            penalty_func = factory_func(**factory_func_keys)
            colname_to_penaltyfunc[colname] = penalty_func
        self.penalty_func = get_combined_penalty_func(colname_to_penaltyfunc)

    def __call__(self, *args, **kwargs):
        return self.penalty_func(*args, **kwargs) 
