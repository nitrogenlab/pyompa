# pyompa
This is a python package for OMP analysis. A preprint describing the mathematical formulation is available [here](https://essopenarchive.org/doi/full/10.1002/essoar.10507053.4). If you encounter issues while using this package, please contact Avanti Shrikumar (avanti [at] cs [dot] stanford [dot] edu) or Professor Karen Casciotti (kcasciot [at] stanford [dot] edu).

For code that applied pyompa to perform water mass analysis for the GP15 cruise, see [https://github.com/nitrogenlab/GP15_watermassanalysis](https://github.com/nitrogenlab/GP15_watermassanalysis)

Note that this implementation has key upgrades relative to the MATLAB OMP implementation by M. Tomczak; among them, this package models mass conservation with a hard constraint (i.e. it strives for a residual of 0 in the mass conservation equations), it has support for flexible redfield ratios when modeling remineralization, it has support for soft penalties to encode prior knowledge about end-member distributions, and it has functionality for navigating the ambiguity in underdetermined systems.

Please see the examples folder for links to Colab notebooks that demonstrate how to run pyompa, including how to replicate the functionality of the previous MATLAB OMP implementation if desired.

## Frequently Asked Questions

### How do you handle mass conservation?

One of the improvements that this implementation has over the original OMP formulation is that conservation if mass is implemented as a *hard constraint*, which means that the solver looks for a solution where mass conservation is satisfied with a residual of 0. The solution this produces is analogous to what would happen if (in theory) one were to set conservation of mass to have a weight of infinity in the original MATLAB OMP formulation. There are strong theoretical reasons to prefer a hard constraint for conservation of mass, which will be elaborated on in our upcoming paper; the gist is that if conservation of mass is allowed to be violated, it can be used by the OMP solver to “hide” large residuals in other variables, thereby producing poor-quality solutions. This is why we use the hard constraint, which disallows violation of conservation of mass.

### Do you normalize/standardize your water type matrix?

There are two aspects to the normalization/standardization that is applied in the original MATLAB OMP formulation. The first is to subtract the mean parameter value across the water types, and the other is to divide by the standard deviation of the parameter value across the water types. We can show that if the mass conservation equation is satisfied exactly (as is ensured in our implementation, due to the hard constraint), then the mean normalization does not change the ideal solution. As for rescaling by the standard deviation: mathematically, this is equivalent to changing the user-specified parameter weights such that they get divided by the standard deviation. We have decided not to do this rescaling in this pyompa implementation because it can make it harder to figure out what the ultimate parameter weights end up being, and it also makes the ultimate parameter weights dependent on the water type matrix in a way that the user might not realize. If you want to recapitulate the effect of this rescaling from the original MATLAB OMP implementation, you can manually adjust the parameter weights by dividing them by the standard deviation of the parameter value across the water types.
