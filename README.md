///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Warning:
This fitting code is designed for use with a shear calculation code.
The three branches, master, fitR, and individual, contain different 
styles of fitting. fitR fits up to r_vir, individual fits each halo
individually, and master contains the original stacking of multiple
sources. I have not had an opportunity to properly merge the branches
together, I was doing a lot of testing and never set up a stable release. 
Due to this, the code may not behave as expected.

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

The lens_fitter code is designed to take sources output from my shear calculation 
code, and generate output of halo masses estimated from NFW, Einasto, and 
truncated NFW profile density profile fits to the reduced shear.

This code locates files in allFiles.dat, and reads sources and halo information
from the files. There is an option to add noise based on the shapeNoise input variable,
though regardless the shear fit is effectively weighted based on the number of halos in
the bins.

The fitter places random probes in the boundaries of the fit parameters, and steps the
probes in the direction of the parameter space that minimized the Chi^2 fit. This is 
repeated for every density profile, and for each jack knife subset to calculate the
uncertainties.

The following parameters are used in lensUserParams.dat: 


N bins       - Default 20, number of radial bins for profile fitting

N binsJac    - Default 10, number of jack knife subsamples to use to calculate uncertainties

N threads    - Default  1, number of omp threads to use while fitting

N fitAttempt - Default 500, maximum number of steps to take while fitting an individual probe

N probes     - Default 100, number of probes to use while fitting

N consistent - Default 20, number of times must meet convergence criteria before exiting fit

fitTolerance - Default 1e-5, fit differential size to count as converged

shapeNoise   - Default 0.3, shape noise to use for source uncertainties

useNoise     - Default   1, turns on/off the addition of shape noise to shear measurements

fox2012F     - Default src/foxH2012.dat, location of file used to interpolate Einasto surface density 

fox2123F     - Defualt src/foxH2123.dat, location of file used to interpolate Einasto average surface density

inputPath    - Directory containing the input data files

firstFile    - ID number of halos to start search for in inputPath

lastFile     - ID number of halos to end search for in inputPath
