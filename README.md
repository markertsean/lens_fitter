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





√
1/ N , it is thus beneficial to optimize the radial bin sizes to hold as many
50sources as possible. Shape noise is Gaussian, with 0 mean and variance
σ s 2 =
σ e 2
N bin
(46)
where N bin is the number of sources in a bin and σ e is the intrinsic shape noise of
the background galaxies, σ e of 0.3 is typical for most ground based observations
(e.g.Okabe et al. 2010, Hoekstra et al. 2012, Heymans et al. 2012).
We fitted our three profiles without noise added to the bins, artificial shape
noise would have provided a more accurate observational analysis, however,
is not the point of this thesis. This thesis is in part a comparison of the
performance of the density profiles in recovering the halo mass, for which the
no-noise fits provides a more direct comparison.
Another source of noise is from the triaxiality of the lens, and correlated
and uncorrelated LSS along the LOS. Due to the nature of extracting RS mea-
surements from simulations, these are intrinsic sources of noise that cannot be
removed, or decreases with an increased number density of sources. Charac-
terizing these components of noise was aligned with the goals of this thesis,
and is performed in Section 5.
Other components of noise are ignored. Cluster member contamination
(Okabe et al. 2010), photometric redshift errors (Mandelbaum et al. 2008),
unknown source redshifts (Okabe et al. 2010), and intrinsic cluster alignment
(Hopkins et al. 2005) are beyond the scope of this thesis. Many of these are
too intrinsically tied to the observational method themselves, or a physical
process unrelated to lensing. We also ignore the effect of magnification and
size bias of background sources, which effect the shape and number density
51of background sources, though can provide an independent measure of mass
(Umetsu et al. 2011, Coupon et al. 2013).
We estimated the mass of the halos by fitting the radially binned RS profiles
of the halos to the RTS of our three density profiles in the thin lens approx-
imation. We performed fitting using both the RTS and TRS profiles, as a
comparison of direct observational methods, and a test of the true tangential
fit. We perform a two parameter fits, varying the mass, and concentration of
the halos. The Einasto profile had an additional fit parameter, α.
We use a χ 2 fitting metric to test our parameter space
N
χ 2 =
i
g i − g x (r i , M 200 , c)
σ s (r i )
2
(47)
where g i is the RS of the i th bin we are fitting, and g x is the predicted RTS
profile for profile x at radius r i , mass M 200 , and concentration c. Note that
due to our selection of no added noise, we effectively weighted by the square
root of the number of sources in a bin. The values r i and g i were computed
by taking the annular average of all the sources in the bin.
Uncertainties were generated by jack knife resampling of the fits. We con-
structed 100 subsamples by removing an equal number of random sources, so
that every source is removed from only one subsample. We use this to generate
the covariance matrix by reanalyzing each subsample in turn.
A common analysis technique is to stack sources to increase the signal to
noise in a bin, and recover an average mass of a sample. We omit this analysis
from our main fitting results, and focused on fits for each individual halo.
52Table 4: Variables accepted by the lens fitting code
