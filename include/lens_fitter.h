#ifndef LENS_FITTER
#define LENS_FITTER

#include <lensing_classes.h>

// Stores the Einasto tables we interpolate over
extern einTable einKappa    ;
extern einTable einKappaAvg ;

void rollBall(        densProfile   &ball ,  // Ball to roll
                      double        &chi2 ,  // Chi2 value
               const  double        *gArr ,  // g values we are fitting
               const  double        *dArr ,  // distance array corresponding to g values
               const  double     *gErrArr ,  // Errors associated with g
               const  double       sigmaC ,  // Critical surface density to use, of sources
               const  userInfo          u );


void rollBall1D(      densProfile   &ball ,  // Ball to roll
                      double        &chi2 ,  // Chi2 value
               const  double        *gArr ,  // g values we are fitting
               const  double        *dArr ,  // distance array corresponding to g values
               const  double     *gErrArr ,  // Errors associated with g
               const  double       sigmaC ,  // Critical surface density to use, of sources
               const  userInfo          u );


void rollingFitDensProfile(
                          densProfile   &  profile ,  // Density profile we are outputting
                    const userInfo               u ,  // Info from the user
                    const double       *      gArr ,  // RTS binned array we "observed"
                    const double       *      dArr ,  // Distance binned array
                    const double       *   gErrArr ); // Error array in RTS


double * generateNFWTruncRTS(
                    const densProfile     &lens ,  // Input density profile to generate profile for
                    const double         N_bins ,  // Actual information from the halo
                    const double          *dist ,  // Projected distances between source and lens
                    const double           SigC ); // Critical surface density of sources


void generateNFWTruncRTS(
                          double          *gArr ,  // RTS array to output
                    const densProfile     &lens ,  // Input density profile to generate profile for
                    const double         N_bins ,  // Actual information from the halo
                    const double          *dist ,  // Projected distances between source and lens
                    const double           SigC ); // Critical surface density of sources


void generateNFWRTS(
                          double          *gArr ,  // RTS array to output
                    const densProfile     &lens ,  // Input density profile to generate profile for
                    const double         N_bins ,  // Actual information from the halo
                    const double          *dist ,  // Projected distances between source and lens
                    const double           SigC ); // Critical surface density of sources

double *generateNFWRTS(
                    const densProfile     &lens ,  // Input density profile to generate profile for
                    const int            N_bins ,  // Actual information from the halo
                    const double          *dist ,  // Projected distances between source and lens
                    const double           SigC ); // Critical surface density of sources


double    SDNFW( const      double               r ,  // Input radius to calc SD at
                 const densProfile      inpProfile ,  // Input NFW profile
                            int              db = 0);

double    SDAvgNFW( const double               r ,  //Input radius to calc SD at
                    const densProfile inpProfile ,  //Input NFW profile
                          int             db = 0 );

double    SDNFWFull( const double     r ,  //Distance to evaluate SD of NFW profile at
                     const double   r_s ,  //Scale radius of profile
                     const double rho_o ); //Initial density of profile

double SDAvgNFWFull( const double     r ,  //Distance to evaluate SD of NFW profile at
                     const double   r_s ,  //Scale radius of profile
                     const double rho_o ); //Initial density of profile


void generateEinRTS(
		          double           *gArr,  //Radially averaged RTS array function will return
		    const densProfile      &lens,  //Input density profile
            const userInfo             u,  //Info from user
		    const double     *sourceDist,  //Projected radial distance of sources to lens centers
		    const double        sourceSc); //Critical surface density of a source

double *generateEinRTS(
		    const densProfile      &lens,  //Input density profile
            const userInfo             u,  //Info from user
		    const double     *sourceDist,  //Projected radial distance of sources to lens centers
		    const double        sourceSc); //Critical surface density of a source
double *generateEinRTS(
		    const densProfile      &lens,  //Input density profile
                    const userInfo             u,  //Info from user
                    const int                  N,
		    const double     *sourceDist,  //Projected radial distance of sources to lens centers
		    const double        sourceSc); //Critical surface density of a source


double interpolateEinRTS(  double        x ,  // r/r_s
                           double        a ,  // alpha
                           einTable  table ); // Table to interpolate on

void fitDensProfile(
                          densProfile   &  profile ,  // Density profile we are outputting
                    const haloInfo      &     halo ,  // Info about parent halo
                    const userInfo               u ,  // Info from the user
                    const double       *      gArr ,  // RTS binned array we "observed"
                    const double       *      dArr ,  // Distance binned array
                    const double       *   gErrArr );


// Concentration estimates from Klypin 2014
double klypinC( double M ) ;

double cosmoRvir( double M, double z ) ; // In Mpch

#endif
