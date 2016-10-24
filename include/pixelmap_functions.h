#ifndef PIXELMAP_FUNCTIONS
#define PIXELMAP_FUNCTIONS

#include <gridmap.h>
#include <image_processing.h>
#include "lensing_classes.h"
#include "my_utilities.h"


// Will take and fill outArr with weighted average of a collapsed mass,
//   setting as our 1D arrays to pass to fitter
int *   avgMArr  (  userInfo        u ,  // User info
                    double   * inpArr ,  // JIMBGR array to collapse into R
                    int      * inpN   ,  // Array containing number of sources to weight by
                    int             i ,  // Integration index
                    int             b ,  // b/a ratio   index
                    int             g ,  // gamma       index
                    int     omitIndex ,  // Index to omit for jack knifing
                    double  ** outArr ); // Array to allocate and populate with average

// Add gaussian error to arrays to fit
void  addgaussUncertaintyArr( double       *inpArr ,
                              double         sigma ,  // Amplitude of shape noise
                              int          * Nsrc  ,  // Array containing number of sources in each bin
                              int       N_elements ); // Number of elements in the array


// Generates gaussian error array to include in error arrays
double * gaussUncertaintyArr( double         sigma ,  // Amplitude of shape noise
                              int          * Nsrc  ,  // Array containing number of sources in each bin
                              int       N_elements ); // Number of elements in the array


// Box-Muller transformation to provide gaussian distribution
double gaussErr( double   sigma ,
                 int       Ngal ); // Number of galaxies


// Collapses average array to match that of a collapsed M
haloInfo  avgMHaloInfo( userInfo    u ,
                        haloInfo  * h ,
                        int       * n ,
                        int         i ,
                        int         b ,
                        int         g );



// Calculates jackknife errors for each profile set
void jacknife( densProfile  *profile ,
               int         N_samples ,
               double       * errArr );


#endif
