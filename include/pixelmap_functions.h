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




// Box-Muller transformation to provide gaussian distribution
double gaussErr( double   sigma ,
                 int       Ngal ); // Number of galaxies

// Calculates jackknife errors for each profile set
void jacknife( densProfile  *profile ,
               int         N_samples ,
               double       * errArr );


#endif
