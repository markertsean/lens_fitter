#ifndef PIXELMAP_FUNCTIONS
#define PIXELMAP_FUNCTIONS

#include <gridmap.h>
#include <image_processing.h>
#include "lensing_classes.h"
#include "my_utilities.h"



// Box-Muller transformation to provide gaussian distribution
double gaussErr( double   sigma ,
                 int       Ngal ); // Number of galaxies

// Calculates jackknife errors for each profile set
void jacknife( densProfile  *profile ,
               int         N_samples ,
               double       * errArr );


#endif
