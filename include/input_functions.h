#ifndef INPUT_FUNCTIONS
#define INPUT_FUNCTIONS

#include <cosmo.h>
#include <iostream>


void readInpFile   (  	      userInfo        &inpInfo ,  // Object we write to, contains parameters governing options
                    const std::string     userFileName ); // Name of the file to read


bool readSources(  haloInfo  & h    ,  // Info on the halo, will be returned for binning
                   userInfo    u    ,  // User input
                   double    * d    ,  // Array of distances
                   double    * gTot ,  // Array of gTot
                   double    * gTan ,  // Array of gTan
                   int       * N    ); // Array counting number in each bin

einTable readFoxH  (          userInfo        &     u  ,
                    const     int             fileType );

#endif
