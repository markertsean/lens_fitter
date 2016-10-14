#ifndef INPUT_FUNCTIONS
#define INPUT_FUNCTIONS

#include <cosmo.h>
#include <iostream>


void readInpFile   (  	      userInfo        &inpInfo ,  // Object we write to, contains parameters governing options
                    const std::string     userFileName ); // Name of the file to read

// Wrapper function for reading source files, will locate files and invoke reader function
int  readSources(  userInfo    u    ,  // User input
                   double    * d    ,  // Array of distances
                   double    * gTot ,  // Array of gTot
                   double    * gTan ,  // Array of gTan
                   int       * N    ,  // Array counting number in each bin
                   int       * N_h  ); // Array containing halo count in each bin


einTable readFoxH  (          userInfo        &     u  ,
                    const     int             fileType );

#endif
