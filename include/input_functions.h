#ifndef INPUT_FUNCTIONS
#define INPUT_FUNCTIONS

#include <cosmo.h>
#include <iostream>


void readInpFile   (  	      userInfo        &inpInfo ,  // Object we write to, contains parameters governing options
                    const std::string     userFileName ); // Name of the file to read


einTable readFoxH  (          userInfo        &     u  ,
                    const     int             fileType );

#endif
