#ifndef OUTPUT_FUNCTIONS
#define OUTPUT_FUNCTIONS

//#include <gridmap.h>
//#include <image_processing.h>

#include "lensing_classes.h"


/*
int  writeAngRTS( haloInfo   & h ,
                  userInfo     u ,
                  PixelMap  gTan ,
                  PixelMap  gSec ,
                  PixelMap  dMap );
//*/

void generateParamfile( std::string haloName );

std::string getHaloFile( int index );


void writeProfileFits( char       fileName[500] ,   // Filename of output file
                       userInfo               u ,   // User input
                       haloInfo               h ,   // Info on our halo
                       densProfile      tot_ein ,   // Einasto   density profile
                       densProfile      tot_nfw ,   // NFW Full  density profile
                       densProfile      tot_nfT ,   // NFW trunc density profile
                       double      * tot_einErr ,   // Einasto   errors
                       double      * tot_nfwErr ,   // NFW Full  errors
                       double      * tot_nfTErr ,   // NFW trunc errors
                       densProfile      tan_ein ,   // Einasto   density profile
                       densProfile      tan_nfw ,   // NFW Full  density profile
                       densProfile      tan_nfT ,   // NFW trunc density profile
                       double      * tan_einErr ,   // Einasto   errors
                       double      * tan_nfwErr ,   // NFW Full  errors
                       double      * tan_nfTErr ,   // NFW trunc errors
                       double      *      dList ,
                       double      *    totList ,
                       double      * totStdList ,
                       double      *    tanList ,
                       double      * tanStdList );



#endif // OUTPUT_FUNCTIONS
