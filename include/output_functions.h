#ifndef OUTPUT_FUNCTIONS
#define OUTPUT_FUNCTIONS

#include <gridmap.h>
#include <image_processing.h>

#include "lensing_classes.h"



int  writeAngRTS( haloInfo   & h ,
                  userInfo     u ,
                  PixelMap  gTan ,
                  PixelMap  gSec ,
                  PixelMap  dMap );


void generateParamfile( std::string haloName );

std::string getHaloFile( int index );


void writeProfileFits( char       fileName[100] ,   // Filename of output file
                       userInfo               u ,   // User input
                       haloInfo               h ,   // Info on our halo
                       densProfile          ein ,   // Einasto   density profile
                       densProfile          nfw ,   // NFW Full  density profile
                       densProfile          nfT ,   // NFW trunc density profile
                       double           *einErr ,   // Einasto   errors
                       double           *nfwErr ,   // NFW Full  errors
                       double           *nfTErr ,   // NFW trunc errors
                       int                  N_h );  // Number of halos written to file


// Writes a short, easy read in file
void writeShort( userInfo      u ,
                 double   * gTot ,
                 double   * gTan ,
                 double   * d    ,
                 int      * n_s  ,
                 int      * n_h  ,
                 haloInfo * bh   );



#endif // OUTPUT_FUNCTIONS
