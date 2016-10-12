#include "output_functions.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/stat.h>
#include <vector>



void checkDir( std::string dirName ){
    struct stat sb;
    if ( stat( dirName.c_str(), &sb ) !=0 ){
      char str[100];
      sprintf( str, "mkdir %s", dirName.c_str());
      system( str );
      logMessage( std::string("Wrote directory: ") + dirName );
    }
}

// Returns true if file exists
bool checkFile( char dirName[] ){
    struct stat sb;
    if ( stat( dirName, &sb ) !=0 ){
      return false;
    }
      return true;
}


bool checkOutputExists( userInfo       u , // If output files exist before first run, abort
                        haloInfo       h ){

  char     fileName[100];
  sprintf( fileName, "%sHalo_%010li_densFits.dat", u.getOutputPath().c_str(), h.getID() );

  return checkFile( fileName );

}



// Reads the halo fits file from a file, to generate file for GLAMER
// If done with file, returns ""
std::string getHaloFile ( int index ){

  std::string   inputLine;
  std::ifstream halo_file("haloList.dat" );

  int counter = 0;

  if ( halo_file.is_open() ){
    while( getline( halo_file, inputLine ) ){
      if ( counter == index ) break;
      if ( halo_file.eof()  ) break;
      ++counter;
    }
  } else{
    std::cout << "Code requires files \"haloList.dat\" in executing directory" << std::endl;
    logMessage("Couldn't open haloList");
    exit(1);
  }

  halo_file.close();

  if ( counter != index ) inputLine = "";

  return inputLine;
}




/*
void writeProfileFits( userInfo        u ,   // User input
                       haloInfo        h ,   // Info on our halo
                       densProfile   ein ,   // Einasto   density profile
                       densProfile   nfw ,   // NFW Full  density profile
                       densProfile   nfT ,   // NFW trunc density profile
                       double    *einErr ,   // Einasto   errors
                       double    *nfwErr ,   // NFW Full  errors
                       double    *nfTErr ,   // NFW trunc errors
                       int       haloNum ){  // How many times we've written, first time we need to write halo info


  checkDir( u.getOutputPath() );

  double integ = u.getIntegLength(); // Need for file name

  if ( integ == -1 ){                // If sphere, just mark as 0 integ length
    integ = 0.0;
  }

  char     fileName[100];
  sprintf( fileName, "%sHalo_%010li_densFits.dat", u.getOutputPath().c_str(), h.getID() );



  if ( haloNum == 1 ){ // First time through, write the general halo info
    FILE *pFile;

    pFile = fopen( fileName, "w" );

    fprintf( pFile , "ID          %10li\n" , h.getID   () );
    fprintf( pFile , "M           %14.6e\n", h.getM    () );
    fprintf( pFile , "C           %10.6f\n", h.getC    () );
    fprintf( pFile , "R_max       %10.6f\n", h.getRmax () );
    fprintf( pFile , "Z           %10.6f\n", h.getZ    () );
    fprintf( pFile , "b/a         %10.6f\n", h.getBA   () );
    fprintf( pFile , "c/a         %10.6f\n", h.getCA   () );
    fprintf( pFile , "phi         %10.6f\n", h.getPhi  () );
    fprintf( pFile , "theta       %10.6f\n", h.getTheta() );
    fprintf( pFile , "alpha       %10.6f\n", h.getAlpha() );
    fprintf( pFile , "gamma       %10.6f\n", h.getGamma() );

    fprintf( pFile , "FOV         %10.6f\n", u.getPhysFOV() );
    fprintf( pFile , "N_pixH      %10i\n"  , u.getNpixH  () );
    fprintf( pFile , "N_pixV      %10i\n"  , u.getNpixV  () );

    fprintf( pFile , "N_src       %10i\n"  , u.getNsrc      () );
    fprintf( pFile , "Z_src       %10.6f\n", u.getSourceZ   () );
    fprintf( pFile , "sigma_shape %10.6f\n", u.getShapeNoise() );

    fclose(  pFile );
  }

  FILE *pFile;

  pFile = fopen( fileName, "a+" );



  fprintf( pFile , "IntegLength %10.6f\n"        , u.getIntegLength() );
  fprintf( pFile , "ImageMass   %14.6e\n"        , u.getImageMass  () );
  fprintf( pFile , "NFW_Full %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n" , log10( nfw.getM_enc() ), nfw.getC(),           -1.0, nfwErr[1], nfwErr[0],      -1.0);
  fprintf( pFile , "NFW_Trnc %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n" , log10( nfT.getM_enc() ), nfT.getC(),           -1.0, nfTErr[1], nfTErr[0],      -1.0);
  fprintf( pFile , "Ein      %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n" , log10( ein.getM_enc() ), ein.getC(), ein.getAlpha(), einErr[1], einErr[0], einErr[2]);


  fclose( pFile );

  std::cout << "Appended file: " << fileName << std::endl << std::endl;

}
//*/


