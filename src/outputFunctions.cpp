#include "output_functions.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/stat.h>
#include <vector>
#include "lens_fitter.h"


void checkDir( std::string dirName ){
    struct stat sb;
    if ( stat( dirName.c_str(), &sb ) !=0 ){
      char str[500];
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

  char     fileName[500];
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
                       double      * tanStdList )
{

  checkDir( u.getOutputPath() );


  FILE *pFile;

  pFile = fopen( fileName, "w" );


  fprintf( pFile , "M           %14.6e\n", h.getM    () );
  fprintf( pFile , "b/a         %10.6f\n", h.getBA   () );
  fprintf( pFile , "gamma       %10.6f\n", h.getGamma() );


  fprintf( pFile , "Tot_NFW_Full %10.6f %10.6f %10.6f %10.6f %10.3e %10.3e %10.3e %10.3e\n" ,
                    log10( tot_nfw.getM_enc() ), tot_nfw.getC(), tot_nfw.getR_max(),               -1.0, tot_nfwErr[1], tot_nfwErr[0], tot_nfwErr[3],          -1.0);
  fprintf( pFile , "Tot_NFW_Trnc %10.6f %10.6f %10.6f %10.6f %10.3e %10.3e %10.3e %10.3e\n" ,
                    log10( tot_nfT.getM_enc() ), tot_nfT.getC(), tot_nfT.getR_max(),               -1.0, tot_nfTErr[1], tot_nfTErr[0], tot_nfTErr[3],          -1.0);
  fprintf( pFile , "Tot_Ein      %10.6f %10.6f %10.6f %10.6f %10.3e %10.3e %10.3e %10.3e\n" ,
                    log10( tot_ein.getM_enc() ), tot_ein.getC(), tot_ein.getR_max(), tot_ein.getAlpha(), tot_einErr[1], tot_einErr[0], tot_einErr[3], tot_einErr[2]);

  fprintf( pFile , "Tan_NFW_Full %10.6f %10.6f %10.6f %10.6f %10.3e %10.3e %10.3e %10.3e\n" ,
                    log10( tan_nfw.getM_enc() ), tan_nfw.getC(), tan_nfw.getR_max(),               -1.0, tan_nfwErr[1], tan_nfwErr[0], tan_nfwErr[3],          -1.0);
  fprintf( pFile , "Tan_NFW_Trnc %10.6f %10.6f %10.6f %10.6f %10.3e %10.3e %10.3e %10.3e\n" ,
                    log10( tan_nfT.getM_enc() ), tan_nfT.getC(), tan_nfT.getR_max(),               -1.0, tan_nfTErr[1], tan_nfTErr[0], tan_nfTErr[3],          -1.0);
  fprintf( pFile , "Tan_Ein      %10.6f %10.6f %10.6f %10.6f %10.3e %10.3e %10.3e %10.3e\n" ,
                    log10( tan_ein.getM_enc() ), tan_ein.getC(), tan_ein.getR_max(), tan_ein.getAlpha(), tan_einErr[1], tan_einErr[0], tan_einErr[3], tan_einErr[2]);

// d tot std NFWt NFWf Ein tan std NFWt NFWf Ein
    double * totEinVals = generateEinRTS( tot_ein, u, dList, u.getSigmaCrit() );
    double * tanEinVals = generateEinRTS( tan_ein, u, dList, u.getSigmaCrit() );

    double * totNFWVals =      generateNFWRTS( tot_nfw, u.getNbins(), dList, u.getSigmaCrit() );
    double * tanNFWVals =      generateNFWRTS( tan_nfw, u.getNbins(), dList, u.getSigmaCrit() );

    double * totNFTVals = generateNFWTruncRTS( tot_nfT, u.getNbins(), dList, u.getSigmaCrit() );
    double * tanNFTVals = generateNFWTruncRTS( tan_nfT, u.getNbins(), dList, u.getSigmaCrit() );


    for ( int i = 0 ; i < u.getNbins(); ++i )
    {
        fprintf(pFile, "%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n", dList[i] ,
                                           totList[i], totStdList[i], totNFTVals[i], totNFWVals[i], totEinVals[i] ,
                                           tanList[i], tanStdList[i], tanNFTVals[i], tanNFWVals[i], tanEinVals[i] );
    }


  fclose( pFile );

    delete [] totNFTVals ;
    delete [] tanNFTVals ;
    delete [] totNFWVals ;
    delete [] tanNFWVals ;
    delete [] totEinVals ;
    delete [] tanEinVals ;

  std::cout << "Wrote file: " << fileName << std::endl ;

}



