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


// Writes a short, easy read in file
void writeShort( userInfo      u ,
                 double   * gTot ,
                 double   * gTan ,
                 double   * d    ,
                 int      * n_s  ,
                 int      * n_h  ,
                 haloInfo * bh   )
{
    std::string storedData = "storedData.dat";

    FILE *pFile ;
    pFile = fopen( storedData.c_str(), "w" );


    // Prints binned halo info, averages for the bins
    for ( int i = 0; i < u.getN_IBin(); ++i ){
    for ( int m = 0; m < u.getN_MBin(); ++m ){
    for ( int b = 0; b < u.getN_BBin(); ++b ){
    for ( int g = 0; g < u.getN_GBin(); ++g ){

        int  H_bin = u.getN_haloBin( i, m, b, g    );
        fprintf( pFile, "%14.6e %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", bh[H_bin].getM     () ,
                                                                                    bh[H_bin].getC     () ,
                                                                                    bh[H_bin].getBA    () ,
                                                                                    bh[H_bin].getCA    () ,
                                                                                    bh[H_bin].getRmax  () ,
                                                                                    bh[H_bin].getPhi   () ,
                                                                                    bh[H_bin].getTheta () ,
                                                                                    bh[H_bin].getGamma () ,
                                                                                    bh[H_bin].getAlpha () );
    }
    }
    }
    }


    // Prints N_halo info
    for ( int i = 0; i < u.getN_IBin(); ++i ){
    for ( int m = 0; m < u.getN_MBin(); ++m ){
    for ( int b = 0; b < u.getN_BBin(); ++b ){

    for ( int g = 0; g < u.getN_GBin(); ++g ){

        int  h = u.getN_haloBin( i, m, b, g    );
        fprintf( pFile, "%10i ", n_h[ h ] );
    }
        fprintf( pFile, "\n" );

    }
    }
    }

    // Prints tables
    for ( int j = 0; j < u.getN_JBin(); ++j ){
    for ( int i = 0; i < u.getN_IBin(); ++i ){
    for ( int m = 0; m < u.getN_MBin(); ++m ){
    for ( int b = 0; b < u.getN_BBin(); ++b ){
    for ( int g = 0; g < u.getN_GBin(); ++g ){

        int k = u.getSrcBin( j, i, m, b, g, 0 );

        // Outputs each row as group
        for ( int r = 0; r < u.getNbins(); ++r ){  fprintf(pFile, "%14i "  ,  n_s[ k + r ] ); } fprintf( pFile, "\n" );
        for ( int r = 0; r < u.getNbins(); ++r ){  fprintf(pFile, "%14.6e ",    d[ k + r ] ); } fprintf( pFile, "\n" );
        for ( int r = 0; r < u.getNbins(); ++r ){  fprintf(pFile, "%14.6e ", gTan[ k + r ] ); } fprintf( pFile, "\n" );
        for ( int r = 0; r < u.getNbins(); ++r ){  fprintf(pFile, "%14.6e ", gTot[ k + r ] ); } fprintf( pFile, "\n" );


    }
    }
    }
    }
    }

    fclose( pFile );

}



void writeProfileFits( char       fileName[100] ,   // Filename of output file
                       userInfo               u ,   // User input
                       haloInfo               h ,   // Info on our halo
                       densProfile          ein ,   // Einasto   density profile
                       densProfile          nfw ,   // NFW Full  density profile
                       densProfile          nfT ,   // NFW trunc density profile
                       double           *einErr ,   // Einasto   errors
                       double           *nfwErr ,   // NFW Full  errors
                       double           *nfTErr ,   // NFW trunc errors
                       int                  N_h ){  // Number of halos written to file


  checkDir( u.getOutputPath() );


  FILE *pFile;

  pFile = fopen( fileName, "a+" );


  fprintf( pFile , "M           %14.6e\n", h.getM    () );
  fprintf( pFile , "C           %10.6f\n", h.getC    () );
  fprintf( pFile , "R_max       %10.6f\n", h.getRmax () );
  fprintf( pFile , "b/a         %10.6f\n", h.getBA   () );
  fprintf( pFile , "gamma       %10.6f\n", h.getGamma() );
  fprintf( pFile , "N_lens      %10i\n"  , N_h          );


  fprintf( pFile , "NFW_Full %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n" ,
                                                                log10( nfw.getM_enc() ), nfw.getC(), nfw.getR_max(),           -1.0, nfwErr[1], nfwErr[0], nfwErr[3],      -1.0);
  fprintf( pFile , "NFW_Trnc %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n" ,
                                                                log10( nfT.getM_enc() ), nfT.getC(), nfT.getR_max(),           -1.0, nfTErr[1], nfTErr[0], nfTErr[3],      -1.0);
  fprintf( pFile , "Ein      %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n" ,
                                                                log10( ein.getM_enc() ), ein.getC(), ein.getR_max(), ein.getAlpha(), einErr[1], einErr[0], einErr[3], einErr[2]);


  fclose( pFile );

  std::cout << "Appended file: " << fileName << std::endl ;

}



