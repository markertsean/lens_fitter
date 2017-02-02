#include <cstring>
#include <cmath>
#include <math.h>
#include <vector>
//#include <CCfits/CCfits>
//#include <slsimlib.h>
#include <stdio.h>
#include "lensing_classes.h"



bool readSourceFile(   FILE      * pFile ,
                       userInfo    u     ,  // User input
                       double    * dist  ,  // Array of distances
                       double    * gTot  ,  // Array of gTot
                       double    * gTan  ,  // Array of gTan
                       int       * N     ,  // Array counting number in each bin
                       haloInfo  * h     )  // Array containing halo info for bins
{

    char   inpC1[35],inpC2[35], inpC3[35];
    double dMax ;

    double M, ba, gamma, c, r_max, ca, phi, theta, alpha, integ;

    int N_src, ID ;

    fscanf(pFile,"%s%s",inpC1,inpC2) ; ID    = atoi( inpC2 );
    fscanf(pFile,"%s%s",inpC1,inpC2) ; M     = atof( inpC2 );
    fscanf(pFile,"%s%s",inpC1,inpC2) ; c     = atof( inpC2 );
    fscanf(pFile,"%s%s",inpC1,inpC2) ; r_max = atof( inpC2 );
    fscanf(pFile,"%s%s",inpC1,inpC2) ; // Z
    fscanf(pFile,"%s%s",inpC1,inpC2) ; ba    = atof( inpC2 );
    fscanf(pFile,"%s%s",inpC1,inpC2) ; ca    = atof( inpC2 );
    fscanf(pFile,"%s%s",inpC1,inpC2) ; phi   = atof( inpC2 );
    fscanf(pFile,"%s%s",inpC1,inpC2) ; theta = atof( inpC2 );
    fscanf(pFile,"%s%s",inpC1,inpC2) ; alpha = atof( inpC2 );
    fscanf(pFile,"%s%s",inpC1,inpC2) ; gamma = atof( inpC2 );
    fscanf(pFile,"%s%s",inpC1,inpC2) ; integ = atof( inpC2 );
    fscanf(pFile,"%s%s",inpC1,inpC2) ; // IntegM
    fscanf(pFile,"%s%s",inpC1,inpC2) ; dMax   = atof( inpC2 ) / 2.0 * std::sqrt(2.); // Go through long axis of image
    fscanf(pFile,"%s%s",inpC1,inpC2) ; // NpH
    fscanf(pFile,"%s%s",inpC1,inpC2) ; // NpV
    fscanf(pFile,"%s%s",inpC1,inpC2) ; N_src = atoi( inpC2 );
    fscanf(pFile,"%s%s",inpC1,inpC2) ; // Z_src
    fscanf(pFile,"%s%s",inpC1,inpC2) ; // sigma_shape

//    if ( std::log10(M) < 14.0 )
//        return false;
    int  I_bin = 0 ;

    if      ( integ ==  0  )  {   I_bin =  0 ;   }
    else if ( integ > 399. )  {   I_bin =  9 ;   }
    else if ( integ > 252. )  {   I_bin =  8 ;   }
    else if ( integ > 159. )  {   I_bin =  7 ;   }
    else if ( integ > 100. )  {   I_bin =  6 ;   }
    else if ( integ >  63. )  {   I_bin =  5 ;   }
    else if ( integ >  39. )  {   I_bin =  4 ;   }
    else if ( integ >  25. )  {   I_bin =  3 ;   }
    else if ( integ >  15. )  {   I_bin =  2 ;   }
    else                      {   I_bin =  1 ;   }


    if ( std::log10(M) < 14.0 && (
         I_bin >= 9 ||
         I_bin < 2 ) )
        return false ;

    int jCounter = 0 ; // Count for jacknife binning

    // For source read in
    while( fscanf( pFile, "%s%s%s", inpC1, inpC2, inpC3 ) != EOF )
    {

        // Dist, tot, tan
        double d = atof( inpC1 );
        double o = atof( inpC2 );
        double a = atof( inpC3 );

        int    d_bin = std::min(std::max(   int( d / dMax * u.getNbins() )    ,0),u.getNbins()-1);
        int    j_bin =    jCounter % u.getJacknifeBins() ;

                        ++jCounter ;

        // Go in all bin spot
        dist[ d_bin ] += d ;
        gTot[ d_bin ] += o ;
        gTan[ d_bin ] += a ;
        N   [ d_bin ] += 1 ;


        // Jack knife
        for ( int i = 0; i<u.getJacknifeBins(); ++i )
        {

            if ( i != j_bin )
            {

                int myBin = u.getSrcBin( j_bin+1, d_bin );

                dist[ myBin ] += d ;
                gTot[ myBin ] += o ;
                gTan[ myBin ] += a ;
                N   [ myBin ] += 1 ;

            }
        }

    }


    // Track halo statistics
    (*h).setID    ( ID    );
    (*h).setM     ( M     );
    (*h).setC     ( c     );
    (*h).setBA    ( ba    );
    (*h).setCA    ( ca    );
    (*h).setRmax  ( r_max );
    (*h).setInteg ( integ );
    (*h).setPhi   ( phi   );
    (*h).setTheta ( theta );
    (*h).setGamma ( gamma );
    (*h).setAlpha ( alpha );


    return true;
}





bool readMapFile(   FILE      * pFile ,
                    userInfo    u     ,  // User input
                    double    * dist  ,  // Array of distances
                    double    * gTot  ,  // Array of gTot
                    double    * gTan  ,  // Array of gTan
                    double    * totS  ,  // Array of gTot
                    double    * tanS  ,  // Array of gTan
                    haloInfo  * h     )  // Array containing halo info for bins
{

    char   inpC1[35], inpC2[35], inpC3[35], inpC4[35], inpC5[35];

    int i = 0 ;
    // For source read in
    while( fscanf( pFile, "%s%s%s%s%s", inpC1, inpC2, inpC3, inpC4, inpC5 ) != EOF )
    {
      
        // Go in all bin spot
        dist[i] = atof( inpC1 ) ;
        gTot[i] = atof( inpC2 ) ;
        gTan[i] = atof( inpC3 ) ;
        totS[i] = atof( inpC4 ) ;
        tanS[i] = atof( inpC5 ) ;

        ++i;
    }

    return true;
}



// Wrapper function for reading source files, will locate files and invoke reader function
int  readSources(  userInfo    u    ,  // User input
                   double    * d    ,  // Array of distances
                   double    * gTot ,  // Array of gTot
                   double    * gTan ,  // Array of gTan
                   double    * totS ,  // Array counting number in each bin
                   double    * tanS ,  // Array counting number in each bin
                   haloInfo  * h    ,  // Array containing averaged halo info
                   int    startLine )  // Line we last read in, need to read back to this line
{
    // Make sure all output arrays 0
    for ( int i = 0; i < u.getNbins(); ++i )
    {
        gTot[i] = 0 ;
        gTan[i] = 0 ;
        totS[i] = 0 ;
        tanS[i] = 0 ;
        d   [i] = 0 ;
    }



    int      halo_c      = 0     ;  // Count number of halos read
    FILE    *inpFileList         ;  // Contains all the input files
    bool     endOfFile   = true  ;  // Track if we reached the end of the file


    inpFileList = fopen( u.getInputFileF().c_str(), "r" ) ;

    char inpFileName[500] ;


    while ( fscanf( inpFileList,"%s",inpFileName) != EOF ) // Go through file
    {
        ++halo_c ;                                         // Count number of lines read
        if ( halo_c > startLine )                          // If num read greater than start num, search for file
        {

          
          FILE *pFile ;
          pFile = fopen( inpFileName, "r");

          if (pFile!=NULL)                           // If file exists, read it, attempt to read others
            {
              

              // Gets info from file name
              std::string myFileName( inpFileName );
              int index = myFileName.find("I"); (*h).setInteg (          std::stod( myFileName.substr( index+1, 5 ) )   );
                  index = myFileName.find("M"); (*h).setM     ( pow( 10, std::stod( myFileName.substr( index+1, 5 ) ) ) );
                  index = myFileName.find("B"); (*h).setBA    (          std::stod( myFileName.substr( index+1, 5 ) )   );
                  index = myFileName.find("G"); 
              if (index != std::string::npos )
                {
                                                (*h).setGamma (          std::stod( myFileName.substr( index+1, 5 ) )   );
                }


              bool validHalo = readMapFile( pFile, u, d, gTot, gTan, totS, tanS, h );

              fclose( pFile );
              
              if (validHalo)                         // If halo above our mass cutoff
                {

                  endOfFile = false ;                // Indicate we did not reach the end of the file
                  break;
                }

            } // File exists
        }     // Jump forward startLine lines
    }         // Halo_id loop

    fclose( inpFileList );


    if ( endOfFile )
        return -1;


    return halo_c;
}





// Max/min/N M, b, g, i's?
// Sc

// Read parameters not included in paramfile (Nbins, Nthreads, etc.)
void readInpFile(          userInfo  &inpInfo  ,   // Info needed for the rest of the code
                  const std::string  inputFile ){  // Input file name
  FILE *pFile;
  char   inpC1[35],inpC2[35];

  // Attempt to open file, if successful go line by line
  //  finding the variable name and value

  pFile = fopen(inputFile.c_str(),"r");

  if (pFile!=NULL){

    logMessage( std::string( "Reading file: ") + inputFile );

    // Scan variables
    while ( fscanf(pFile,"%s%s",inpC1,inpC2) != EOF ){
      std::string inpS = std::string(inpC1);
           if ( inpS=="N_bins"      ){        inpInfo.setNbins           (        atoi(inpC2) ); }  // Number of bins for profile fitting
      else if ( inpS=="N_binsJac"   ){        inpInfo.setJacknifeBins    (        atoi(inpC2) ); }  // Number of jacknife bins to use
      else if ( inpS=="N_threads"   ){        inpInfo.setNthreads        (        atoi(inpC2) ); }  // Number of omp threads

      else if ( inpS=="N_fitAttempt"){        inpInfo.setMaxFitNum       (        atoi(inpC2) ); }  // Maximum number of times to step the fitter
      else if ( inpS=="N_probes"    ){        inpInfo.setNchrome         (        atoi(inpC2) ); }  // Number of probes to use in fitter
      else if ( inpS=="N_consistent"){        inpInfo.setNConsistent     (        atoi(inpC2) ); }  // Number of times to converge before exiting fit
      else if ( inpS=="fitTolerance"){        inpInfo.setTolerance       (        atof(inpC2) ); }  // How small step size should be before being consistent

      else if ( inpS=="useNoise"    ){        inpInfo.setUseNoise        (        atoi(inpC2) ); }  // Sets whether or not to add artificial noise
      else if ( inpS=="sigmaCrit"   ){        inpInfo.setSigmaCrit       (        atof(inpC2) ); }  // Sigma crit, need to set a default
      else if ( inpS=="shapeNoise"  ){        inpInfo.setShapeNoise      (        atof(inpC2) ); }  // Shape noise to use for source errors, default 0.3

      else if ( inpS=="storedData"  ){        inpInfo.setStoredData      ( std::string(inpC2) ); }  // Location of saved data file, containing bulk data of halos
      else if ( inpS=="fox2012F"    ){        inpInfo.setFoxH2012F       ( std::string(inpC2) ); }  // Location of foxH files we will interpolate over
      else if ( inpS=="fox2123F"    ){        inpInfo.setFoxH2123F       ( std::string(inpC2) ); }  // Location of foxH files we will interpolate over
      else if ( inpS=="outputPath"  ){        inpInfo.setOutputPath      ( std::string(inpC2) ); }  // Directory to place output files in
      else if ( inpS=="inputPath"   ){        inpInfo.setInputPath       ( std::string(inpC2) ); }  // Directory to read input files from

      else if ( inpS=="firstFile"   ){        inpInfo.setFirstFile       (        atoi(inpC2) ); }  // First file halo ID to use
      else if ( inpS== "lastFile"   ){        inpInfo.setLastFile        (        atoi(inpC2) ); }  // Last  file halo ID to use
      else if ( inpS==  "fitType"   ){        inpInfo.setFitType         (        atoi(inpC2) ); }  // Number of free parameters to fit

      else if ( inpS=="doM"         ){        inpInfo.setDoM             (        atoi(inpC2) ); }  // Collapse B&G bins to study M
      else if ( inpS=="doG"         ){        inpInfo.setDoG             (        atoi(inpC2) ); }  // Collapse B and M bins to study G
      else if ( inpS=="doBG"        ){        inpInfo.setDoBG            (        atoi(inpC2) ); }  // Collapse M bins to study B&G


      else{

          // Abort if unrecognized variables

          std::cout << " Couldn't recognize input from " << inputFile <<
                     ": " << inpS << std::endl << std::endl;
          logMessage( std::string("Unrecognized input: ") + inpS );
        exit(1);
      }
    }
    fclose(pFile);
  }
  // Abort if couldn't open file
  else{
                std::cout << "Couldn't find file: " << inputFile << std::endl;
    logMessage( std::string( "Couldn't find file: ") + inputFile );
    logMessage( std::string( "Aborting." ) );
    exit(1);
  }

  logMessage( std::string(  "N_bins    = ") + std::to_string((long long  ) inpInfo.getNbins           () ) +
              std::string("\nN_binsJac = ") + std::to_string((long long  ) inpInfo.getJacknifeBins    () ) +
              std::string("\nN_threads = ") + std::to_string((long long  ) inpInfo.getNthreads        () ) +
              std::string("\nN_fitAtte = ") + std::to_string((long long  ) inpInfo.getMaxFitNum       () ) +
              std::string("\nN_probes  = ") + std::to_string((long long  ) inpInfo.getNchrome         () ) +
              std::string("\nN_consist = ") + std::to_string((long long  ) inpInfo.getNConsistent     () ) +
              std::string("\nfitTolera = ") + std::to_string((long double) inpInfo.getTolerance       () ) +
              std::string("\nshapeNois = ") + std::to_string((long double) inpInfo.getShapeNoise      () ) +
              std::string("\nuseNoise  = ") + std::to_string((long long  ) inpInfo.getUseNoise        () ) +
              std::string("\nfitType   = ") + std::to_string((long long  ) inpInfo.getFitType         () ) +
              std::string("\nFoxH2012F = ") + std::string   (              inpInfo.getFoxH2012F       () ) +
              std::string("\nFoxH2123F = ") + std::string   (              inpInfo.getFoxH2123F       () ) +
              std::string("\noutputPat = ") + std::string   (              inpInfo.getOutputPath      () ) );

}



// Reads the fox H tables, saved in log
einTable readFoxH( userInfo &u, const int fileType ){

  double minX, maxX;
  double minA, maxA;

  int x_bins, a_bins;

  einTable einKappa;

  std::string          myFile = u.getFoxH2012F();//"src/foxH2012.dat";
  if ( fileType == 2 ) myFile = u.getFoxH2123F();//"src/foxH2123.dat";

  FILE *pFile;

  pFile = fopen(myFile.c_str(),"r");

  if (pFile!=NULL){

    logMessage( std::string( "Reading file: ") + myFile );

    fscanf( pFile, "%16lf%16lf%4i",&minX,&maxX,&x_bins);
    fscanf( pFile, "%16lf%16lf%4i",&minA,&maxA,&a_bins);

    // Allocate the files
    einKappa.setX_min( minX );
    einKappa.setX_max( maxX );
    einKappa.setA_min( minA );
    einKappa.setA_max( maxA );

    einKappa   .setBins( a_bins, x_bins );


    logMessage( std::string( "Number of alpha bins: ") + std::to_string( (long long) a_bins ) );
    logMessage( std::string( "Number of x     bins: ") + std::to_string( (long long) x_bins ) );


    for ( int i = 0; i < a_bins; ++i ){ // Each row is a new alpha
    for ( int j = 0; j < x_bins; ++j ){ // Each column is a different x

      double inpVal;

      fscanf( pFile, "%16lf", &inpVal); // Goes across rows, then down columns

      einKappa.setVal( i, j, inpVal );

    }
    }

  } else {

    logMessage( std::string( "Cannot open FoxH file: ") + myFile );

    std::cout << "Couldn't open FoxH file: " << myFile << std:: endl;
    exit(0);

  }

  logMessage( std::string( "FoxH Read in complete" ) );

  return einKappa;
}
