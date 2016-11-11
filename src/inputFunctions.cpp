#include <cstring>
#include <cmath>
#include <vector>
#include <CCfits/CCfits>
#include <slsimlib.h>
#include <stdio.h>
#include "lensing_classes.h"


bool readShortFile( userInfo      u ,
                    double   * gTot ,
                    double   * gTan ,
                    double   * d    ,
                    int      * n_s  ,
                    int      * n_h  ,
                    haloInfo * bh   )
{
    std::string storedData = "storedData.dat";

    std::ifstream file( storedData.c_str() );

    if ( ! file.is_open() )
        return false;



    // Reads binned halo info, averages for the bins
    for ( int i = 0; i < u.getN_IBin(); ++i ){
    for ( int m = 0; m < u.getN_MBin(); ++m ){
    for ( int b = 0; b < u.getN_BBin(); ++b ){
    for ( int g = 0; g < u.getN_GBin(); ++g ){

        int    H_bin = u.getN_haloBin( i, m, b, g    );
        double  temp ;

        file >> temp ; bh[ H_bin ].setM     ( temp );
        file >> temp ; bh[ H_bin ].setC     ( temp );
        file >> temp ; bh[ H_bin ].setBA    ( temp );
        file >> temp ; bh[ H_bin ].setCA    ( temp );
        file >> temp ; bh[ H_bin ].setRmax  ( temp );
        file >> temp ; bh[ H_bin ].setPhi   ( temp );
        file >> temp ; bh[ H_bin ].setTheta ( temp );
        file >> temp ; bh[ H_bin ].setGamma ( temp );
        file >> temp ; bh[ H_bin ].setAlpha ( temp );

    }
    }
    }
    }



    // Reads N_halo info
    for ( int i = 0; i < u.getN_IBin(); ++i ){
    for ( int m = 0; m < u.getN_MBin(); ++m ){
    for ( int b = 0; b < u.getN_BBin(); ++b ){
    for ( int g = 0; g < u.getN_GBin(); ++g ){
        int k = u.getN_haloBin( i, m, b, g );
        file >>  n_h[ k ];
    }
    }
    }
    }

    // Reads tables
    for ( int j = 0; j < u.getN_JBin(); ++j ){
    for ( int i = 0; i < u.getN_IBin(); ++i ){
    for ( int m = 0; m < u.getN_MBin(); ++m ){
    for ( int b = 0; b < u.getN_BBin(); ++b ){
    for ( int g = 0; g < u.getN_GBin(); ++g ){

        int k = u.getSrcBin( j, i, m, b, g, 0 );

        // Outputs each row as group
        for ( int r = 0; r < u.getNbins(); ++r ){  file >>  n_s[ k + r ] ; }
        for ( int r = 0; r < u.getNbins(); ++r ){  file >>    d[ k + r ] ; }
        for ( int r = 0; r < u.getNbins(); ++r ){  file >> gTan[ k + r ] ; }
        for ( int r = 0; r < u.getNbins(); ++r ){  file >> gTot[ k + r ] ; }


    }
    }
    }
    }
    }


    file.close() ;

    return true;
}

void readSourceFile(   FILE      * pFile ,
                       userInfo    u     ,  // User input
                       double    * dist  ,  // Array of distances
                       double    * gTot  ,  // Array of gTot
                       double    * gTan  ,  // Array of gTan
                       int       * N     ,  // Array counting number in each bin
                       int       * N_h   ,  // Array counting number of halos in bins
                       haloInfo  * h     ,  // Array containing halo info for bins
                       int         I_bin )
{

    char   inpC1[35],inpC2[35], inpC3[35];
    double dMax ;
    int    N_src;

    double M, ba, gamma, c, r_max, ca, phi, theta, alpha;


    fscanf(pFile,"%s%s",inpC1,inpC2) ; // ID
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
    fscanf(pFile,"%s%s",inpC1,inpC2) ; // Integ
    fscanf(pFile,"%s%s",inpC1,inpC2) ; // IntegM
    fscanf(pFile,"%s%s",inpC1,inpC2) ; dMax   = atof( inpC2 ) / 2.0 * std::sqrt(2.); // Go through long axis of image
    fscanf(pFile,"%s%s",inpC1,inpC2) ; // NpH
    fscanf(pFile,"%s%s",inpC1,inpC2) ; // NpV
    fscanf(pFile,"%s%s",inpC1,inpC2) ; N_src = atoi( inpC2 );
    fscanf(pFile,"%s%s",inpC1,inpC2) ; // Z_src
    fscanf(pFile,"%s%s",inpC1,inpC2) ; // sigma_shape

    if ( phi < -10 ) return;


         I_bin   = I_bin + 1;
    int  M_bin   = std::min(std::max(   int( ( log10( M )    - u.getM_minBin() ) / ( u.getM_maxBin() - u.getM_minBin() ) * u.getN_MBin() )    ,0),u.getN_MBin()-1);
    int  B_bin   = std::min(std::max(   int( (        ba     - u.getB_minBin() ) / ( u.getB_maxBin() - u.getB_minBin() ) * u.getN_BBin() )    ,0),u.getN_BBin()-1);
    int  G_bin   = std::min(std::max(   int( (        gamma  - u.getG_minBin() ) / ( u.getG_maxBin() - u.getG_minBin() ) * u.getN_GBin() )    ,0),u.getN_GBin()-1);
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

        int myBin = u.getSrcBin( j_bin, I_bin, M_bin, B_bin, G_bin, d_bin );

        dist[ myBin ] += d ;
        gTot[ myBin ] += o ;
        gTan[ myBin ] += a ;
        N   [ myBin ] += 1 ;
    }

    // Counts the number of halos in this bin
    int  H_bin   = u.getN_haloBin( I_bin, M_bin, B_bin, G_bin ) ;
    N_h[ H_bin ] = N_h[ H_bin ] + 1 ;

    // Generates running average of bin parameters
    // First time through, just save
    // Otherwise, weight previous average to generate running
    if ( N_h[H_bin] == 1 )
    {
        h[H_bin].setM     ( M     );
        h[H_bin].setC     ( c     );
        h[H_bin].setBA    ( ba    );
        h[H_bin].setCA    ( ca    );
        h[H_bin].setRmax  ( r_max );
        h[H_bin].setPhi   ( phi   );
        h[H_bin].setTheta ( theta );
        h[H_bin].setGamma ( gamma );
        h[H_bin].setAlpha ( alpha );
    } else
    {
        int N = N_h[ H_bin ] ;

        h[H_bin].setM     (   M       +  h[H_bin].getM     ()  );
        h[H_bin].setC     (   c       +  h[H_bin].getC     ()  );
        h[H_bin].setBA    (   ba      +  h[H_bin].getBA    ()  );
        h[H_bin].setCA    (   ca      +  h[H_bin].getCA    ()  );
        h[H_bin].setRmax  (   r_max   +  h[H_bin].getRmax  ()  );
        h[H_bin].setPhi   (   phi     +  h[H_bin].getPhi   ()  );
        h[H_bin].setTheta (   theta   +  h[H_bin].getTheta ()  );
        h[H_bin].setGamma (   gamma   +  h[H_bin].getGamma ()  );
        h[H_bin].setAlpha (   alpha   +  h[H_bin].getAlpha ()  );
    }


}


// Wrapper function for reading source files, will locate files and invoke reader function
int  readSources(  userInfo    u    ,  // User input
                   double    * d    ,  // Array of distances
                   double    * gTot ,  // Array of gTot
                   double    * gTan ,  // Array of gTan
                   int       * N    ,  // Array counting number in each bin
                   int       * N_h  ,  // Array containing halo count in each bin
                   haloInfo  * h    )  // Array containing averaged halo info

{
    // Make sure all output arrays 0
    for ( int i = 0; i < u.getN_srcBin(); ++i )
    {
        gTot[i] = 0 ;
        gTan[i] = 0 ;
        d   [i] = 0 ;
        N   [i] = 0 ;
    }
    for ( int i = 0; i < u.getN_srcBin()/u.getNbins() ; ++i )
    {
        N_h [i] = 0 ;
    }



    int     halo_id = u.getFirstFile() ;  // Halo ids of files we are looking for
    int     halo_c  = 0                ;  // Count number of halos read

    char inputFile[100];
    FILE *pFile;


    while ( halo_id < u.getLastFile() ) // Check the id numbers for at least base halo file
    {

        sprintf(       inputFile, "%sHalo_%010li_%06.1f_Sources.dat", u.getInputPath().c_str(), halo_id, 0.0 );
        pFile = fopen( inputFile, "r");

        if (pFile!=NULL)               // If file exists, read it, attempt to read others
        {

                readSourceFile( pFile, u, d, gTot, gTan, N, N_h, h, -1 ); // -1 indicates no integration

                fclose( pFile );

            for ( int i = 0; i < u.getN_IBin()-1; ++ i )
            {
                sprintf(       inputFile,  "%sHalo_%010li_%06.1f_Sources.dat", u.getInputPath().c_str(), halo_id, pow( 10, u.getI_bin( i ) ) );
                pFile = fopen( inputFile,  "r");


                if ( pFile != NULL ){
                readSourceFile( pFile, u, d, gTot, gTan, N, N_h, h, i );
                fclose( pFile );
                }
            }

            ++halo_c;
        } // File exists

        ++halo_id;
    }  // Halo_id loop

    // Current gTot, tan, values are sums, makes them averages
    for ( int i = 0; i < u.getN_srcJackBin(); ++i )
    {
        if ( N[i] > 0 )
        {
            gTot[i] = gTot[i] / N[i] ;
            gTan[i] = gTan[i] / N[i] ;
            d   [i] = d   [i] / N[i] ;
        }
    }

    for ( int i = 0; i < u.getN_srcBin() / u.getNbins() ; ++i )
    {
        if ( N_h[i] > 0 )
        {
            h[i].setM     (   h[i].getM     () / N_h[i] );
            h[i].setC     (   h[i].getC     () / N_h[i] );
            h[i].setBA    (   h[i].getBA    () / N_h[i] );
            h[i].setCA    (   h[i].getCA    () / N_h[i] );
            h[i].setRmax  (   h[i].getRmax  () / N_h[i] );
            h[i].setPhi   (   h[i].getPhi   () / N_h[i] );
            h[i].setTheta (   h[i].getTheta () / N_h[i] );
            h[i].setGamma (   h[i].getGamma () / N_h[i] );
            h[i].setAlpha (   h[i].getAlpha () / N_h[i] );

        }
    }


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

      else if ( inpS=="sigmaCrit"   ){        inpInfo.setSigmaCrit       (        atof(inpC2) ); }  // Sigma crit, need to set a default
      else if ( inpS=="shapeNoise"  ){        inpInfo.setShapeNoise      (        atof(inpC2) ); }  // Shape noise to use for source errors, default 0.3

      else if ( inpS=="fox2012F"    ){        inpInfo.setFoxH2012F       ( std::string(inpC2) ); }  // Location of foxH files we will interpolate over
      else if ( inpS=="fox2123F"    ){        inpInfo.setFoxH2123F       ( std::string(inpC2) ); }  // Location of foxH files we will interpolate over
      else if ( inpS=="outputPath"  ){        inpInfo.setOutputPath      ( std::string(inpC2) ); }  // Directory to place output files in
      else if ( inpS=="inputPath"   ){        inpInfo.setInputPath       ( std::string(inpC2) ); }  // Directory to read input files from

      else if ( inpS=="firstFile"   ){        inpInfo.setFirstFile       (        atoi(inpC2) ); }  // First file halo ID to use
      else if ( inpS== "lastFile"   ){        inpInfo.setLastFile        (        atoi(inpC2) ); }  // Last  file halo ID to use

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
