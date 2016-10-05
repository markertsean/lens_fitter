/**
This code is to analyze source files output by the lensing_calculator code, with format dist gtot gtan

Read in file governing bins & such

Search for possible source files
Read each file found, and place in bin, i, m, b, g, r ?

once binned, run fitter

output fits to output files

**/


#include <ctime>
#include <slsimlib.h>
#include <iostream>
#include <cstring>
#include <simpleTree.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <omp.h>
#include <thread>
#include <mutex>
#include <cosmo.h>
#include "grid_maintenance.h"
#include "gridmap.h"

#include <CCfits/CCfits>

//My files
#include "astro_constants.h"
#include "lensing_classes.h"
#include "my_utilities.h"
#include "input_functions.h"
#include "lens_fitter.h"
#include "output_functions.h"
#include "pixelmap_functions.h"


// Name of the code log file
std::string logFileName = "";

// Tables for Einasto interpolation
einTable einKappa    ;
einTable einKappaAvg ;

int main(int arg,char **argv){


  // Initializes the log file, generates logfiles directory
  //  and a file name based on current time
  initLogFile();

  // Default values that will not be changing
  long    seed    = time(NULL); //-1827674;

  logMessage( std::string("Seed = ") + std::to_string( (long long) seed ) );

  srand(seed); // Sets random seed


  //////////////////////////////////
  ///////READ IN USERINFO///////////
  //////////////////////////////////

/*

  // Generates Fox H tables for interpolating over the Einasto profiles

  std::cout << "Reading foxH tables: "                  << std::endl;

  einKappa    = readFoxH( userInput, 1 );
  einKappaAvg = readFoxH( userInput, 2 );

  std::cout << "              Done."  << std::endl << std::endl;

    ////////////////////////////////////////////////////////////
    ///////////////////Generate sources/////////////////////////
    ////////////////////////////////////////////////////////////


    // Generates random source positions, distances, and errors for the shape measurements
    // Random positions are stored in a 1D array, with pixelmap index
    // Distance is also stored in a 1D array

    std::cout << "Generating sources..." << std::endl;

    double *srcErrArr ;
    double *srcDArr   ;
    int    *indexes   ;


    double *srcErrArrTemp;
    double *srcDArrTemp  ;
    int    *indexesTemp  ;

    srcErrArrTemp = new double[ userInput.getNsrc() ];
    srcDArrTemp   = new double[ userInput.getNsrc() ];
    indexesTemp   = new int   [ userInput.getNsrc() ];


    logMessage( std::string("Allocated src arrays of size: ") + std::to_string((long long) userInput.getNsrc()) );

    for (int i = 0; i < userInput.getNsrc(); ++i ){
      srcErrArrTemp[i] = userInput.getShapeNoise();                           // Errors all a temporary 0.3
    }



    // Place in random index locations

    std::cout << "  generating indexes..." << std::endl;

    getRandomSourcesIndexes( indexesTemp, userInput );

    logMessage( std::string("Sources placed") );



    // Generate distances based on index

    std::cout << "  calculating distances..." << std::endl;


    // Comb through sources and only keep those outside our halo radius
    {

      int new_Nsrc = 0;

      // Determine distances of the sources from center of cluster, if within Rmax flag as -1
      new_Nsrc = distArrCalc( srcDArrTemp, indexesTemp, &distMap, userInput.getPhysFOV() / userInput.getAngFOV() * 180 / M_PI, userInput.getNsrc(), myHalo.getRmax() );

      logMessage( std::string("Source distances from center found") );


      srcErrArr = new double [ new_Nsrc ];
      srcDArr   = new double [ new_Nsrc ];
      indexes   = new int    [ new_Nsrc ];

      int validCounter = 0;

      for ( int i = 0; i < userInput.getNsrc(); ++i ){
        if ( srcDArrTemp[i] != -1 ){
          srcErrArr[ validCounter ] = srcErrArrTemp[i];
          srcDArr  [ validCounter ] = srcDArrTemp  [i];
          indexes  [ validCounter ] = indexesTemp  [i];
          ++validCounter;
        }
      }
      userInput.setNsrc( new_Nsrc );

    }

    delete[] srcErrArrTemp ;
    delete[] srcDArrTemp   ;
    delete[] indexesTemp   ;


    std::cout << "Done." << std::endl << std::endl;

    // Write the sources to file
    writeDataPoints( userInput, myHalo, srcDArr, indexes, g_totMap, g_tanMap );

// Remove fitting, it's slow and can do later
/*

    ////////////////////////////////////////////////////////////
    ///////////////////To calculate errors//////////////////////
    //////////////Loop over and jacknife sample/////////////////
    ////////////////////////////////////////////////////////////

    int N_jackbins = userInput.getJacknifeBins() * userInput.getJacknifeBins() ;


                std::cout <<"Using " <<     N_jackbins    << " subsets to calculate errors" << std::endl;
    logMessage( std::string   (            "N_jackbins = ") +
                std::to_string( (long long) N_jackbins    ) );


    densProfile nfwFits[ N_jackbins + 1 ];
    densProfile nfTFits[ N_jackbins + 1 ];
    densProfile einFits[ N_jackbins + 1 ];


                std::cout <<"Generating fits... "  << std::endl;
    logMessage( std::string("Generating fits... ") );


    for ( int  omitIndex = -1; omitIndex < N_jackbins ; ++ omitIndex ) {

      nfwFits[ omitIndex + 1 ].setR_max( myHalo.getRmax() );
      nfTFits[ omitIndex + 1 ].setR_max( myHalo.getRmax() );
      einFits[ omitIndex + 1 ].setR_max( myHalo.getRmax() );
      einFits[ omitIndex + 1 ].setType( 2 );
      nfTFits[ omitIndex + 1 ].setType( 0 );


      ////////////////////////////////////////////////////////////
      ///////////////////Determine radial average/////////////////
      ///////////////////Bin different source vals////////////////
      ////////////////////////////////////////////////////////////



      // Here loop over sources, indicating which index to omit (start at -1)
      // Pass to dist and shear calculators, to indicate ommited index
      if ( omitIndex > -1 ){
                  std::cout <<"  Omitting subset " << omitIndex << std::endl;
      logMessage( std::string("  Omitting subset ") +
                  std::to_string( (long long) omitIndex ) );
      }

      // Array to bin distances, RTS, and their errors
      double    distArr[ userInput.getNbins() ],
                gTanArr[ userInput.getNbins() ],
                gErrArr[ userInput.getNbins() ];


      radialDistAverage( distArr, srcDArr, userInput, indexes, center, omitIndex );

      logMessage( std::string("  Source distances averaged") );

      radialShearAverage( gTanArr, gErrArr, indexes, g_totMap, srcErrArr, srcDArr, userInput, center, omitIndex );

      logMessage( std::string("  Shear values averaged") );

      std::cout << "  Sources averaged" << std::endl;


      //////////////////////////////////////////////////////////
      ////////////////////////FIT PROFILE///////////////////////
      //////////////////////////////////////////////////////////

      // Attempts to fit the density using the radial averages of distance and RTS

                  std::cout <<"  Calculating NFW fit..." << std::endl;
      logMessage( std::string("  Calculating NFW fit...") );
      rollingFitDensProfile( nfwFits[ omitIndex + 1], myHalo, userInput, gTanArr, distArr, gErrArr, cosmo );


                  std::cout <<"  Calculating NFW trunc fit..." << std::endl;
      logMessage( std::string("  Calculating NFW trunc fit...") );
      rollingFitDensProfile( nfTFits[ omitIndex + 1], myHalo, userInput, gTanArr, distArr, gErrArr, cosmo );


                  std::cout <<"  Calculating EIN fit..." << std::endl;
      logMessage( std::string("  Calculating EIN fit...") );
      rollingFitDensProfile( einFits[ omitIndex + 1], myHalo, userInput, gTanArr, distArr, gErrArr, cosmo );

      std::cout << std::endl ;

    }


    // Calculate the jacknife errors

    double nfwErr[3]; // 0 is C, 1 is log(M), 2 is alpha
    double nfTErr[3];
    double einErr[3];

    jacknife( nfwFits, N_jackbins , nfwErr );
    jacknife( nfTFits, N_jackbins , nfTErr );
    jacknife( einFits, N_jackbins , einErr );


                std::cout <<"Done.              " << std::endl;
    logMessage( std::string("Fitting complete"   ));




    writeProfileFits( userInput, myHalo, einFits[0], nfwFits[0], nfTFits[0], einErr, nfwErr, nfTErr, halo_index );


//*/

  exit(0);
  return 0;
}
