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

    userInfo userInput;

    std::string userFile = "lensUserParams.dat";
    std::cout << "Reading user input" << userFile         << std::endl;

    readInpFile( userInput, userFile );


    // Generates Fox H tables for interpolating over the Einasto profiles

    std::cout << "Reading foxH tables: "                  << std::endl;

    einKappa    = readFoxH( userInput, 1 );
    einKappaAvg = readFoxH( userInput, 2 );

    std::cout << "              Done."  << std::endl << std::endl;

    int N_jackbins = userInput.getJacknifeBins() ;


                std::cout <<"Using " <<     N_jackbins    << " subsets to calculate errors" << std::endl;
    logMessage( std::string   (            "N_jackbins = ") +
                std::to_string( (long long) N_jackbins    ) );


    ////////////////////////////////////////////////////////////
    ///////////////////Read in the sources//////////////////////
    ////////////////////////////////////////////////////////////

    std::cout << "Reading sources..." << std::endl;


    // Arrays containing avgs, binned
    double *  gTotJackArr ;
    double *  gTanJackArr ;
    double *     dJackArr ;
    int    *     nJackArr ; // Source counts in a I, M, B, G, R bin
    int    *       ninArr ; // Halo   counts in a I, M, B, G    bin

    haloInfo * binnedHalo ;


    gTotJackArr = new   double[ userInput.getN_srcJackBin()                        ] () ;
    gTanJackArr = new   double[ userInput.getN_srcJackBin()                        ] () ;
       dJackArr = new   double[ userInput.getN_srcJackBin()                        ] () ;
       nJackArr = new      int[ userInput.getN_srcJackBin()                        ] () ;
         ninArr = new      int[ userInput.getN_srcBin    () / userInput.getNbins() ] () ;

    binnedHalo  = new haloInfo[ userInput.getN_srcBin    () / userInput.getNbins() ] () ;



// read haloinfo into files
// Only reading in halo files
    if (           !readShortFile( userInput, gTotJackArr, gTanJackArr,    dJackArr, nJackArr, ninArr, binnedHalo ) )
    {
        int N_files = readSources( userInput,    dJackArr, gTotJackArr, gTanJackArr, nJackArr, ninArr, binnedHalo ) ;

        if( N_files == 0 )
        {
            std::cout << "Did not read any files" << std::endl;
            exit(1);
        }

            std::cout << "Read " << N_files << " files" << std::endl;

        writeShort(                userInput, gTotJackArr, gTanJackArr,    dJackArr, nJackArr, ninArr, binnedHalo );

            std::cout << "Outputted short file  " << std::endl;


    } else
    {
            std::cout << "Read stored file" << std::endl;
    }


    ////////////////////////////////////////////////////////////
    ///////////////////Loops over source bins///////////////////
    /////////////////////Calculate the fits/////////////////////
    ////////////////////////////////////////////////////////////

/*
Loops:
    j errs
    i=0, i=1
    All M
    b&g looked at

    j errs
    Collapsed b&g into M
    All i
    All M

absorb all halos data into new structure, haloinfo maybe, binned
//*/


    // Loop over each bin, doing jacknife analysis
    // Will generate fits for gTot, gTan, 3 profiles
    for ( int i = 0; i < userInput.getN_IBin(); ++i ){
    for ( int b = 0; b < userInput.getN_BBin(); ++b ){
    for ( int g = 0; g < userInput.getN_GBin(); ++g ){


        // Count number of halos in collapsed bin
        int N_inbin = 0;
        for ( int m = 0; m < userInput.getN_MBin(); ++m )
        {
            int index = userInput.getN_haloBin( i, m, b, g );
            N_inbin += ninArr[index] ;
        }


        // Only analyze bin in halos are present
        if ( N_inbin > 0 )
        {

            std::cout << std::endl << "Collapsing M bins, i = " << i << " b = " << b << " g = " << g << std::endl;


            // Will store density profile fits for each jacknife bin, + full sample as 0th index
            densProfile nfwFits_tot[ N_jackbins + 1 ];
            densProfile nfTFits_tot[ N_jackbins + 1 ];
            densProfile einFits_tot[ N_jackbins + 1 ];

            densProfile nfwFits_tan[ N_jackbins + 1 ];
            densProfile nfTFits_tan[ N_jackbins + 1 ];
            densProfile einFits_tan[ N_jackbins + 1 ];


            haloInfo avgHalo;


            // Loop over jacknife bins, ommiting from the sum whichever jack knife bins
            // if omitindex == -1, full sample
            for ( int j = 0; j < N_jackbins + 1; ++j )
            {

                int omitIndex = j - 1 ;

                // Here loop over sources, indicating which index to omit (start at -1)
                // Pass to dist and shear calculators, to indicate ommited index
                if ( omitIndex > -1 )
                {
                                std::cout <<"  Omitting subset " << omitIndex << std::endl;
                    logMessage( std::string("  Omitting subset ") +
                                std::to_string(    (long long)      omitIndex ) );
                }


                // Arrays to populate by calculating averages at index w/ collapsed M
                double *gTot ;
                double *gTan ;
                double *dArr ;

                             avgMArr( userInput, gTotJackArr, nJackArr, i, b, g, omitIndex, &gTot );
                             avgMArr( userInput, gTanJackArr, nJackArr, i, b, g, omitIndex, &gTan );
                int* N_arr = avgMArr( userInput,    dJackArr, nJackArr, i, b, g, omitIndex, &dArr );

                // Uncertainty array
                double *eArr = gaussUncertaintyArr(       userInput.getShapeNoise(), N_arr, userInput.getNbins() );

                            addgaussUncertaintyArr( gTot, userInput.getShapeNoise(), N_arr, userInput.getNbins() ); // Adds random amount of noise to bin
                            addgaussUncertaintyArr( gTan, userInput.getShapeNoise(), N_arr, userInput.getNbins() );


                // Generates average info we will be using to compare against
                avgHalo = avgMHaloInfo( userInput, binnedHalo, ninArr, i, b, g ) ;

                nfwFits_tot[ omitIndex + 1 ].setR_max( avgHalo.getRmax() );
                nfTFits_tot[ omitIndex + 1 ].setR_max( avgHalo.getRmax() );
                einFits_tot[ omitIndex + 1 ].setR_max( avgHalo.getRmax() );
                nfwFits_tan[ omitIndex + 1 ].setR_max( avgHalo.getRmax() );
                nfTFits_tan[ omitIndex + 1 ].setR_max( avgHalo.getRmax() );
                einFits_tan[ omitIndex + 1 ].setR_max( avgHalo.getRmax() );

                nfwFits_tot[ omitIndex + 1 ].setType( 1 ); // Sets as full  NFW
                nfTFits_tot[ omitIndex + 1 ].setType( 0 ); // Sets as trunc NFW
                einFits_tot[ omitIndex + 1 ].setType( 2 ); // Sets as Einasto
                nfwFits_tan[ omitIndex + 1 ].setType( 1 ); // Sets as full  NFW
                nfTFits_tan[ omitIndex + 1 ].setType( 0 ); // Sets as trunc NFW
                einFits_tan[ omitIndex + 1 ].setType( 2 ); // Sets as Einasto



                // Attempts to fit the density using the radial averages of distance and RTS

                            std::cout <<"  Calculating NFW fit..." << std::endl;
                logMessage( std::string("  Calculating NFW fit...") );
                rollingFitDensProfile( nfwFits_tot[ omitIndex + 1], userInput, gTot, dArr, eArr );
                rollingFitDensProfile( nfwFits_tan[ omitIndex + 1], userInput, gTan, dArr, eArr );


                            std::cout <<"  Calculating NFW trunc fit..." << std::endl;
                logMessage( std::string("  Calculating NFW trunc fit...") );
                rollingFitDensProfile( nfTFits_tot[ omitIndex + 1], userInput, gTot, dArr, eArr );
                rollingFitDensProfile( nfTFits_tan[ omitIndex + 1], userInput, gTan, dArr, eArr );


                            std::cout <<"  Calculating EIN fit..." << std::endl;
                logMessage( std::string("  Calculating EIN fit...") );
                rollingFitDensProfile( einFits_tot[ omitIndex + 1], userInput, gTot, dArr, eArr );
                rollingFitDensProfile( einFits_tan[ omitIndex + 1], userInput, gTan, dArr, eArr );

                std::cout << std::endl ;




                delete [] eArr  ;
                delete [] gTot  ;
                delete [] gTan  ;
                delete [] dArr  ;
                delete [] N_arr ;

            }
                        std::cout << " Generating jack knife errors..." << std::endl;

            // Calculate the jacknife errors from above fits
            double nfwErr_tot[4]; // 0 is C, 1 is log(M), 2 is alpha, 3 is Rvir
            double nfTErr_tot[4];
            double einErr_tot[4];

            jacknife( nfwFits_tot, N_jackbins , nfwErr_tot );
            jacknife( nfTFits_tot, N_jackbins , nfTErr_tot );
            jacknife( einFits_tot, N_jackbins , einErr_tot );

            double nfwErr_tan[4]; // 0 is C, 1 is log(M), 2 is alpha, 3 is Rvir
            double nfTErr_tan[4];
            double einErr_tan[4];

            jacknife( nfwFits_tan, N_jackbins , nfwErr_tan );
            jacknife( nfTFits_tan, N_jackbins , nfTErr_tan );
            jacknife( einFits_tan, N_jackbins , einErr_tan );


                        std::cout <<"Done.              " << std::endl;
            logMessage( std::string("Fitting complete"   ) ) ;

            char     fileName[100] ;

            double        iBin = -6;
            if ( i != 0 ) iBin = userInput.getI_bin( i-1 );

            sprintf( fileName, "%sDensFitTot_M_B%05.3f_G%05.3f_I%06.1f.dat" ,  userInput.getOutputPath().c_str(),
                                                                               userInput.getB_bin  ( b ),
                                                                               userInput.getG_bin  ( g ),
                                                                               pow( 10,      iBin  )   );

            writeProfileFits(       fileName   ,
                                    userInput  ,
                                    avgHalo    ,
                                einFits_tot[0] ,
                                nfwFits_tot[0] ,
                                nfTFits_tot[0] ,
                                 einErr_tot    ,
                                 nfwErr_tot    ,
                                 nfTErr_tot    ,
                                N_inbin        );

            sprintf( fileName, "%sDensFitTan_M_B%05.3f_G%05.3f_I%06.1f.dat" ,  userInput.getOutputPath().c_str(),
                                                                               userInput.getB_bin  ( b ),
                                                                               userInput.getG_bin  ( g ),
                                                                               pow( 10,      iBin  )   );

            writeProfileFits(       fileName   ,
                                    userInput  ,
                                    avgHalo    ,
                                einFits_tan[0] ,
                                nfwFits_tan[0] ,
                                nfTFits_tan[0] ,
                                 einErr_tan    ,
                                 nfwErr_tan    ,
                                 nfTErr_tan    ,
                                N_inbin        );

            std::cout << std::endl;

        }
    }
    }
    }


/*

      //////////////////////////////////////////////////////////
      ////////////////////////FIT PROFILE///////////////////////
      //////////////////////////////////////////////////////////

      nfwFits[ omitIndex + 1 ].setR_max( myHalo.getRmax() );
      nfTFits[ omitIndex + 1 ].setR_max( myHalo.getRmax() );
      einFits[ omitIndex + 1 ].setR_max( myHalo.getRmax() );
      einFits[ omitIndex + 1 ].setType( 2 );
      nfTFits[ omitIndex + 1 ].setType( 0 );


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
