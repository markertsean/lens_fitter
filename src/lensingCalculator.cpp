/**
This code is to analyze source files output by the lensing_calculator code, with format dist gtot gtan

Read in file governing bins & such

Search for possible source files
Read each file found, and place in bin, i, m, b, g, r ?

once binned, run fitter

output fits to output files

**/


#include <ctime>
//#include <slsimlib.h>
#include <iostream>
#include <cstring>
//#include <simpleTree.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <omp.h>
#include <thread>
#include <mutex>
//#include <cosmo.h>
//#include "grid_maintenance.h"
//#include "gridmap.h"

//#include <CCfits/CCfits>

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

    std::cout << "Reading files..." << std::endl;

    int N_lines_read = 0 ; // Tracks the number of files read from list of files, upon returning -1, the code exits

    while ( N_lines_read > -1 )
    {

        // Arrays containing avgs, binned
        double *  gTotArr ;
        double *  gTanArr ;
        double *  gTotStd ;
        double *  gTanStd ;
        double *     dArr ;


        haloInfo *     myHalo ;


        gTotArr = new   double[ userInput.getNbins() ] () ;
        gTanArr = new   double[ userInput.getNbins() ] () ;
        gTotStd = new   double[ userInput.getNbins() ] () ;
        gTanStd = new   double[ userInput.getNbins() ] () ;

           dArr = new   double[ userInput.getNbins() ] () ;

         myHalo = new haloInfo ;

         N_lines_read = readSources( userInput, dArr, gTotArr, gTanArr, gTotStd, gTanStd, myHalo, N_lines_read );

        if ( N_lines_read == -1 )
            break;
        /*
        delete gTotArr;
delete gTanArr;
delete gTotStd;
delete gTanStd;
delete dArr;
        //*/
    }
        /*

            ///////////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////
            ////////////////Perform fitting of sources/////////////////////
            ///////////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////

            // Will store density profile fits for each jacknife bin, + full sample as 0th index
            densProfile nfwFits_tot[ N_jackbins + 1 ];
            densProfile nfTFits_tot[ N_jackbins + 1 ];
            densProfile einFits_tot[ N_jackbins + 1 ];

            densProfile nfwFits_tan[ N_jackbins + 1 ];
            densProfile nfTFits_tan[ N_jackbins + 1 ];
            densProfile einFits_tan[ N_jackbins + 1 ];

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
                double *gTot = new double[ userInput.getNbins() ] ();
                double *gTan = new double[ userInput.getNbins() ] ();
                double *dArr = new double[ userInput.getNbins() ] ();
                int    *N_arr= new int   [ userInput.getNbins() ] ();

                for ( int r = 0; r < userInput.getNbins() ; ++r )
                {
                    int index = userInput.getSrcBin( j, r );
                    gTot[r]   = gTotJackArr[index];
                    gTan[r]   = gTanJackArr[index];
                    dArr[r]   =    dJackArr[index];
                   N_arr[r]   =    nJackArr[index];
                }

                // Uncertainty array
                double *eArr = gaussUncertaintyArr(       userInput.getShapeNoise(), N_arr, userInput.getNbins() );

                if ( userInput.getUseNoise() == 1 )
                {
                            addgaussUncertaintyArr( gTot, userInput.getShapeNoise(), N_arr, userInput.getNbins() ); // Adds random amount of noise to bin
                            addgaussUncertaintyArr( gTan, userInput.getShapeNoise(), N_arr, userInput.getNbins() );
                }


                nfwFits_tot[ omitIndex + 1 ].setR_max( (*myHalo).getRmax() );
                nfTFits_tot[ omitIndex + 1 ].setR_max( (*myHalo).getRmax() );
                einFits_tot[ omitIndex + 1 ].setR_max( (*myHalo).getRmax() );
                nfwFits_tan[ omitIndex + 1 ].setR_max( (*myHalo).getRmax() );
                nfTFits_tan[ omitIndex + 1 ].setR_max( (*myHalo).getRmax() );
                einFits_tan[ omitIndex + 1 ].setR_max( (*myHalo).getRmax() );

                nfwFits_tot[ omitIndex + 1 ].setM_enc( (*myHalo).getM   () );
                nfTFits_tot[ omitIndex + 1 ].setM_enc( (*myHalo).getM   () );
                einFits_tot[ omitIndex + 1 ].setM_enc( (*myHalo).getM   () );
                nfwFits_tan[ omitIndex + 1 ].setM_enc( (*myHalo).getM   () );
                nfTFits_tan[ omitIndex + 1 ].setM_enc( (*myHalo).getM   () );
                einFits_tan[ omitIndex + 1 ].setM_enc( (*myHalo).getM   () );

                nfwFits_tot[ omitIndex + 1 ].setC    ( (*myHalo).getC   () );
                nfTFits_tot[ omitIndex + 1 ].setC    ( (*myHalo).getC   () );
                einFits_tot[ omitIndex + 1 ].setC    ( (*myHalo).getC   () );
                nfwFits_tan[ omitIndex + 1 ].setC    ( (*myHalo).getC   () );
                nfTFits_tan[ omitIndex + 1 ].setC    ( (*myHalo).getC   () );
                einFits_tan[ omitIndex + 1 ].setC    ( (*myHalo).getC   () );

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


            char     fileName[500] ;

            sprintf( fileName, "%sDensFitTot_ID%010i_I%06.1f.dat" ,  userInput.getOutputPath().c_str(),
                                                                     (*myHalo).getID()    ,
                                                                     (*myHalo).getInteg() );

            writeProfileFits(       fileName   ,
                                    userInput  ,
                                    *myHalo    ,
                                einFits_tot[0] ,
                                nfwFits_tot[0] ,
                                nfTFits_tot[0] ,
                                 einErr_tot    ,
                                 nfwErr_tot    ,
                                 nfTErr_tot    );

            sprintf( fileName, "%sDensFitTan_ID%010i_I%06.1f.dat" ,  userInput.getOutputPath().c_str(),
                                                                     (*myHalo).getID()    ,
                                                                     (*myHalo).getInteg() );

            writeProfileFits(       fileName   ,
                                    userInput  ,
                                    *myHalo    ,
                                einFits_tan[0] ,
                                nfwFits_tan[0] ,
                                nfTFits_tan[0] ,
                                 einErr_tan    ,
                                 nfwErr_tan    ,
                                 nfTErr_tan    );

            std::cout << std::endl;

        delete[] gTotJackArr ;
        delete[] gTanJackArr ;
        delete[]    dJackArr ;
        delete[]    nJackArr ;
    }
        //*/
  exit(0);
  return 0;

}
