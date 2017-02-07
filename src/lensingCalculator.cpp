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
#include <cmath>
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
        double *  gTotArr    ;
        double *  gTanArr    ;
        double *  gTotStdArr ;
        double *  gTanStdArr ;
        double *     dArr    ;


        haloInfo *     myHalo ;


        gTotArr    = new double[ userInput.getN_list() ] () ;
        gTanArr    = new double[ userInput.getN_list() ] () ;
        gTotStdArr = new double[ userInput.getN_list() ] () ;
        gTanStdArr = new double[ userInput.getN_list() ] () ;

           dArr    = new double[ userInput.getN_list() ] () ;

         myHalo    = new haloInfo ;

         N_lines_read = readSources( userInput, dArr, gTotArr, gTanArr, gTotStdArr, gTanStdArr, myHalo, N_lines_read );

        if ( N_lines_read == -1 )
            break;

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
                double *gTotS= new double[ userInput.getNbins() ] ();
                double *gTanS= new double[ userInput.getNbins() ] ();
                double *   d = new double[ userInput.getNbins() ] ();

                // Get size of the arrays
                for ( int r = 0; r < userInput.getNbins() ; ++r )
                {

                    if (( r * userInput.getJacknifeBins() / userInput.getNbins() ) != omitIndex )
                    {
                        gTot [r]   =     gTotArr   [r];
                        gTan [r]   =     gTanArr   [r];
                        gTotS[r]   = (gTotStdArr[r]/gTotArr[r]);
                        gTanS[r]   = (gTanStdArr[r]/gTanArr[r]);
                        d    [r]   =        dArr   [r];
                    }
                }

                nfwFits_tot[ omitIndex + 1 ].setType( 1 ); // Sets as full  NFW
                nfTFits_tot[ omitIndex + 1 ].setType( 0 ); // Sets as trunc NFW
                einFits_tot[ omitIndex + 1 ].setType( 2 ); // Sets as Einasto
                nfwFits_tan[ omitIndex + 1 ].setType( 1 ); // Sets as full  NFW
                nfTFits_tan[ omitIndex + 1 ].setType( 0 ); // Sets as trunc NFW
                einFits_tan[ omitIndex + 1 ].setType( 2 ); // Sets as Einasto


                // Attempts to fit the density using the radial averages of distance and RTS

                            std::cout <<"  Calculating NFW fit..." << std::endl;
                logMessage( std::string("  Calculating NFW fit...") );
                rollingFitDensProfile( nfwFits_tot[ omitIndex + 1], userInput, gTot, d, gTotStdArr );
                rollingFitDensProfile( nfwFits_tan[ omitIndex + 1], userInput, gTan, d, gTanStdArr );


                            std::cout <<"  Calculating NFW trunc fit..." << std::endl;
                logMessage( std::string("  Calculating NFW trunc fit...") );
                rollingFitDensProfile( nfTFits_tot[ omitIndex + 1], userInput, gTot, d, gTotStdArr );
                rollingFitDensProfile( nfTFits_tan[ omitIndex + 1], userInput, gTan, d, gTanStdArr );

                            std::cout <<"  Calculating EIN fit..." << std::endl;
                logMessage( std::string("  Calculating EIN fit...") );
                rollingFitDensProfile( einFits_tot[ omitIndex + 1], userInput, gTot, d, gTotStdArr );
                rollingFitDensProfile( einFits_tan[ omitIndex + 1], userInput, gTan, d, gTanStdArr );

                std::cout << std::endl ;

                delete [] gTot  ;
                delete [] gTan  ;
                delete [] gTotS ;
                delete [] gTanS ;
                delete [] d     ;


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


            if ( (*myHalo).getGamma() < 0.001 )
            {

                sprintf( fileName, "%sDensFits_I%05.1f_M%05.2f_B%05.3f.dat",  userInput.getOutputPath().c_str(),
                                                                              (*myHalo).getInteg() ,
                                                                   std::log10((*myHalo).getM())    ,
                                                                              (*myHalo).getBA()    );
            }else
            {
                sprintf( fileName, "%sDensFits_I%05.1f_M%05.2f_B%05.3f_G%05.3f.dat",  userInput.getOutputPath().c_str(),
                                                                                      (*myHalo).getInteg() ,
                                                                           std::log10((*myHalo).getM())    ,
                                                                                      (*myHalo).getBA()    ,
                                                                                      (*myHalo).getGamma() );
            }

            writeProfileFits(       fileName   ,
                                    userInput  ,
                                    *myHalo    ,
                                einFits_tot[0] ,
                                nfwFits_tot[0] ,
                                nfTFits_tot[0] ,
                                 einErr_tot    ,
                                 nfwErr_tot    ,
                                 nfTErr_tot    ,
                                einFits_tan[0] ,
                                nfwFits_tan[0] ,
                                nfTFits_tan[0] ,
                                 einErr_tan    ,
                                 nfwErr_tan    ,
                                 nfTErr_tan    ,
                                          dArr ,
                                       gTotArr ,
                                    gTotStdArr ,
                                       gTanArr ,
                                    gTanStdArr );

            std::cout << std::endl;

        delete[]    gTotArr ;
        delete[]    gTanArr ;
        delete[] gTotStdArr ;
        delete[] gTanStdArr ;
        delete[]       dArr ;
        //*/
    }
  exit(0);
  return 0;

}
