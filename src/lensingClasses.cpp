#include <math.h>
#include "lensing_classes.h"


haloInfo::haloInfo(){

       z=-1.0;        // Redshift
       m=-1.0;        // Mass
       c=-1.0;        // Concentration
    rmax=-1.0;        // Rvir, may be R200
      ba=-1.0;        // b/a axis ratio
      ca=-1.0;        // c/a axis ratio
     phi=-1.0;        // Orientation (out of the page?)
   theta=-1.0;        // Orientation (on the page?)
   alpha=-1.0;        // Orientation xy plane
   gamma=-1.0;        // Orientation along z
      id=-1  ;        // Halo ID number
}



userInfo::userInfo(){

  I_step   =  0.20026     ; // Step in I to make in log
  I_minBin =  1.0        ; // Bin min and max values for input sources
  I_maxBin = 2.6020599913279625  ;
  M_minBin = 13.0        ;
  M_maxBin = 16.0        ;
  B_minBin =  0.5        ;
  B_maxBin =  1.0        ;
  G_minBin =  0.0        ;
  G_maxBin =  0.5 * M_PI ;

  I_Nbins  =  9 ;
  M_Nbins  =  6 ;
  B_Nbins  =  5 ;
  G_Nbins  =  4 ;

  firstFile = 0    ;
  lastFile  = 1000000000 ;

  fox2012F = "src/foxH2012.dat"; // FoxH files to read
  fox2123F = "src/foxH2123.dat";
  inputFileFile = "recoveredFiles.txt";

  outputPath = "data/";
  storedFile = "storedData.dat";

    sigmaC     =  2.991203e+15 ; // Critical surface density in M_sun/Mpc^2
   shapeNoise  =  0.3          ; // Intrinsic shape noise in the sources
   N_bins      = 20            ; // Number of bins for radial averaging
   num_threads =  1            ; // Number of threads for parallel processing


      cMin =  2.5;     // Range of concentration values to fit
      cMax =  7.5;
      rMin =  0.3;     // Range of R_max values to fit, Mpc
      rMax =  2.0;
      mMin = 12.5;     // Range of mass values to fit
      mMax = 16.5;
  alphaMin =  0.15;    // Range of alpha values to fit
  alphaMax =  0.23;


  maxFitAttempts = 2e2   ; // Maximum number of steps to roll ball, or times to reproduce
   N_chromosomes = 1e3   ; // Number of chromosomes or balls

      consistent = 2e1   ; // Number of steps to converge before accepting
       tolerance = 1e-4  ; // If difference between steps less than this, converged


      N_chiTrack = 1e1   ; // Number of chromosomes to check for convergence
       mutChance = 1e-2  ; // Likelihood of converging
      avgTestVal = 0.7   ; // Criteria for reproducing

  jacknifeBins = 10;

    useNoise   =  1;

    fitType    =  1;
    doM        =  1;
    doG        =  1;
    doBG       =  1;
}



densProfile::densProfile(){
  concentration = -1.0;
  alpha         = -1.0;
  r_max         = -1.0;
  M_enc         = -1.0;
  type          =    1;
}


// If passed an argument, Einasto profile
// Alternatively can change type later
densProfile::densProfile( double inpA ){
  concentration = -1.0;
  alpha         = inpA;
  r_max         = -1.0;
  M_enc         = -1.0;
  type          =    2;
}


// Msun/Mpc^3
double densProfile::getRho_o() const {

  if ( M_enc > 0 && r_max > 0 && concentration > 0){ // All needed parameters def

    if ( type != 2 ){ // NFW
                                         return  M_enc / (4. * M_PI * r_max*r_max*r_max )  *
                                                            concentration   * concentration       * concentration /
                                                 ( log( 1 + concentration ) - concentration / ( 1 + concentration) );
    } else {          // Einasto

      double rs = r_max / concentration;

                                         return  M_enc * alpha / (   4. * M_PI    *    rs * rs * rs   *
                                                    exp(             2. / alpha ) *
                                                    pow( alpha / 2., 3. / alpha ) *
                                                    tgamma(          3. / alpha ) );
    }
  }

  return -1;
}








// Log file stuff



// Initializes log file, making directory and file name
void initLogFile(){

  // Sets up the times
  int execution_start = clock();
  time_t      nowTime =      time(        0 );
  tm       *startTime = localtime( &nowTime );

  // Logfile name & directory
  char logFileNameC[100];
  struct stat sb;
  char str[] = "mkdir logfiles";

  // Create logfiles directory, if one does not exist
  if ( stat( "logfiles/", &sb) != 0 ){
    system( str );
  }

  // Log file name is mostly the date
  sprintf( logFileNameC, "logfiles/lensCalc.%4i.%02i.%02i.%02i.%02i.%02i.log",
    (*startTime).tm_year+1900,
    (*startTime).tm_mon ,
    (*startTime).tm_mday,
    (*startTime).tm_hour,
    (*startTime).tm_min ,
    (*startTime).tm_sec );

  logFileName= std::string(logFileNameC) ;

  // First line of the log file is the time
  logMessage( (std::string(                "Code initialized at ")+
               std::to_string( (long long) (*startTime).tm_hour  )+
               std::string(                ":"                   )+
               std::to_string( (long long) (*startTime).tm_min   )+
               std::string(                ":"                   )+
               std::to_string( (long long) (*startTime).tm_sec   )));


}


// Generates log files, printing text to file
void logMessage( const std::string &text ){

    std::ofstream log_file(  logFileName, std::ios_base::out | std::ios_base::app );
    log_file << text << std::endl;

}

