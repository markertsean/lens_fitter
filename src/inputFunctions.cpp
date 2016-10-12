#include <cstring>
#include <CCfits/CCfits>
#include <slsimlib.h>
#include <lensing_classes.h>






/*
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

  fprintf( pFile , "IntegLength %10.6f\n", u.getIntegLength() );
  fprintf( pFile , "IntegMass   %14.8e\n", u.getImageMass() );
  fprintf( pFile , "FOV         %10.6f\n", u.getPhysFOV() );
  fprintf( pFile , "N_pixH      %10i\n"  , u.getNpixH  () );
  fprintf( pFile , "N_pixV      %10i\n"  , u.getNpixV  () );

  fprintf( pFile , "N_src       %10i\n"  , u.getNsrc      () );
  fprintf( pFile , "Z_src       %10.6f\n", u.getSourceZ   () );
  fprintf( pFile , "sigma_shape %10.6f\n", u.getShapeNoise() );


avgVal += gMap.getValue( k * u.getNpixH() + j );
avgVal2+= tMap.getValue( k * u.getNpixH() + j );
n_srcs += 1;

fprintf( pFile , "%10.6f %14.8e %14.8e\n", dArr[i], avgVal / n_srcs, avgVal2 / n_srcs );
//*/

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
