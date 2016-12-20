#ifndef MY_LENSING_CLASSES
#define MY_LENSING_CLASSES

#include <iostream>
#include <fstream>
#include <cstring>
#include <ctime>
#include <cosmo.h>
#include <math.h>
#include <sys/stat.h>
#include <astro_constants.h>




// Log file name
extern std::string logFileName;


// Generates log files, printing text to file
void logMessage( const std::string &text );

// Initializes log file, making directory and file name
void initLogFile();



// Holds info on lens, source, any halo for easy access
class haloInfo{
  public:

//    inline haloInfo();
    haloInfo();

    void setZ     ( double inpZ ){ z     = inpZ; }
    void setM     ( double inpM ){ m     = inpM; }
    void setC     ( double inpC ){ c     = inpC; }
    void setRmax  ( double inpR ){ rmax  = inpR; }
    void setBA    ( double inpF ){ ba    = inpF; }
    void setCA    ( double inpF ){ ca    = inpF; }
    void setPhi   ( double inpF ){ phi   = inpF; }
    void setTheta ( double inpF ){ theta = inpF; }
    void setAlpha ( double inpF ){ alpha = inpF; }
    void setGamma ( double inpF ){ gamma = inpF; }
    void setID    ( long   inpI ){ id    = inpI; }

    double      getZ     () const { return z     ; }
    double      getM     () const { return m     ; }
    double      getC     () const { return c     ; }
    double      getRmax  () const { return rmax  ; }
    double      getBA    () const { return ba    ; }
    double      getCA    () const { return ca    ; }
    double      getPhi   () const { return phi   ; }
    double      getTheta () const { return theta ; }
    double      getAlpha () const { return alpha ; }
    double      getGamma () const { return gamma ; }
    long        getID    () const { return id    ; }


  private:
    double       z;             // redshift
    double       m;             // mass Msun
    double       c;             // concentration
    double    rmax;             // max radius Mpc
    double      ba;             // b/a axis ratio
    double      ca;             // c/a axis ratio
    double     phi;             // y/z orientation angle
    double   theta;             // x/z orientation angle

    double   alpha;             // Orientation on xy plane, +y 0
    double   gamma;             // Orientation along z, +z 0

    long        id;             // id number of halo


    void   unDefVar( std::string inpS ){
      std::cerr<<"WARNING: variable "<<inpS<<" undefined in <haloInfo>"<<std::endl;
    }
};



// Holds info on user input values
class userInfo{
  public:

    userInfo();

    void setJacknifeBins    ( int    inpI ) {   jacknifeBins = inpI ; }
    void setShapeNoise      ( double inpF ) {     shapeNoise = inpF ; }

    void setMaxFitNum       ( int    inpI ) { maxFitAttempts = inpI ; }
    void setNConsistent     ( int    inpI ) {  consistent    = inpI ; }
    void setTolerance       ( double inpF ) {      tolerance = inpF ; }
    void setMutChance       ( double inpF ) {      mutChance = inpF ; }
    void setTestVal         ( double inpF ) {     avgTestVal = inpF ; }
    void setNtrack          ( int    inpI ) {  N_chiTrack    = inpI ; }
    void setNchrome         ( int    inpI ) {  N_chromosomes = inpI ; }

    void setChiMin          ( double inpF ) {           cMin = inpF ; }
    void setChiMax          ( double inpF ) {           cMax = inpF ; }
    void setMassMin         ( double inpF ) {           mMin = inpF ; }
    void setMassMax         ( double inpF ) {           mMax = inpF ; }
    void setRMinFit         ( double inpF ) {           rMin = inpF ; }
    void setRMaxFit         ( double inpF ) {           rMax = inpF ; }
    void setConMin          ( double inpF ) {           cMin = inpF ; }
    void setConMax          ( double inpF ) {           cMax = inpF ; }
    void setAlphaMin        ( double inpF ) {       alphaMin = inpF ; }
    void setAlphaMax        ( double inpF ) {       alphaMax = inpF ; }

    void setNbins           ( int    inpI ) {   N_bins       = inpI ; }
    void setNthreads        ( int    inpI ) {   num_threads  = inpI ; }
    void setUseNoise        ( int    inpI ) {     useNoise   = inpI ; }
    void setSigmaCrit       ( double inpF ) {      sigmaC    = inpF ; }

    void setI_BinMin        ( double inpF ) {     I_minBin   = inpF ; }
    void setM_BinMin        ( double inpF ) {     M_minBin   = inpF ; }
    void setB_BinMin        ( double inpF ) {     B_minBin   = inpF ; }
    void setG_BinMin        ( double inpF ) {     G_minBin   = inpF ; }
    void setI_step          ( double inpF ) {     I_step     = inpF ; }
    void setM_BinMax        ( double inpF ) {     M_maxBin   = inpF ; }
    void setB_BinMax        ( double inpF ) {     B_maxBin   = inpF ; }
    void setG_BinMax        ( double inpF ) {     G_maxBin   = inpF ; }

    void setN_IBin          ( int    inpI ) {     M_Nbins    = inpI ; }
    void setN_MBin          ( int    inpI ) {     M_Nbins    = inpI ; }
    void setN_BBin          ( int    inpI ) {     M_Nbins    = inpI ; }
    void setN_GBin          ( int    inpI ) {     M_Nbins    = inpI ; }

    void setFoxH2012F       ( std::string inpS ) {  fox2012F = inpS ; }
    void setFoxH2123F       ( std::string inpS ) {  fox2123F = inpS ; }
    void setOutputPath      ( std::string inpS ) {outputPath = inpS ; }
    void setInputPath       ( std::string inpS ) { inputPath = inpS ; }

    void setFirstFile       ( int    inpI ) {     firstFile  = inpI ; }
    void setLastFile        ( int    inpI ) {      lastFile  = inpI ; }

    void setFitType         ( int    inpI ) {        fitType = inpI ; }

    int    getFitType         () const { return fitType        ; }
    int    getFirstFile       () const { return firstFile      ; }
    int    getLastFile        () const { return lastFile       ; }
    int    getJacknifeBins    () const { return jacknifeBins   ; }
    int    getNbins           () const { return N_bins         ; }
    int    getNthreads        () const { return num_threads    ; }
    int    getNchrome         () const { return N_chromosomes  ; }
    int    getNtrack          () const { return N_chiTrack     ; }
    int    getMaxFitNum       () const { return maxFitAttempts ; }
    int    getNConsistent     () const { return consistent     ; }
    int    getN_JBin          () const { return jacknifeBins   ; }
    int    getN_IBin          () const { return I_Nbins+1      ; }
    int    getN_MBin          () const { return M_Nbins        ; }
    int    getN_BBin          () const { return B_Nbins        ; }
    int    getN_GBin          () const { return G_Nbins        ; }
    int    getN_srcBin        () const { return N_bins * ( I_Nbins+1 ) * M_Nbins * B_Nbins * G_Nbins ; }
    int    getN_srcJackBin    () const { return N_bins * ( I_Nbins+1 ) * M_Nbins * B_Nbins * G_Nbins * jacknifeBins ; }


    int    getUseNoise        () const { return  useNoise      ; }

    double getShapeNoise      () const { return  shapeNoise    ; }
    double getSigmaCrit       () const { return  sigmaC        ; }

    double getAlphaMin        () const { return       alphaMin ; }
    double getAlphaMax        () const { return       alphaMax ; }
    double getConMin          () const { return           cMin ; }
    double getConMax          () const { return           cMax ; }
    double getRMinFit         () const { return           rMin ; }
    double getRMaxFit         () const { return           rMax ; }
    double getMassMin         () const { return           mMin ; }
    double getMassMax         () const { return           mMax ; }
    double getChiMin          () const { return           cMin ; }
    double getChiMax          () const { return           cMax ; }

    double getTestVal         () const { return     avgTestVal ; }
    double getTolerance       () const { return      tolerance ; }
    double getMutChance       () const { return      mutChance ; }

    std::string getInputPath  () const { return  inputPath     ; }
    std::string getOutputPath () const { return  outputPath    ; }
    std::string getFoxH2012F  () const { return  fox2012F      ; }
    std::string getFoxH2123F  () const { return  fox2123F      ; }

    double getM_minBin () const { return M_minBin; }
    double getM_maxBin () const { return M_maxBin; }
    double getB_minBin () const { return B_minBin; }
    double getB_maxBin () const { return B_maxBin; }
    double getG_minBin () const { return G_minBin; }
    double getG_maxBin () const { return G_maxBin; }

    double getI_bin    ( int i ) const { return I_minBin + i *   I_step                          ; } // Returns the bin value at the index
    double getM_bin    ( int i ) const { return M_minBin + i * ( M_maxBin - M_minBin ) / M_Nbins ; } // Returns the bin value at the index
    double getB_bin    ( int i ) const { return B_minBin + i * ( B_maxBin - B_minBin ) / B_Nbins ; } // Returns the bin value at the index
    double getG_bin    ( int i ) const { return G_minBin + i * ( G_maxBin - G_minBin ) / G_Nbins ; } // Returns the bin value at the index


    double getN_haloBin ( int i, int m, int b, int g ) { return   i * G_Nbins * B_Nbins * M_Nbins
                                                                + m * G_Nbins * B_Nbins
                                                                + b * G_Nbins
                                                                + g ;                                      }

    double getSrcBin   (  int j, int i, int m, int b, int g, int bin ) { return   j * N_bins * G_Nbins * B_Nbins * M_Nbins * ( I_Nbins + 1 )
                                                                                + i * N_bins * G_Nbins * B_Nbins * M_Nbins
                                                                                + m * N_bins * G_Nbins * B_Nbins
                                                                                + b * N_bins * G_Nbins
                                                                                + g * N_bins
                                                                                + bin;                                      }


    double getSrcBin   ( int i, int m, int b, int g, int bin ) { return   i * N_bins * G_Nbins * B_Nbins * M_Nbins
                                                                        + m * N_bins * G_Nbins * B_Nbins
                                                                        + b * N_bins * G_Nbins
                                                                        + g * N_bins
                                                                        + bin;                                      }

    int    getSrcMBin  (  int j, int i, int b, int g, int bin ) { return   j * N_bins * G_Nbins * B_Nbins * ( I_Nbins + 1 )
                                                                         + i * N_bins * G_Nbins * B_Nbins
                                                                         + b * N_bins * G_Nbins
                                                                         + g * N_bins
                                                                         + bin;                                      }



  private:

    // Stuff to read in

    double shapeNoise ;

    int    N_bins     ;  // Number of bins for radial averaging
    int    num_threads;  // Number of threads for parallel processing

    std::string fox2012F;
    std::string fox2123F;
    std::string  inputPath;
    std::string outputPath;

    double sigmaC   ; // For all sources, same critical surface density

    double M_minBin ; // Min and max bin values for sources to be grouped into
    double M_maxBin ;
    double B_minBin ;
    double B_maxBin ;
    double G_minBin ;
    double G_maxBin ;
    double I_minBin ;
    double I_step   ; // No max value for I, just a step


    int    M_Nbins  ; // Number of bins for each
    int    B_Nbins  ;
    int    G_Nbins  ;
    int    I_Nbins  ;


    // Chi2 & genetic algorithm fitting boundaries
    double     cMin ;
    double     cMax ;
    double     mMin ;
    double     mMax ;
    double     rMin ;
    double     rMax ;
    double alphaMin ;
    double alphaMax ;

    int  useNoise      ; // Flag to add noise to simulation
    int jacknifeBins   ; // Number of bins used for jackknife errors
    int maxFitAttempts ; // Max attempts at fitting before abort
    int  N_chromosomes ; // Number of chromosomes in population
    int     N_chiTrack ; // Number of chi's to track for convergence
    int     consistent ; // Number of times need avg below tolerance
    double   tolerance ; // Average residual must be below tolerance
    double   mutChance ; // Likelihood of mutation
    double  avgTestVal ; // chiAvg*this is random range
    int  fitType       ; // 1 for 1 free param (in NFW), 3 for 3 free parameters

    int      firstFile ; // First file haloID to read in
    int       lastFile ;
};



// Stores values of FoxH functions used for generating
//  convergence of Einasto RTS
class einTable {

  public:

    einTable(){

    }


    void setBins  ( int a, int x ) { x_bins = x ;
                                     a_bins = a ;
                                     N_bins =a*x;
                                  initVals()    ; }

    void setX_min ( double     d ) { minX   = d ; }  // Boundaries of the table
    void setX_max ( double     d ) { maxX   = d ; }
    void setA_min ( double     d ) { minA   = d ; }
    void setA_max ( double     d ) { maxA   = d ; }


    double getX_min () { return   minX ; }
    double getX_max () { return   maxX ; }
    double getA_min () { return   minA ; }
    double getA_max () { return   maxA ; }

    int    getA_bins() { return a_bins ; }
    int    getX_bins() { return x_bins ; }
    int    getN_bins() { return N_bins ; }



    // Populate the array
    void setVal   ( int     aBin ,
                    int     xBin ,
                    double value ){

      if ( val != NULL ){

        val[ xBin + aBin * x_bins ] = value;

      }
    }


    double getVal ( int    aBin ,
                    int    xBin ){
      if ( val != NULL ){
        return val[ xBin + aBin * x_bins ];
      }
        return 0;
    }


  private:

    int    a_bins;
    int    x_bins;
    int    N_bins;
    double minA, minX;
    double maxA, maxX;

    double *val = NULL;

    void initVals () { val = new double[ N_bins ]; }

};


// Holds info on density profile, concentration, mass, shape param, etc
class densProfile{
  public:

    densProfile();
    densProfile( double inpA );

    // Modifiers, as parameters are modified need to adjust dependant values
    void setAlpha( double inpA ){         alpha = inpA; }
    void setR_max( double inpR ){         r_max = inpR; } // Mpc
    void setC    ( double inpC ){ concentration = inpC; }
    void setM_enc( double inpM ){         M_enc = inpM; }
    void setType ( int    inpT ){          type = inpT; }


    // If when getting an undefined variable, need to spit out a warning
    double getR_s  () const {  if ( concentration==-1.0 ) unDefVar("\"concentration\"");
                               if (         r_max==-1.0 ) unDefVar("\"R_max  \""      );   return r_max / concentration ;  }
    double getC    () const {  if ( concentration==-1.0 ) unDefVar("\"concentration\"");   return         concentration ;  }
    double getR_max() const {  if (         r_max==-1.0 ) unDefVar("\"R_max\""        );   return                 r_max ;  }
    double getM_enc() const {  if (         M_enc==-1.0 ) unDefVar("\"M_enc\""        );   return                 M_enc ;  }
    double getAlpha() const {  if (         alpha==-1.0 ) unDefVar("\"alpha\""        );   return                 alpha ;  }
    int    getType () const {                                                              return                  type ;  }

    double getRho_o() const ;

  private:
    double concentration;
    double alpha        ;
    double r_max        ;
    double M_enc        ;
    int    type         ; //1 NFW, 2 Einasto

    void   unDefVar( std::string inpS ) const {
      std::cerr<<"WARNING: variable "<<inpS<<" undefined in <lensProfile>"<<std::endl;
      exit(0);
    }
};


#endif
