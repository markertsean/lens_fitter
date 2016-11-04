#include "astro_constants.h"
#include "lensing_classes.h"
#include "pixelmap_functions.h"




// Will take and fill outArr with weighted average of a collapsed mass,
//   setting as our 1D arrays to pass to fitter
int *  avgBGArr  (  userInfo        u ,  // User info
                    double   * inpArr ,  // JIMBGR array to collapse into R
                    int      * inpN   ,  // Array containing number of sources to weight by
                    int             i ,  // Integration index
                    int             m ,  // Mass        index
                    int     omitIndex ,  // Index to omit for jack knifing
                    double  ** outArr )  // Array to allocate and populate with average
{

           * outArr = new double[ u.getNbins() ] (); // Allocate our output array
    int    * outN   = new int   [ u.getNbins() ] (); // Allocate our return array

    for ( int b = 0 ; b < u.getN_BBin() ; ++b ){     // We are collapsing B
    for ( int g = 0 ; g < u.getN_GBin() ; ++g ){     // We are collapsing G
    for ( int r = 0 ; r < u.getNbins () ; ++r ){
    for ( int j = 0 ; j < u.getN_JBin() ; ++j ){

        if ( j != omitIndex )                        // Will omit an index from jacknifing, to avoid use -1
        {
            int k = u.getSrcBin( j, i, m, b, g, r ) ;
            (*outArr)[ r ] += inpArr[k] * inpN[k] ;  // Weighted sum
              outN   [ r ] +=             inpN[k] ;

        }
    }
    }
    }
    }


    for ( int r = 0 ; r < u.getNbins() ; ++r )       // Normalize the sum
    {
        if ( outN[r] > 0 )
        (*outArr)[ r ] /= outN[ r ];
    }
    return outN;
}



// Will take and fill outArr with weighted average of a collapsed mass,
//   setting as our 1D arrays to pass to fitter
int *  avgMBArr  (  userInfo        u ,  // User info
                    double   * inpArr ,  // JIMBGR array to collapse into R
                    int      * inpN   ,  // Array containing number of sources to weight by
                    int             i ,  // Integration index
                    int             g ,  // gamma       index
                    int     omitIndex ,  // Index to omit for jack knifing
                    double  ** outArr )  // Array to allocate and populate with average
{

           * outArr = new double[ u.getNbins() ] (); // Allocate our output array
    int    * outN   = new int   [ u.getNbins() ] (); // Allocate our return array

    for ( int m = 0 ; m < u.getN_MBin() ; ++m ){     // We are collapsing M
    for ( int b = 0 ; b < u.getN_BBin() ; ++b ){     // Collapsing B
    for ( int r = 0 ; r < u.getNbins () ; ++r ){
    for ( int j = 0 ; j < u.getN_JBin() ; ++j ){

        if ( j != omitIndex )                        // Will omit an index from jacknifing, to avoid use -1
        {
            int k = u.getSrcBin( j, i, m, b, g, r ) ;
            (*outArr)[ r ] += inpArr[k] * inpN[k] ;  // Weighted sum
              outN   [ r ] +=             inpN[k] ;

        }
    }
    }
    }
    }


    for ( int r = 0 ; r < u.getNbins() ; ++r )       // Normalize the sum
    {
        if ( outN[r] > 0 )
        (*outArr)[ r ] /= outN[ r ];
    }
    return outN;
}



// Will take and fill outArr with weighted average of a collapsed mass,
//   setting as our 1D arrays to pass to fitter
int *   avgMArr  (  userInfo        u ,  // User info
                    double   * inpArr ,  // JIMBGR array to collapse into R
                    int      * inpN   ,  // Array containing number of sources to weight by
                    int             i ,  // Integration index
                    int             b ,  // b/a ratio   index
                    int             g ,  // gamma       index
                    int     omitIndex ,  // Index to omit for jack knifing
                    double  ** outArr )  // Array to allocate and populate with average
{

           * outArr = new double[ u.getNbins() ] (); // Allocate our output array
    int    * outN   = new int   [ u.getNbins() ] (); // Allocate our return array

    for ( int m = 0 ; m < u.getN_MBin() ; ++m ){     // We are collapsing M
    for ( int r = 0 ; r < u.getNbins () ; ++r ){
    for ( int j = 0 ; j < u.getN_JBin() ; ++j ){

        if ( j != omitIndex )                        // Will omit an index from jacknifing, to avoid use -1
        {
            int k = u.getSrcBin( j, i, m, b, g, r ) ;
            (*outArr)[ r ] += inpArr[k] * inpN[k] ;  // Weighted sum
              outN   [ r ] +=             inpN[k] ;

        }
    }
    }
    }


    for ( int r = 0 ; r < u.getNbins() ; ++r )       // Normalize the sum
    {
        if ( outN[r] > 0 )
        (*outArr)[ r ] /= outN[ r ];
    }
    return outN;
}



// Bon-Muller transformation to provide gaussian distribution
double gaussErr( double sigma ,
                 int     Ngal )
{

  float x1, x2, w;

  do {

    x1 = 2.0 * randVal( 0.0, 1.0 ) - 1.0;
    x2 = 2.0 * randVal( 0.0, 1.0 ) - 1.0;

    w  = x1 * x1 + x2 * x2;

  } while ( w >= 1.0 );

  w = std::sqrt( ( -2.0 * std::log( w ) ) / w );

  return x1 * w * sigma / std::sqrt( Ngal );
}



// Add gaussian error to arrays to fit
void  addgaussUncertaintyArr( double       *inpArr ,
                              double         sigma ,  // Amplitude of shape noise
                              int          * Nsrc  ,  // Array containing number of sources in each bin
                              int       N_elements )  // Number of elements in the array
{
    for ( int i = 0; i < N_elements; ++i )
        inpArr[i] = inpArr[i] + gaussErr( sigma, Nsrc[i] );
}


// Generates gaussian error to include in error arrays
double * gaussUncertaintyArr( double         sigma ,  // Amplitude of shape noise
                              int          * Nsrc  ,  // Array containing number of sources in each bin
                              int       N_elements )  // Number of elements in the array
{
    double * outArr = new double [ N_elements ] ();

    for ( int i = 0; i < N_elements; ++i )
        outArr[i] = std::sqrt( sigma*sigma / Nsrc[i]);

    return outArr;
}




// Collapses average array to match that of a collapsed M
haloInfo avgBGHaloInfo( userInfo    u ,
                        haloInfo  * h ,
                        int       * n ,
                        int         i ,
                        int         m )
{

    haloInfo outH;

    double M  = 0 ;
    double c  = 0 ;
    double ba = 0 ;
    double ca = 0 ;
    double rm = 0 ;
    double p  = 0 ;
    double t  = 0 ;
    double ga = 0 ;
    double a  = 0 ;

    int    N  = 0 ;


    for ( int b = 0; b < u.getN_BBin(); ++b ){
    for ( int g = 0; g < u.getN_GBin(); ++g ){

        int index = u.getN_haloBin( i, m, b, g );

        if ( n[ index ] > 0 )
        {

            M   = ( M  * N + h[ index ].getM     () * n[ index ] ) / ( N + n[ index ] ) ;
            c   = ( c  * N + h[ index ].getC     () * n[ index ] ) / ( N + n[ index ] ) ;
            ba  = ( ba * N + h[ index ].getBA    () * n[ index ] ) / ( N + n[ index ] ) ;
            ca  = ( ca * N + h[ index ].getCA    () * n[ index ] ) / ( N + n[ index ] ) ;
            rm  = ( rm * N + h[ index ].getRmax  () * n[ index ] ) / ( N + n[ index ] ) ;
            p   = ( p  * N + h[ index ].getPhi   () * n[ index ] ) / ( N + n[ index ] ) ;
            t   = ( t  * N + h[ index ].getTheta () * n[ index ] ) / ( N + n[ index ] ) ;
            ga  = ( ga * N + h[ index ].getGamma () * n[ index ] ) / ( N + n[ index ] ) ;
            a   = ( a  * N + h[ index ].getAlpha () * n[ index ] ) / ( N + n[ index ] ) ;

            N  +=            n[ index ] ;
        }
    }
    }

    outH.setM       (  M );
    outH.setC       (  c );
    outH.setBA      ( ba );
    outH.setCA      ( ca );
    outH.setRmax    ( rm );
    outH.setPhi     (  p );
    outH.setTheta   (  t );
    outH.setGamma   ( ga );
    outH.setAlpha   (  a );

    return outH;
}





// Collapses average array to match that of a collapsed M
haloInfo avgMBHaloInfo( userInfo    u ,
                        haloInfo  * h ,
                        int       * n ,
                        int         i ,
                        int         g )
{

    haloInfo outH;

    double M  = 0 ;
    double c  = 0 ;
    double ba = 0 ;
    double ca = 0 ;
    double rm = 0 ;
    double p  = 0 ;
    double t  = 0 ;
    double ga = 0 ;
    double a  = 0 ;

    int    N  = 0 ;

    for ( int m = 0; m < u.getN_MBin(); ++m ){
    for ( int b = 0; b < u.getN_BBin(); ++b ){

        int index = u.getN_haloBin( i, m, b, g );

        if ( n[ index ] > 0 )
        {

            M   = ( M  * N + h[ index ].getM     () * n[ index ] ) / ( N + n[ index ] ) ;
            c   = ( c  * N + h[ index ].getC     () * n[ index ] ) / ( N + n[ index ] ) ;
            ba  = ( ba * N + h[ index ].getBA    () * n[ index ] ) / ( N + n[ index ] ) ;
            ca  = ( ca * N + h[ index ].getCA    () * n[ index ] ) / ( N + n[ index ] ) ;
            rm  = ( rm * N + h[ index ].getRmax  () * n[ index ] ) / ( N + n[ index ] ) ;
            p   = ( p  * N + h[ index ].getPhi   () * n[ index ] ) / ( N + n[ index ] ) ;
            t   = ( t  * N + h[ index ].getTheta () * n[ index ] ) / ( N + n[ index ] ) ;
            ga  = ( ga * N + h[ index ].getGamma () * n[ index ] ) / ( N + n[ index ] ) ;
            a   = ( a  * N + h[ index ].getAlpha () * n[ index ] ) / ( N + n[ index ] ) ;

            N  +=            n[ index ] ;
        }
    }
    }

    outH.setM       (  M );
    outH.setC       (  c );
    outH.setBA      ( ba );
    outH.setCA      ( ca );
    outH.setRmax    ( rm );
    outH.setPhi     (  p );
    outH.setTheta   (  t );
    outH.setGamma   ( ga );
    outH.setAlpha   (  a );

    return outH;
}






// Collapses average array to match that of a collapsed M
haloInfo  avgMHaloInfo( userInfo    u ,
                        haloInfo  * h ,
                        int       * n ,
                        int         i ,
                        int         b ,
                        int         g )
{

    haloInfo outH;

    double M  = 0 ;
    double c  = 0 ;
    double ba = 0 ;
    double ca = 0 ;
    double rm = 0 ;
    double p  = 0 ;
    double t  = 0 ;
    double ga = 0 ;
    double a  = 0 ;

    int    N  = 0 ;

    for ( int m = 0; m < u.getN_MBin(); ++m ){

        int index = u.getN_haloBin( i, m, b, g );

        if ( n[ index ] > 0 )
        {

            M   = ( M  * N + h[ index ].getM     () * n[ index ] ) / ( N + n[ index ] ) ;
            c   = ( c  * N + h[ index ].getC     () * n[ index ] ) / ( N + n[ index ] ) ;
            ba  = ( ba * N + h[ index ].getBA    () * n[ index ] ) / ( N + n[ index ] ) ;
            ca  = ( ca * N + h[ index ].getCA    () * n[ index ] ) / ( N + n[ index ] ) ;
            rm  = ( rm * N + h[ index ].getRmax  () * n[ index ] ) / ( N + n[ index ] ) ;
            p   = ( p  * N + h[ index ].getPhi   () * n[ index ] ) / ( N + n[ index ] ) ;
            t   = ( t  * N + h[ index ].getTheta () * n[ index ] ) / ( N + n[ index ] ) ;
            ga  = ( ga * N + h[ index ].getGamma () * n[ index ] ) / ( N + n[ index ] ) ;
            a   = ( a  * N + h[ index ].getAlpha () * n[ index ] ) / ( N + n[ index ] ) ;

            N  +=            n[ index ] ;
        }
    }

    outH.setM       (  M );
    outH.setC       (  c );
    outH.setBA      ( ba );
    outH.setCA      ( ca );
    outH.setRmax    ( rm );
    outH.setPhi     (  p );
    outH.setTheta   (  t );
    outH.setGamma   ( ga );
    outH.setAlpha   (  a );

    return outH;
}




// Calculates jackknife errors for each profile set
void jacknife( densProfile  *profile ,
               int         N_samples ,
               double       * errArr )
{

  // Variances

  double M_var(0);
  double C_var(0);
  double A_var(0);
  double R_var(0);

  double Mbias(0);
  double Cbias(0);
  double Abias(0);
  double Rbias(0);

  double Mavg(0);
  double Cavg(0);
  double Aavg(0);
  double Ravg(0);

  double Mbar_i[ N_samples ];
  double Cbar_i[ N_samples ];
  double Abar_i[ N_samples ];
  double Rbar_i[ N_samples ];


  // Generate xbar_i
  for ( int i = 0; i < N_samples; ++i ){
      Abar_i[ i ]  = 0;
      Mbar_i[ i ]  = 0;
      Cbar_i[ i ]  = 0;
      Rbar_i[ i ]  = 0;
  for ( int j = 0; j < N_samples; ++j ){

    if ( j != i ){
      if ( profile[j+1].getType() == 2 )
      Abar_i[ i ] +=             profile[j+1].getAlpha()  ;
      Mbar_i[ i ] += std::log10( profile[j+1].getM_enc() );
      Cbar_i[ i ] +=             profile[j+1].getC    ()  ;
      Rbar_i[ i ] +=             profile[j+1].getR_max()  ;
    }

  }
      Abar_i[ i ]  = Abar_i[ i ] / ( N_samples - 1 );
      Cbar_i[ i ]  = Cbar_i[ i ] / ( N_samples - 1 );
      Mbar_i[ i ]  = Mbar_i[ i ] / ( N_samples - 1 );
      Rbar_i[ i ]  = Rbar_i[ i ] / ( N_samples - 1 );

      if ( profile[i+1].getType() == 2 )
      Aavg        +=             profile[i+1].getAlpha()   / N_samples ;
      Mavg        += std::log10( profile[i+1].getM_enc() ) / N_samples ;
      Cavg        +=             profile[i+1].getC    ()   / N_samples ;
      Ravg        +=             profile[i+1].getR_max()   / N_samples ;

  }


  double coeff = ( N_samples - 1.0 ) / N_samples;

  // Generate variance
  for ( int i = 0; i < N_samples; ++i ){
    A_var += ( Abar_i[i] - Aavg ) * ( Abar_i[i] - Aavg ) ;
    C_var += ( Cbar_i[i] - Cavg ) * ( Cbar_i[i] - Cavg ) ;
    M_var += ( Mbar_i[i] - Mavg ) * ( Mbar_i[i] - Mavg ) ;
    R_var += ( Rbar_i[i] - Ravg ) * ( Rbar_i[i] - Ravg ) ;
  }

  A_var = A_var * ( N_samples - 1.0 ) / N_samples ;
  C_var = C_var * ( N_samples - 1.0 ) / N_samples ;
  M_var = M_var * ( N_samples - 1.0 ) / N_samples ;
  R_var = R_var * ( N_samples - 1.0 ) / N_samples ;

  errArr[0] = C_var ;
  errArr[1] = M_var ;
  errArr[2] = A_var ;
  errArr[3] = R_var ;

}
