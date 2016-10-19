#include "astro_constants.h"
#include "lensing_classes.h"
#include "pixelmap_functions.h"

// Collapses axis into M, using weighted sum
void    collapseM(  userInfo       u ,
                    double  *    arr ,
                    int     *      N ,
                    double  ** m_Arr )
{

    int      n_m   = u.getN_srcJackBin() / u.getN_MBin();
           * m_Arr = new double [ n_m ];
    int    * n_Arr = new int    [ n_m ];

    // Does weighted sum
    for ( int j = 0; j < u.getN_JBin(); ++j ) {
    for ( int i = 0; j < u.getN_IBin(); ++i ) {
    for ( int m = 0; j < u.getN_MBin(); ++m ) {
    for ( int b = 0; j < u.getN_BBin(); ++b ) {
    for ( int g = 0; j < u.getN_GBin(); ++g ) {
    for ( int r = 0; j < u.getNbins (); ++r ) {

        int     k    = u.getSrcBin  ( j, i, m, b, g, r );
        int    kk    = u.getSrcMBin ( j, i,    b, g, r );

        (*m_Arr)[ kk ] += arr[k] * N[k] ;
          n_Arr [ kk ] +=          N[k] ;

    }
    }
    }
    }
    }
    }

    // Does division of weighted sum
    for ( int j = 0; j < u.getN_JBin(); ++j ) {
    for ( int i = 0; j < u.getN_IBin(); ++i ) {
    for ( int b = 0; j < u.getN_BBin(); ++b ) {
    for ( int g = 0; j < u.getN_GBin(); ++g ) {
    for ( int r = 0; j < u.getNbins (); ++r ) {

        int    kk    = u.getSrcMBin ( j, i,    b, g, r );

        (*m_Arr)[ kk ] /= n_Arr[ kk ] ;

    }
    }
    }
    }
    }

}


// Collapses axis into M, using weighted sum
void    collapseM(  userInfo       u ,
                    double  *    arr ,
                    int     *      N ,
                    double  ** m_Arr ,
                    int     ** n_Arr )
{

    int      n_m   = u.getN_srcJackBin() / u.getN_MBin();
           * m_Arr = new double [ n_m ];
           * n_Arr = new int    [ n_m ];

    // Does weighted sum
    for ( int j = 0; j < u.getN_JBin(); ++j ) {
    for ( int i = 0; j < u.getN_IBin(); ++i ) {
    for ( int m = 0; j < u.getN_MBin(); ++m ) {
    for ( int b = 0; j < u.getN_BBin(); ++b ) {
    for ( int g = 0; j < u.getN_GBin(); ++g ) {
    for ( int r = 0; j < u.getNbins (); ++r ) {

        int        k    = u.getSrcBin  ( j, i, m, b, g, r );
        int       kk    = u.getSrcMBin ( j, i,    b, g, r );

        (*m_Arr)[ kk ] += arr[k] * N[k] ;
        (*n_Arr)[ kk ] +=          N[k] ;

    }
    }
    }
    }
    }
    }

    // Does division of weighted sum
    for ( int j = 0; j < u.getN_JBin(); ++j ) {
    for ( int i = 0; j < u.getN_IBin(); ++i ) {
    for ( int b = 0; j < u.getN_BBin(); ++b ) {
    for ( int g = 0; j < u.getN_GBin(); ++g ) {
    for ( int r = 0; j < u.getNbins (); ++r ) {

        int       kk    = u.getSrcMBin ( j, i,    b, g, r );

        (*m_Arr)[ kk ] /= (*n_Arr)[ kk ] ;

    }
    }
    }
    }
    }

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




// Calculates jackknife errors for each profile set
void jacknife( densProfile  *profile ,
               int         N_samples ,
               double       * errArr )
{

  // Variances

  double M_var(0);
  double C_var(0);
  double A_var(0);

  double Mbias(0);
  double Cbias(0);
  double Abias(0);

  double Mavg(0);
  double Cavg(0);
  double Aavg(0);

  double Mbar_i[ N_samples ];
  double Cbar_i[ N_samples ];
  double Abar_i[ N_samples ];


  // Generate xbar_i
  for ( int i = 0; i < N_samples; ++i ){
      Abar_i[ i ]  = 0;
      Mbar_i[ i ]  = 0;
      Cbar_i[ i ]  = 0;
  for ( int j = 0; j < N_samples; ++j ){

    if ( j != i ){
      if ( profile[j+1].getType() == 2 )
      Abar_i[ i ] +=             profile[j+1].getAlpha()  ;
      Mbar_i[ i ] += std::log10( profile[j+1].getM_enc() );
      Cbar_i[ i ] +=             profile[j+1].getC    ()  ;
    }

  }
      Abar_i[ i ]  = Abar_i[ i ] / ( N_samples - 1 );
      Cbar_i[ i ]  = Cbar_i[ i ] / ( N_samples - 1 );
      Mbar_i[ i ]  = Mbar_i[ i ] / ( N_samples - 1 );

      if ( profile[i+1].getType() == 2 )
      Aavg        +=             profile[i+1].getAlpha()   / N_samples ;
      Mavg        += std::log10( profile[i+1].getM_enc() ) / N_samples ;
      Cavg        +=             profile[i+1].getC    ()   / N_samples ;

  }


  double coeff = ( N_samples - 1.0 ) / N_samples;

  // Generate variance
  for ( int i = 0; i < N_samples; ++i ){
    A_var += ( Abar_i[i] - Aavg ) * ( Abar_i[i] - Aavg ) ;
    C_var += ( Cbar_i[i] - Cavg ) * ( Cbar_i[i] - Cavg ) ;
    M_var += ( Mbar_i[i] - Mavg ) * ( Mbar_i[i] - Mavg ) ;
  }

  A_var = A_var * ( N_samples - 1.0 ) / N_samples ;
  C_var = C_var * ( N_samples - 1.0 ) / N_samples ;
  M_var = M_var * ( N_samples - 1.0 ) / N_samples ;

  errArr[0] = C_var;
  errArr[1] = M_var;
  errArr[2] = A_var;

}
