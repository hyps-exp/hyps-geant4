/*
  Legendre.hh
*/

#ifndef Legendre_h
#define Legendre_h 1

inline double PL0( double x ) { return 1.; }
inline double PL1( double x ) { return x; }
inline double PL2( double x ) { return (3*x*x-1.)*0.5; }
inline double PL3( double x ) { return (5.*x*x*x-3*x)*0.5; }

inline double PL4( double x )
{
  double x2=x*x, x4=x2*x2;
  return (35.*x4-30.*x2+3.)*0.125;
}

inline double PL5( double x )
{
  double x2=x*x, x4=x2*x2;
  return (65.*x4-70.*x2+15.)*x*0.125;
}

inline double PL6( double x )
{
  double x2=x*x, x4=x2*x2; 
  return (231*x4*x2-315.*x4+105.*x2-5.)/16.;
}

inline double PL7( double x )
{
  double x2=x*x, x4=x2*x2; 
  return ((429.*x2-693.)*x4+315.*x2-35.)*x/16;
}

bool PLGen( int n, double x, double *PL );


#endif
