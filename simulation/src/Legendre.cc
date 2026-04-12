/*
  Legendre.cc
*/

#include "Legendre.hh"

bool PLGen( int n, double x, double *PL )
{
  if(n<0) return false;
  PL[0]=1.;
  if(n==0) return true;
  PL[1]=x;
  if(n==1) return true;

  for( int i=2; i<=n; ++i ){
    PL[i]=(PL[i-1]*x*double(2*i-1)-PL[i-2]*double(i-1))/double(i);
  }
  return true;
}
