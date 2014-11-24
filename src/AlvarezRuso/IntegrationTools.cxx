//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steve Boyd ( s.b.boyd \at warwick.ac.uk)
   University of Warwick

*/
//____________________________________________________________________________

#include <cstdlib>
#include <complex>

#include <TMath.h>

#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "AlvarezRuso/IntegrationTools.h"

namespace genie {
namespace alvarezruso {

typedef std::complex<double> cdouble;

//
// Routine to work out the evaluation points for a Gaussian integral given the
// end points, and the number of sampling points
//
void integrationtools::SG20R(const double a, const double b, const unsigned int n, const unsigned int nsamp, 
           double* x, unsigned int& np, double* w)
{
  static const double y[10] = {.9931285991, .9639719272, .9122344282, .8391169718, .7463319064, .6360536807, .5108670019, .3737060887, .2277858511, .0765265211 };
  np = 20 * n;
  double dint = (b - a) / double(n);
  double delt = dint * 0.5;
  double orig = a - delt;
  int i1 = -nsamp;
  int i2, j1, j2;
  double dorig;
  for(unsigned int i = 1; i <= n; i++)
  {
    orig += dint;
    dorig = orig + orig;
    i1 += nsamp;
    i2 = i1 + nsamp+1;
    for(unsigned int j = 1; j <= 10; j++)
    {
      j1 = i1 + j;
      j2 = i2 - j;
      x[j1-1] = orig - delt * y[j-1];
      x[j2-1] = dorig - x[j1-1];
    }
  }
}

//-----------------------------------------------------------------------------------------------------------
// Gaussian-Legendre integration of the function defined by CF
cdouble integrationtools::RG201D(const double A, const double B, const unsigned int N, const unsigned int nsamp, const cdouble CF[])
{
  // Gaussian-Legendre integration of the function defined by CF
  const double W[10] = {.0176140071,.0406014298,.0626720483,.0832767415,.1019301198,.1181945319,.1316886384,.1420961093,.1491729864,.1527533871};
  cdouble CR(0.0, 0.0);
  int I1 = -nsamp;
  int I2, J1, J2;
  for(unsigned int i = 1; i <= N; ++i)
  {
    I1 += nsamp;
    I2 = I1 + nsamp+1;
    for(unsigned int j = 1; j <= 10; ++j)
    {
      J1 = I1 + j;
      J2 = I2 - j;
      CR += W[j-1] * (CF[J1-1]+CF[J2-1]);
    }
  } 
  //CRES=CR*0.5*(B-A)/float(N)
  cdouble CRES = CR*0.5*(B-A)/Double_t(N);
  return CRES;
}
//-----------------------------------------------------------------------------------------------------------
// Gaussian-Legendre integration of the function defined by CF
void integrationtools::RG202D(const double a, const double b, unsigned int n, unsigned int l, 
              unsigned int m, std::vector< std::vector<cdouble> >& cf, 
              const unsigned int nsamp, std::vector<cdouble>& cres)
{
  // This is a fast integrator based on a Gauss-Legendre method. This only support two-dimensional integration
  n = 2; l = 0; m = 3;

  static const double w[10] = {.0176140071,.0406014298,.0626720483,.0832767415,.1019301198,.1181945319,.1316886384,.1420961093,.1491729864,.1527533871};

  std::vector<cdouble> cr(4, cdouble(0.0,0.0));

  int i1 = -nsamp;
  int i2;
  int j1;
  int j2;

  for(unsigned int i = 0; i != n; ++i)
  {
    i1 += nsamp; 
    i2 = i1 + nsamp-1;
    for(unsigned int j = 0; j != 10; ++j)
    {
      j1 = i1 + j;
      j2 = i2 - j;

      for(unsigned int ll = l; ll <= m; ++ll)
      {
        cr[ll] += w[j] * ( cf[ll][j1] + cf[ll][j2] );
      }
    }
  }

  for(unsigned int i = 0; i != 4; ++i)
  {
    cres[i] = cr[i] * 0.5 * (b-a) / static_cast<double>(n);
  }
}

} // alvarezruso namespace
} // genie namespace
