//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Daniel Scully ( d.i.scully \at warwick.ac.uk)
   University of Warwick

*/
//____________________________________________________________________________

#include "ARSampledNucleus.h"
#include "Messenger/Messenger.h"
#include "Utils/StringUtils.h"

#include <string>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <istream>
#include <complex>

#include <TSystem.h>
#include <TF1.h>
#include <TMath.h>
#include <Math/GaussLegendreIntegrator.h>

#include "Utils/NuclearUtils.h"
#include "AlvarezRuso/IntegrationTools.h"
#include "AlvarezRuso/ARConstants.h"

//
// Equation/Table numbers refer to:
// J. Nieves and E. Oset
// A Theoretical approach to pionic atoms and the problem of anomalies
// Nuclear Physics A 554 (1993) 509-553
//

namespace genie {
namespace alvarezruso {

ARSampledNucleus::ARSampledNucleus()
{
  densities        = NULL;
  radii            = NULL;
  sample_points_1  = NULL;
  sample_points_2  = NULL;
  sample_weights_1 = NULL;
  sample_weights_2 = NULL;
}

ARSampledNucleus::ARSampledNucleus(const unsigned int& ZNumber, const unsigned int& ANumber, bool debug) : debug_(debug)
{
  sampling = 20;
  n_densities = 2*sampling;
  a = ANumber;
  z = ZNumber;
  
  densities = NULL;
  radii = NULL;
  sample_points_1  = NULL;
  sample_points_2  = NULL;
  sample_weights_1 = NULL;
  sample_weights_2 = NULL;
  
  this->Fill();
}

ARSampledNucleus::~ARSampledNucleus()
{
  for(unsigned int i = 0; i != n_densities; ++i)
  {
    if (densities && densities[i] ) delete[] densities[i];
    if (radii     && radii    [i] ) delete[] radii    [i];
  }
  
  if (densities       ) delete[] densities;
  if (radii           ) delete[] radii    ;
  if (sample_points_1 ) delete[] sample_points_1;
  if (sample_points_2 ) delete[] sample_points_2;
  if (sample_weights_1) delete[] sample_weights_1;
  if (sample_weights_2) delete[] sample_weights_2;
}

void ARSampledNucleus::SetSampling(const unsigned int & nsamp)
{
  sampling = nsamp;
  n_densities = 2*sampling;
  this->FillSamplePoints();
}
void ARSampledNucleus::Fill()
{
  // Values from Table 1 of Nieves et al.
  
  r_max = 3.0 * TMath::Power(this->A(), (1.0/3.0));
  
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("AR_PiWFunction_Table", pDEBUG)<< "N:: r_max = " << r_max
  << "N:: z = " << z 
  << "N:: a = " << a 
#endif

  this->FillSamplePoints();
  this->FillDensities();
}

double ARSampledNucleus::Density(const int i, const int j) const
{
  return densities[i][j];
}

double ARSampledNucleus::Radius(const int i, const int j) const
{
  return radii[i][j];
}

void ARSampledNucleus::FillSamplePoints()
{
  if (sample_points_1)  delete[] sample_points_1;
  if (sample_points_2)  delete[] sample_points_2;
  if (sample_weights_1) delete[] sample_weights_1;
  if (sample_weights_2) delete[] sample_weights_2;
  
  sample_points_1 = new double[n_densities];
  sample_points_2 = new double[n_densities];
  
  sample_weights_1 = new double[n_densities];
  sample_weights_2 = new double[n_densities];
  unsigned int decoy;
  integrationtools::SG20R(0.0, r_max, 2, sampling, sample_points_1, decoy, sample_weights_1);
  integrationtools::SG20R(-r_max, r_max, 2, sampling, sample_points_2, decoy, sample_weights_2);
}

void ARSampledNucleus::FillDensities()
{
  double r;
  
  for(unsigned int i = 0; i != n_densities; ++i)
  {
    if (densities && densities[i]) delete[] densities[i];
    if (radii     && radii    [i]) delete[] radii    [i];
  }
  if (densities) delete[] densities;
  if (radii    ) delete[] radii    ;
  
  densities = new double*[n_densities];
  radii     = new double*[n_densities];
  
  for(unsigned int i = 0; i != n_densities; ++i)
  {
    densities[i] = new double[n_densities];
    radii    [i] = new double[n_densities];
    
    for(unsigned int j = 0; j != n_densities; ++j)
    {
      r = TMath::Sqrt( sample_points_1[i]*sample_points_1[i] + sample_points_2[j]*sample_points_2[j] );
      radii[i][j] = r;
      densities[i][j]=utils::nuclear::Density(r, a);
    }
  }
}

} //namespace alvarezruso
} //namespace genie
