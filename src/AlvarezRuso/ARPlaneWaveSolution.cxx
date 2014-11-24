//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Daniel Scully ( d.i.scully \at warwick.ac.uk)
   University of Warwick

*/
//____________________________________________________________________________

#include <iostream>
#include <sstream>
#include <string>
#include <complex>

#include "AlvarezRuso/ARPlaneWaveSolution.h"
#include "AlvarezRuso/AlvarezRusoCOHPiPDXSec.h"
#include "AlvarezRuso/ARSampledNucleus.h"
#include "AlvarezRuso/ARConstants.h"
#include "Messenger/Messenger.h"

namespace genie
{
namespace alvarezruso
{
 
ARPlanewaveSolution::ARPlanewaveSolution(bool debug): ARWFSolution(debug)
{
  if(debug_) std::cerr << "PwS:WFS@ constructor" << std::endl;
}

ARPlanewaveSolution::~ARPlanewaveSolution()
{
}

std::complex<double> ARPlanewaveSolution::Element(const double radius, const double cosine_rz, const double p_pion,const pUnits_t piUnit)
{
  if( radius == 0.0 )    LOG("AlvarezRusoXSecPi",pWARN) << "PlanewaveSolution/Element >> Warning: radius = " << radius   ;
  if( cosine_rz == 0.0 ) LOG("AlvarezRusoXSecPi",pWARN) << "PlanewaveSolution/Element >> Warning: cosine = " << cosine_rz;
  if( p_pion == 0.0 )    LOG("AlvarezRusoXSecPi",pWARN) << "PlanewaveSolution/Element >> Warning: p_pion = " << p_pion   ;
  std::complex<double> result = exp( std::complex<double>(0,1) * p_pion * radius * cosine_rz);
  if( real(result) == 0.0 ) LOG("AlvarezRusoXSecPi",pWARN) << "PlanewaveSolution/Element >> Warning: result = " << result;
  return result;
}

void ARPlanewaveSolution::Solve()
{
  if(0) std::cout << "Hi!" << std::endl;
}

} //namespace alvarezruso
} //namespace genie

