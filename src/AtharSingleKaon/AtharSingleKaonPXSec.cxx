//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Chris Marshall and Martti Nirkko

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Base/XSecIntegratorI.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/KineVar.h"
#include "AtharSingleKaon/AtharSingleKaonPXSec.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/NuclearUtils.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

//____________________________________________________________________________
AtharSingleKaonPXSec::AtharSingleKaonPXSec() :
XSecAlgorithmI("genie::AtharSingleKaonPXSec")
{

}
//____________________________________________________________________________
AtharSingleKaonPXSec::AtharSingleKaonPXSec(string config) :
XSecAlgorithmI("genie::AtharSingleKaonPXSec", config)
{

}
//____________________________________________________________________________
AtharSingleKaonPXSec::~AtharSingleKaonPXSec()
{

}
//____________________________________________________________________________
double AtharSingleKaonPXSec::XSec(
                  const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //const InitialState & init_state = interaction -> InitState();
  //const Kinematics &   kinematics = interaction -> Kine();
  //const Target &       target     = init_state.Tgt();

  //double E  = init_state.ProbeE(kRfLab);
  //double Q2 = kinematics.Q2();
  //int    Z  = target.Z(); 
  //int    N  = target.N();

  double xsec  = 1.0;

  //-- The algorithm computes dxsec/d?
  //   Check whether variable tranformation is needed
  //if(kps!=kPSQ2fE) { // change this
  //  double J = utils::kinematics::Jacobian(interaction,kPSQ2fE,kps); //change this
  //  xsec *= J;
  //}

  return xsec;
}
//____________________________________________________________________________
double AtharSingleKaonPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool AtharSingleKaonPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info = interaction->ProcInfo();  
  if(!proc_info.IsDeepInelastic() || !proc_info.IsWeakCC()) return false;

  return true;
}
//____________________________________________________________________________
void AtharSingleKaonPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AtharSingleKaonPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AtharSingleKaonPXSec::LoadConfig(void)
{
  //AlgConfigPool * confp = AlgConfigPool::Instance();
  //const Registry * gc = confp->GlobalParameterList();

  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________
