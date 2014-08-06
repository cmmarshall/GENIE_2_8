//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Chris Marshall <marshall \at pas.rochester.edu>
          University of Rochester
          Martti Nirkko
          University of Berne

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TMath.h>
#include <TF3.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Conventions/Units.h"
#include "AtharSingleKaon/ASKKinematicsGenerator.h"
#include "Conventions/KinePhaseSpace.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGCore/EventGeneratorI.h"
#include "EVGCore/RunningThreadInfo.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepFlags.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
ASKKinematicsGenerator::ASKKinematicsGenerator() :
KineGeneratorWithCache("genie::ASKKinematicsGenerator")
{
  //fEnvelope = 0;
}
//___________________________________________________________________________
ASKKinematicsGenerator::ASKKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::ASKKinematicsGenerator", config)
{
  //fEnvelope = 0;
}
//___________________________________________________________________________
ASKKinematicsGenerator::~ASKKinematicsGenerator()
{
  //if(fEnvelope) delete fEnvelope;
}
//___________________________________________________________________________
void ASKKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  if(fGenerateUniformly) {
    LOG("ASKKinematics", pNOTICE)
          << "Generating kinematics uniformly over the allowed phase space";
  }

  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();
  CalculateKin_AtharSingleKaon(evrec);
}
//___________________________________________________________________________
void ASKKinematicsGenerator::CalculateKin_AtharSingleKaon(GHepRecord * evrec) const
{
  // Get the Primary Interacton object
  Interaction * interaction = evrec->Summary();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  // Initialise a random number generator 
  RandomGen * rnd = RandomGen::Instance();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
  double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);

  // maximum kinetic energy for kaon or muon is just Enu - mk - ml
  int leppdg = interaction->FSPrimLeptonPdg();
  const TLorentzVector P4_nu = *(interaction->InitStatePtr()->GetProbeP4(kRfLab));
  double enu = interaction->InitState().ProbeE(kRfLab); // Enu in lab frame
  double mk = (interaction->ExclTag().StrangeHadronPdg()==kPdgK0) ? 0.497614 : 0.493677; // kaon mass
  double ml = 0.0;
  if( leppdg == kPdgMuon ) ml = kMuonMass;
  else if( leppdg == kPdgElectron ) ml = kElectronMass;
  else if( leppdg == kPdgTau ) ml = kTauMass;

  const double Tkmax = enu - mk - ml;
  const double Tlmax = enu - mk - ml;

  // Tkmax must be > 0, otherwise we're below threshold
  if( Tkmax <= 0.0 ) {
    LOG("ASKKinematics", pWARN) << "No available phase space";
    evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
    genie::exceptions::EVGThreadException exception;
    exception.SetReason("No available phase space");
    exception.SwitchOnFastForward();
    throw exception;
  }

  // For speed we will draw from x = log(1-costheta)

  const double Tkmin = kASmallNum;
  const double Tlmin = kASmallNum;
  const double xmin = -20.0; // this is artificial cut-off and should be a parameter
  const double xmax =  0.69314718056; // log(2) is physical boundary
  const double phikqmin = 0.0;
  const double phikqmax = 2.0 * kPi;
  const double dtk = Tkmax - Tkmin;
  const double dtl = Tlmax - Tlmin;
  const double dx = xmax - xmin;
  const double dphikq = phikqmax - phikqmin;  

  // for uniform generation only
  const double cosThetalMin = -1.0;
  const double cosThetalMax = 1.0;
  const double dcosthetal = cosThetalMax - cosThetalMin;

  //------ Try to select a valid tk, tl, costhetal, phikq quadruplet

  unsigned int iter = 0;
  bool accept=false;
  double xsec=-1, tk = -1, tl = -1, costhetal = -1, phikq = -1;

  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
        LOG("ASKKinematics", pWARN)
             << "*** Could not select a valid (tk, tl, costhetal) triplet after "
                                               << iter << " iterations";
        evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
     }

     if(fGenerateUniformly) {
        //-- Generate a x,y pair uniformly in the kinematically allowed range.
        tk = Tkmin + dtk * rnd->RndKine().Rndm();
        tl = Tlmin + dtl * rnd->RndKine().Rndm();
        costhetal = cosThetalMin + dcosthetal * rnd->RndKine().Rndm();
        phikq = phikqmin + dphikq * rnd->RndKine().Rndm();
     } else {
        tk = Tkmin + dtk * rnd->RndKine().Rndm();
        tl = Tlmin + dtl * rnd->RndKine().Rndm();
        double x  = xmin + dx * rnd->RndKine().Rndm(); // log(1-costheta)
        phikq = phikqmin + dphikq * rnd->RndKine().Rndm();
        //costhetal = cosThetalMin + dcosthetal * rnd->RndKine().Rndm();

        costhetal = 1.0 - TMath::Exp(x);

     }

     LOG("ASKKinematics", pDEBUG) << "Trying: Tk = " << tk << ", Tl = " << tl << ", cosThetal = " << costhetal << ", phikq = " << phikq;

     interaction->KinePtr()->SetKV(kKVTk, tk);
     interaction->KinePtr()->SetKV(kKVTl, tl);
     interaction->KinePtr()->SetKV(kKVctl, costhetal);
     interaction->KinePtr()->SetKV(kKVphikq, phikq);

     // Set Q2 and W
     double el = tl + ml;
     double pl = TMath::Sqrt(el*el - ml*ml);
     double Mf = (interaction->ExclTag().NProtons()==1) ? kProtonMass : kNeutronMass; 
     TVector3 lepton_3vector = TVector3(0,0,0);
     lepton_3vector.SetMagThetaPhi(pl,TMath::ACos(costhetal),0.0); //phi_l doesn't affect q2, we'll choose it later
     TLorentzVector P4_lep    = TLorentzVector(lepton_3vector , el );
     TLorentzVector q = P4_nu - P4_lep;
     double Q2 = -q.Mag2();
     double nu = P4_nu.E() - P4_lep.E(); // hadronic system energy nu

     interaction->KinePtr()->SetQ2(Q2);
     double W2 = Mf*Mf + 2*Mf*nu - Q2;

     // computing cross section for the current kinematics
     xsec = fXSecModel->XSec(interaction, kPSTkTlctl);

     //-- decide whether to accept the current kinematics
     if(!fGenerateUniformly) {
        // for now, just use a full box, up to the max cross section
        // if its slow then later we might use an envelope
        double max = xsec_max; //fEnvelope->Eval(tk,tl,costhetal);
        double t   = max * rnd->RndKine().Rndm();
        double J   = TMath::Abs(1. - costhetal);

        this->AssertXSecLimits(interaction, xsec*J, max);

        LOG("ASKKinematics", pDEBUG) 
            << "xsec= " << xsec << ", J= " << J << ", xsec*J = " << xsec*J << " Rnd= " << t;

        if( xsec*J > xsec_max ) {
          LOG("ASKKinematics", pWARN)
             << "!!!!!!!!!!!!! xsec= " << xsec << ", J= " << J << ", xsec*J = " << xsec*J << " max= " << xsec_max;
        }

        accept = (t< J*xsec);
     }
     else { 
        accept = (xsec>0);
     }

     //-- If the generated kinematics are accepted, finish-up module's job
     if(accept) {

        // calculate the stuff

        // for uniform kinematics, compute an event weight as
        // wght = (phase space volume)*(differential xsec)/(event total xsec)
        if(fGenerateUniformly) {
          double wght = 1.0; // change this
          wght *= evrec->Weight();
          LOG("ASKKinematics", pNOTICE) << "Current event wght = " << wght;
          evrec->SetWeight(wght);
        }

        // reset bits
        interaction->ResetBit(kISkipProcessChk);
        interaction->ResetBit(kISkipKinematicChk);

        interaction->KinePtr()->SetKV(kKVSelTk, tk);
        interaction->KinePtr()->SetKV(kKVSelTl, tl);
        interaction->KinePtr()->SetKV(kKVSelctl, costhetal);
        interaction->KinePtr()->SetKV(kKVSelphikq, phikq);
        interaction->KinePtr()->SetQ2(Q2, true);
        interaction->KinePtr()->SetW(TMath::Sqrt(W2), true);
        interaction->KinePtr()->Setx( Q2/(2*Mf*nu), true );
        interaction->KinePtr()->Sety( nu/P4_nu.E(), true );
        interaction->KinePtr()->ClearRunningValues();

        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec*TMath::Abs(1. - costhetal),kPSTkTlctl);

        return;
     }
  }// iterations
}
//___________________________________________________________________________
double ASKKinematicsGenerator::ComputeMaxXSec(const Interaction * in) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  SLOG("ASKKinematics", pDEBUG)
          << "Scanning the allowed phase space {K} for the max(dxsec/d{K})";
#endif

  double max_xsec = 0;
  //double Ev = in->InitState().ProbeE(kRfLab);

  const int Ntk = 50;
  const int Ntl = 50;
  const int Nctl = 50;
  //const int Nphikq = 10;

  int leppdg = in->FSPrimLeptonPdg();
  double enu = in->InitState().ProbeE(kRfLab); // Enu in lab frame
  double mk = (in->ExclTag().StrangeHadronPdg()==kPdgK0) ? 0.497614 : 0.493677; // kaon mass
  double ml = 0.0;
  if( leppdg == kPdgMuon ) ml = kMuonMass;
  else if( leppdg == kPdgElectron ) ml = kElectronMass;
  else if( leppdg == kPdgTau ) ml = kTauMass;

  const double Tkmax = enu - mk - ml;
  const double Tlmax = enu - mk - ml;
  const double Tkmin = kASmallNum;
  const double Tlmin = kASmallNum;
  const double xmin = -20.0;
  const double xmax =  0.69314718056;
  //const double phikqmin = kASmallNum;
  //const double phikqmax = 2. * kPi;
  const double dtk = (Tkmax - Tkmin)/Ntk;
  const double dtl = (Tlmax - Tlmin)/Ntl;
  const double dx = (xmax - xmin)/Nctl;
  //const double dphikq = (phikqmax-phikqmin)/Nphikq;

  for(int i=0; i<Ntk; i++) {
    double tk = Tkmin + dtk*i;
    for(int j=0; j<Ntl; j++) {
      double tl = Tlmin + dtl*j;
      for(int k=0; k<Nctl; k++) {
        double logoneminuscosthetal = xmin + dx*k;
        double ctl = 1.0 - TMath::Exp(logoneminuscosthetal);
        //for(int m=0; m<Nphikq; m++) {
        // in lieu of loop over phi_kq, just set equal to pi, where differential XSec is always maximized
        // physically, this means the phase space is better for hadron to be opposite lepton in lab frame?
          double phikq = kPi; //phikqmin + dphikq*m;

          in->KinePtr()->SetKV(kKVTk, tk);
          in->KinePtr()->SetKV(kKVTl, tl);
          in->KinePtr()->SetKV(kKVctl, ctl);
          in->KinePtr()->SetKV(kKVphikq, phikq);

          double xsec = fXSecModel->XSec(in, kPSTkTlctl);
          // want maximum of d4sigma / dtk dtl dphikq dx where x is log(1 - costheta)
          // that is this xsec * dcostheta/dx = 1-costheta
          xsec *= (1.0 - ctl);

          LOG("ASKKinematics", pDEBUG)  
         	 << "Max diff XS serach: xsec(tk= " << tk << ", tl= " << tl << ", ctl= " << ctl << ", phiqk= " << phikq << ") = " << xsec;

          max_xsec = TMath::Max(max_xsec, xsec);
        //}//phiqk
      }//ctl
    }//tl
  }//tk

  LOG("ASKKinmatics", pINFO) << "Max XSec is " << max_xsec;

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy.
  max_xsec *= fSafetyFactor;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  SLOG("ASKKinematics", pDEBUG) << in->AsString();
  SLOG("ASKKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("ASKKinematics", pDEBUG) << "Computed using alg = " << fXSecModel->Id();
#endif



  return max_xsec;
}
/*
//___________________________________________________________________________
double ASKKinematicsGenerator::MaxXSec(const GHepRecord * ev) const
{
  Interaction * in = ev->Summary();
  return this->MaxXSec(in);
}

//___________________________________________________________________________
double ASKKinematicsGenerator::MaxXSec(const Interaction * in) const
{
  double max_xsec = 0;
  //double Ev = in->InitState().ProbeE(kRfLab);

  const int Ntk = 50;
  const int Ntl = 50;
  const int Nctl = 50;
  const int Nphikq = 10;

  int leppdg = in->FSPrimLeptonPdg();
  double enu = in->InitState().ProbeE(kRfLab); // Enu in lab frame
  double mk = (in->ExclTag().StrangeHadronPdg()==kPdgK0) ? 0.497614 : 0.493677; // kaon mass
  double ml = 0.0;
  if( leppdg == kPdgMuon ) ml = kMuonMass;
  else if( leppdg == kPdgElectron ) ml = kElectronMass;
  else if( leppdg == kPdgTau ) ml = kTauMass;

  const double Tkmax = enu - mk - ml;
  const double Tlmax = enu - mk - ml;
  const double Tkmin = kASmallNum;
  const double Tlmin = kASmallNum;
  const double xmin = -10.0;
  const double xmax =  0.69314718056;
  const double phikqmin = kASmallNum;
  const double phikqmax = 2. * kPi;
  const double dtk = (Tkmax - Tkmin)/Ntk;
  const double dtl = (Tlmax - Tlmin)/Ntl;
  const double dx = (xmax - xmin)/Nctl;
  const double dphikq = (phikqmax-phikqmin)/Nphikq;

  for(int i=0; i<Ntk; i++) {
    double tk = Tkmin + dtk*i;
    for(int j=0; j<Ntl; j++) {
      double tl = Tlmin + dtl*j;
      for(int k=0; k<Nctl; k++) {
        double logoneminuscosthetal = xmin + dx*k;
        double ctl = 1.0 - TMath::Exp(logoneminuscosthetal);
        for(int m=0; m<Nphikq; m++) {
          double phikq = phikqmin + dphikq*m;

          in->KinePtr()->SetKV(kKVTk, tk);
          in->KinePtr()->SetKV(kKVTl, tl);
          in->KinePtr()->SetKV(kKVctl, ctl);
          in->KinePtr()->SetKV(kKVphikq, phikq);

          double xsec = fXSecModel->XSec(in, kPSTkTlctl);

          LOG("ASKKinematics", pDEBUG)  
         	 << "Max diff XS serach: xsec(tk= " << tk << ", tl= " << tl << ", ctl= " << ctl << ", phiqk= " << phikq << ") = " << xsec;

          max_xsec = TMath::Max(max_xsec, xsec);
        }//phiqk
      }//ctl
    }//tl
  }//tk

  LOG("ASKKinmatics", pINFO) << "Max XSec is " << max_xsec;

  //this->CacheMaxXSec(in, max_xsec);  

  return max_xsec;
}

*/
//___________________________________________________________________________
double ASKKinematicsGenerator::Energy(const Interaction * interaction) const
{
// Override the base class Energy() method to cache the max xsec for the
// neutrino energy in the LAB rather than in the hit nucleon rest frame.

  const InitialState & init_state = interaction->InitState();
  double E = init_state.ProbeE(kRfLab);
  return E;
}
//___________________________________________________________________________
void ASKKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ASKKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ASKKinematicsGenerator::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  //-- max xsec safety factor (for rejection method) and min cached energy
  fSafetyFactor = fConfig->GetDoubleDef("MaxXSec-SafetyFactor", 1.5);
  fEMin         = fConfig->GetDoubleDef("Cache-MinEnergy",      0.5);

  //-- Generate kinematics uniformly over allowed phase space and compute
  //   an event weight?
  fGenerateUniformly = fConfig->GetBoolDef("UniformOverPhaseSpace", false);

  //-- Maximum allowed fractional cross section deviation from maxim cross
  //   section used in rejection method
  fMaxXSecDiffTolerance = 
         fConfig->GetDoubleDef("MaxXSec-DiffTolerance",999999.);
  assert(fMaxXSecDiffTolerance>=0);

  //-- Envelope employed when    sampling is used 
  //   (initialize with dummy range)
  //   not using this right now
//  if(fEnvelope) delete fEnvelope;
//  fEnvelope = new TF3("envelope",
//    	  kinematics::COHImportanceSamplingEnvelope,0.,1,0.,1,2);
}
//____________________________________________________________________________

