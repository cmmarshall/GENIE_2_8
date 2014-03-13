//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors:

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TVector3.h>

#include "Conventions/Constants.h"
#include "AtharSingleKaon/ASKHadronicSystemGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/PrintUtils.h"
#include "Base/XSecAlgorithmI.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGCore/EventGeneratorI.h"
#include "EVGCore/RunningThreadInfo.h"


using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
ASKHadronicSystemGenerator::ASKHadronicSystemGenerator() :
HadronicSystemGenerator("genie::ASKHadronicSystemGenerator")
{

}
//___________________________________________________________________________
ASKHadronicSystemGenerator::ASKHadronicSystemGenerator(string config) :
HadronicSystemGenerator("genie::ASKHadronicSystemGenerator", config)
{

}
//___________________________________________________________________________
ASKHadronicSystemGenerator::~ASKHadronicSystemGenerator()
{

}
//___________________________________________________________________________
void ASKHadronicSystemGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// Access cross section algorithm for running thread
  //RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  //const EventGeneratorI * evg = rtinfo->RunningThread();
  //const XSecAlgorithmI *fXSecModel = evg->CrossSectionAlg();
  CalculateHadronicSystem_AtharSingleKaon(evrec);
}
//___________________________________________________________________________
void ASKHadronicSystemGenerator::CalculateHadronicSystem_AtharSingleKaon(GHepRecord * evrec) const
{
//
// This method generates the final state hadronic system (kaon + nucleus) in 
// ASK interactions
//
  RandomGen * rnd = RandomGen::Instance();

  Interaction * interaction = evrec->Summary();
  Kinematics * kinematics = interaction->KinePtr();
  const XclsTag & xcls_tag  = interaction->ExclTag();

  //-- Access neutrino, initial nucleus and final state prim. lepton entries
  GHepParticle * nu  = evrec->Probe();
  GHepParticle * Ni  = evrec->HitNucleon();
  GHepParticle * fsl = evrec->FinalStatePrimaryLepton();
  assert(nu);
  assert(Ni);
  assert(fsl);

  const TLorentzVector & vtx   = *(nu->X4());
  const TLorentzVector & p4nu  = *(nu ->P4());
  const TLorentzVector & p4fsl = *(fsl->P4());

  //-- Determine the pdg code of the final state nucleon
  int nuc_pdgc = Ni->Pdg(); // same as the initial nucleus
  int kaon_pdgc = xcls_tag.StrangeHadronPdg();

  const Target & tgt = interaction->InitState().Tgt();
  GHepStatus_t ist = (tgt.IsNucleus()) ?
                          kIStHadronInTheNucleus : kIStStableFinalState;

  //-- basic kinematic inputs
  double E    = nu->E();  
  double M    = (xcls_tag.NProtons()==1) ? kProtonMass : kNeutronMass;
  double mk   = PDGLibrary::Instance()->Find(kaon_pdgc)->Mass();
  double mk2  = TMath::Power(mk,2);

  //-- specific kinematic quantities
  double kaon_T = kinematics->GetKV(kKVSelTk);
  double kaon_E = kaon_T + mk;
  double pk = sqrt( kaon_E*kaon_E - mk2 );

  TLorentzVector q = p4nu - p4fsl;
  TVector3 qvec = q.Vect();
  double Q2 = q.Mag2();
  double q3 = qvec.Mag();
  // Get pN from energy conservation
  // q.E() + M = sqrt(pk*pk + mk2) + sqrt(pN*pN + M*M)
  double eN = q.E() + M - sqrt(pk*pk + mk2);
  double pN = sqrt( eN*eN - M*M );

  // now we get the angles using:
  // pk*sin(thetaK) = pN*sin(thetaN) and q3=pk*cos(thetaK) + pN*cos(thetaN)
  double pkl = (pk*pk - pN*pN + q3*q3) / (2*q3);
  double pkt = sqrt(pk*pk - pkl*pkl);

  double pNt = pkt;
  double pNl = sqrt(pN*pN - pNt*pNt);

  // These are longitudinal and transverse w.r.t. q vector
  TVector3 kaon = TVector3(0., pkt, pkl);
  TVector3 nucleon = TVector3(0., pNt, pNl);

  // Pick a phi for the kq plane relative to nul plane (random)
  double phi_kq = 2. * kPi * rnd->RndKine().Rndm();

  // kaon and nucleon are in plane with opposite phis
  kaon.RotateZ(phi_kq);
  nucleon.RotateZ(-phi_kq);

  // Now we need to rotate to the q vector direction
  kaon.RotateUz(qvec.Unit());
  nucleon.RotateUz(qvec.Unit());

  double pxNf = nucleon.Px();
  double pyNf = nucleon.Py();
  double pzNf = nucleon.Pz();
  double ENf = sqrt(pxNf*pxNf + pyNf*pyNf + pzNf*pzNf + M*M);  

  double pxKf = kaon.Px();
  double pyKf = kaon.Py();
  double pzKf = kaon.Pz();
  double EKf = sqrt(pxKf*pxKf + pyKf*pyKf + pzKf*pzKf + mk2);  

  //-- Save the particles at the GHEP record

  // mom is mother, not momentum
  int mom = evrec->TargetNucleusPosition();

  // AddParticle (int pdg, GHepStatus_t ist, int mom1, int mom2, int dau1, int dau2, double px, double py, double pz, double E, double x, double y, double z, double t)
  
  evrec->AddParticle(
     nuc_pdgc, ist, mom,-1,-1,-1, 
     pxNf, pyNf, pzNf, ENf, vtx.X(), vtx.Y(), vtx.Z(), vtx.T());

  evrec->AddParticle(
     kaon_pdgc,ist, mom,-1,-1,-1, 
     pxKf, pyKf, pzKf, EKf, vtx.X(), vtx.Y(), vtx.Z(), vtx.T());
}

