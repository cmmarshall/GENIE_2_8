//____________________________________________________________________________
/*!

\class    genie::alvarezruso::AlvarezRusoCOHPiPDXsec

\brief    5d differential cross section for Alvarez-Ruso Coherent Pion Production xsec

\ref      

\author   Steve Dennis
          University of Warwick, Rutherford Appleton Laboratory

\created  05/12/2013

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________
#ifndef _AR_COH_MULTIDIFF_H_
#define _AR_COH_MULTIDIFF_H_

#include <TMath.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <Math/LorentzVector.h>

#include "AlvarezRuso/ARSampledNucleus.h"
#include "AlvarezRuso/ARConstants.h"
#include "AlvarezRuso/ARWavefunction.h"
#include "Utils/NuclearUtils.h"

#include <complex>

namespace genie
{
namespace alvarezruso
{

class ARWFSolution;

enum current_t{kCC, kNC};
enum flavour_t{kE, kMu, kTau};
enum nutype_t{kNu,kAntiNu};
enum formfactors_t{kNieves, kGarcia};

class AlvarezRusoCOHPiPDXSec
{
  public:
    
    AlvarezRusoCOHPiPDXSec(int Z_, int A_, const current_t current_, 
          const flavour_t flavour_ = kE, const nutype_t nutype = kNu, 
          const formfactors_t ff_ = kNieves);
    ~AlvarezRusoCOHPiPDXSec();
    double DXSec(const double E_nu_, const double E_l_, const double theta_l_, 
           const double phi_l_, const double theta_pi_, const double phi_pi_);
    void SetDebug(bool debug)  {  debug_ = debug;  }
    
    // Fill the ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >s based on the values from the kinematics
    void SetKinematics();
    
    // Fill values based on the flavour
    void SetFlavour();
    
    // Fill values based on the current
    void SetCurrent();

    int GetSampling() {
      return sampling;
     }

    std::complex<double> DeltaCouplingInMed(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > delta_momentum, 
         ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > pion_momentum, double density_p, double density_n);
    double PiDecayVertex(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > pion_momentum, double mass);
    std::complex<double>  DeltaPropagatorInMed(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > delta_momentum);
    double DeltaWidthPauliBlocked(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > delta_momentum, double density);
    double DeltaWidthFree(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > delta_momentum);
    std::complex<double>  H(unsigned int i, unsigned int j) const;
    double DifferentialCrossSection();
    double PionMomentumCM(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > delta_momentum);
    double PNVertexFactor(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > momentum, double mass);
    double DeltaSelfEnergyRe(double density);
    double DeltaSelfEnergyIm(double density);
    double DeltaSelfEnergyConstant(double a, double b, double c, double E);
    std::complex<double> NucleonPropagator(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > nucleon_momentum);
    
    void NuclearCurrent(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > q, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > pdir, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > pcrs, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > ppi, std::complex<double>  *jPtr);

    // Fill the wavefunctions
    void SolveWavefunctions();
    
    //______________________________________________________________
    // Properties
    
    bool debug_;
    int fZ;
    int fA;
    unsigned int sampling;
    int n_elements;
    // Choice of current
    current_t current;
    // Choice of neutrino flavour
    flavour_t flavour;
    // Chocie of initial neutrino type
    nutype_t nutype;
    // Choice of form-factor approximation
    formfactors_t formfactors;
    // Nuclear values
    ARSampledNucleus nucleus;
    // Wavefunction calculator
    ARWFSolution* wfsolution;
    
    // Kinematics of the event
    double E_nu;      // initial neutrino energy [GeV]
    double E_l;      // scattered lepton energy [GeV]
    double theta_l;      // scattered lepton angle
    double theta_pi;    // pion angle
    double phi;      // angle between lepton and pion
    
    double fLastE_pi;
    
    // Four-momenta of particles and transfers involved
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > q;    // momentum-transfer
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p_nu;    // incoming neutrino
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p_l;    // outgoing lepton/neutrino
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p_pi;    // outgoing pion
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p_n_i;    // incoming (stationary) nucleon
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p_n_o;    // outgoing nucleon
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p_direct;  // intermediary particle (Delta/nucleon) in direct diagrams
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p_cross;    // intermediary particle (Delta/nucleon) in crossed diagrams
    
    ARConstants constant;
    
    // Values and constants which will be used during the running
    // Get set to one of the mass values above in the constructor
    double m_pi;      // mass of the pion
    double m_l;      // mass of the lepton
    double g_factor;    // Current factor. if(NC) = 1/2 G^2. if(CC) = 1/2 G^2 cos^2(theta_c)
    
    // Form factors for PiNDelta decay vertex
    double F_direct_delta;
    double F_direct_nucleon;
    double F_cross_delta;
    double F_cross_nucleon;
    
    // Wavefunction
    ARWavefunction* uwave;
    ARWavefunction* uwaveDr;
    ARWavefunction* uwaveDtheta;
    
    std::complex<double>  J_hadronic[4];
};

} //namespace alvarezruso
} //namespace genie
#endif
