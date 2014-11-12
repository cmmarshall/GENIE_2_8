//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Daniel Scully ( d.i.scully \at warwick.ac.uk)
   University of Warwick
 Based on Fortran code by Luis Alvarez-Ruso.
*/
//____________________________________________________________________________

// C++
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <complex>

// Root
#include <TVector3.h>
#include <TMath.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <Math/LorentzVector.h>

//Genie
#include "Algorithm/AlgConfigPool.h"
#include "Messenger/Messenger.h"
#include "Conventions/Constants.h"

//AlvarezRuso
#include "AlvarezRuso/AlvarezRusoCOHPiPDXSec.h"
#include "AlvarezRuso/ARSampledNucleus.h"
#include "AlvarezRuso/AREikonalSolution.h"
#include "AlvarezRuso/IntegrationTools.h"
#include "AlvarezRuso/ARWavefunction.h"

using namespace genie::constants;

typedef std::complex<double> cdouble;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef ROOT::Math::SVector< cdouble , 4> CVector;

namespace genie {
namespace alvarezruso {

AlvarezRusoCOHPiPDXSec::AlvarezRusoCOHPiPDXSec(int Z_, int A_, const current_t current_, 
   const flavour_t flavour_, const nutype_t nutype_,const formfactors_t ff_)
  : debug_(false),
  fZ(Z_),
  fA(A_),
  sampling(20),
  n_elements(sampling * 2),
  current( current_ ),
  flavour( flavour_ ),
  nutype( nutype_ ),
  formfactors( ff_ ),
  nucleus( fZ, fA, false ),
  wfsolution( new AREikonalSolution(debug_, this) ),
  fLastE_pi(-1.),
  uwave( new ARWavefunction(sampling, debug_) ),
  uwaveDr( new ARWavefunction(sampling, debug_) ),
  uwaveDtheta( new ARWavefunction(sampling, debug_) )
{
  SetCurrent();
  SetFlavour();
  nucleus.SetSampling(sampling);
  
}

AlvarezRusoCOHPiPDXSec::~AlvarezRusoCOHPiPDXSec()
{
  delete this->wfsolution;
  delete this->uwave;
  delete this->uwaveDr;
  delete this->uwaveDtheta;
}


double AlvarezRusoCOHPiPDXSec::DXSec( const double E_nu_, const double E_l_, const double theta_l_, const double phi_l_, const double theta_pi_, const double phi_pi_)
{

  E_nu = E_nu_;
  E_l = E_l_;
  theta_l = theta_l_;
  theta_pi = theta_pi_;
  phi = phi_pi_ - phi_l_;
  
  if( (E_nu/constant.HBar()) < (m_pi + m_l) )
  {
    return 0.0;
  }
  else if( (E_l / constant.HBar()) < m_l )
  {
    return 0.0;
  }
  
  SetKinematics();
  
  if( p_pi.E() < m_pi )
  {
    return 0.0;
  }
  
  if( TMath::Abs((q - p_pi).M()) > (2.0 / constant.HBar()) )
  {
    // Comment from original fortran:
    //    This is to eliminate very high momentum transfers to the nucleus, which
    //    have negligible impact on observables but might create numerical instabilities
    return 0.0;
  }

  F_direct_delta = PiDecayVertex(p_pi, constant.DeltaPMass());
  F_cross_delta = PiDecayVertex(-p_pi, constant.DeltaPMass());
  F_direct_nucleon = PiDecayVertex(p_pi, constant.NucleonMass());
  F_cross_nucleon = PiDecayVertex(-p_pi, constant.NucleonMass());
  
  
  // Only need to resolve wave funtions if Epi changes
  if ( (TMath::Abs(fLastE_pi-p_pi.E()) > 1E-10) ){ 
    SolveWavefunctions();
  }

  LorentzVector pni = p_pi - q;
  pni *= 0.5;
  pni.SetE( constant.NucleonMass() );
  p_direct = q + pni;
  p_cross = pni - p_pi;
  NuclearCurrent(q, p_direct, p_cross, p_pi, J_hadronic);

  double dxsec = DifferentialCrossSection();
  
  fLastE_pi = p_pi.E();

  return dxsec;
}



cdouble AlvarezRusoCOHPiPDXSec::H(unsigned int i, unsigned int j) const
{
  cdouble H = ( conj(J_hadronic[i]) * J_hadronic[j] );
  return H;
}


double AlvarezRusoCOHPiPDXSec::DifferentialCrossSection()
{
  const cdouble i(0,1);

  cdouble term1,term2,term3,term4;
  if (nutype == kNu) {
     term1 = q.X() * ( H(0,1) - i*H(0,2) + H(1,0) - H(1,3) + i*H(2,0) - i*H(2,3) -H(3,1) + i*H(3,2) );
     term2 = q.Z() * ( -H(0,0) + H(0,3) + H(1,1) - i*H(1,2) + i*H(2,1) + H(2,2) + H(3,0) - H(3,3) );
     term3 = 2.0 * p_nu.E() * ( H(0,0) - H(0,3) - H(3,0) + H(3,3) );
     term4 = -q.E() * ( H(0,0) - H(0,3) + H(1,1) - i*H(1,2) + i*H(2,1) + H(2,2) - H(3,0) + H(3,3) );
  } else {
     term1 = q.X() * ( H(0,1) + i*H(0,2) + H(1,0) - H(1,3) - i*H(2,0) + i*H(2,3) -H(3,1) - i*H(3,2) );
     term2 = q.Z() * ( -H(0,0) + H(0,3) + H(1,1) + i*H(1,2) - i*H(2,1) + H(2,2) + H(3,0) - H(3,3) );
     term3 = 2.0 * p_nu.E() * ( H(0,0) - H(0,3) - H(3,0) + H(3,3) );
     term4 = -q.E() * ( H(0,0) - H(0,3) + H(1,1) + i*H(1,2) - i*H(2,1) + H(2,2) - H(3,0) + H(3,3) );
  }
  
  cdouble complex_amp2 = 8.0 * p_nu.E() * (term1 + term2 + term3 + term4);
  
  double amp2 = real(complex_amp2);
  
  double d5 = g_factor / 8.0 * p_l.P() * p_pi.P() / p_nu.E() / (32 * kPi4 * kPi ) * 
                           amp2 * constant.cm38Conversion() / constant.HBar();
  
  return d5;
}


/* *********************************************************************
 * Vertex factor for Delta/Nucleon decaying to Pi + Nucleon
 * formerly offshell()
 */
double AlvarezRusoCOHPiPDXSec::PiDecayVertex(LorentzVector momentum, double mass)
{  
  // s-channel form factor for the piNDelta vertex as in M. Post thesis (pg 35)
  // Calculated for a NUCLEON AT REST
  
  double Lam = 1.0 / constant.HBar();
  double Lam4 = Lam*Lam*Lam*Lam;
  
  double mass2 = mass*mass;
  
  LorentzVector decaying_momentum(momentum.x(),momentum.y(),momentum.z(), (momentum.E()+constant.NucleonMass()));
  
  double factor = decaying_momentum.mag2() - mass2;
  factor *= factor;

  double ofshel = Lam4 / ( Lam4 + factor );
  return ofshel;
}

/* *********************************************************************
 * Fill the LorentzVectors with values based on the values from the
 * kinematics
 */
void AlvarezRusoCOHPiPDXSec::SetKinematics()
{
  // Initial neutrino momentum
  p_nu = LorentzVector(0, 0, (E_nu/constant.HBar()), (E_nu/constant.HBar()) );
  
  // Final lepton momentum
  double mod_p_l = TMath::Sqrt( (E_l/constant.HBar())*(E_l/constant.HBar()) - m_l*m_l );
  double p_l_x = mod_p_l * TMath::Sin(theta_l);
  double p_l_z = mod_p_l * TMath::Cos(theta_l);
  p_l = LorentzVector(p_l_x, 0, p_l_z, (E_l/constant.HBar()));
  
  // Momentum transfer
  q = p_nu - p_l;
  
  // Pion momentum
  double E_pi = q.E();
  double mod_p_pi = TMath::Sqrt( E_pi*E_pi - m_pi*m_pi );
  double p_pi_x = mod_p_pi * TMath::Sin(theta_pi) * TMath::Cos(phi);
  double p_pi_y = mod_p_pi * TMath::Sin(theta_pi) * TMath::Sin(phi);
  double p_pi_z = mod_p_pi * TMath::Cos(theta_pi);
  p_pi = LorentzVector(p_pi_x, p_pi_y, p_pi_z, E_pi);
}


/* *********************************************************************
 * 
 */
void AlvarezRusoCOHPiPDXSec::SetFlavour()
{
  if(current == kNC)
  {
    m_l = 0.0;
  }
  else if(current == kCC)
  {
    switch(flavour)
    {
      case kE:
        m_l = constant.ElectronMass();
        break;
      case kMu:
        m_l = constant.MuonMass();
        break;
      case kTau:
        m_l = constant.TauMass();
        break;
      default:
        LOG("AlvarezRusoCohPi",pERROR) << "[ERROR] Unknown lepton flavour";
        exit(1);
    }
  }
  else
  {
    LOG("AlvarezRusoCohPi",pERROR) << "[ERROR] Unknown current"; 
    exit(1);
  }
}





/* *********************************************************************
 * Fill values based on the current
 */
void AlvarezRusoCOHPiPDXSec::SetCurrent()
{
  switch(this->current)
  {
    case kCC:
      m_pi = constant.PiPMass();
      g_factor = 0.5 * constant.GFermi()*constant.GFermi()*constant.CosCabibboAngle()*constant.CosCabibboAngle();
      break;
    case kNC:
      m_pi = constant.Pi0Mass();
      g_factor = 0.5 * constant.GFermi()*constant.GFermi();
      break;
    default:
      LOG("AlvarezRusoCohPi",pERROR) << "[ERROR] Unknown current type";
      exit(1);
  }
}




/*
 * Solve the wavefunctions
 */

/// This is only a function of the nucleus and pion momentum/energy
/// so if neither of those have changed there is no need to re-calculate
/// the wavefunction values.
/// Such a caching has not been implemented here yet!

void AlvarezRusoCOHPiPDXSec::SolveWavefunctions()
{
  unsigned int n_points = n_elements;
  
  double x1;
  double x2;
  double radius;
  double cosine_rz; // angle w.r.t the pion momentum
  
  // for calculating derivatives
  double delta;
  cdouble uwave_plus;
  cdouble uwave_minus;

  // Loop over grid of points in the nuclear potential
  for(unsigned int i = 0; i != n_points; ++i)
  {
    for(unsigned int j = 0; j != n_points; ++j)
    {    
      x1 = nucleus.SamplePoint1(i);
      x2 = nucleus.SamplePoint2(j);

      // radius of position in potential from centre
      radius = nucleus.Radius(i,j);
      // angle of sampling point wrt to neutrino direction
      cosine_rz = x2 / radius;
    
      // Calculate wavefunction
      uwave->set(i, j, wfsolution->Element(radius, -cosine_rz, 
                          p_pi.P(), kInNatural ));
      delta = 0.0001;
      if( radius < delta ) delta = radius;
    
      // Calculate derivative of wavefunction in the radial direction
      uwave_plus = wfsolution->Element( (radius+delta), -cosine_rz, 
                                     p_pi.P(), kInNatural);
      uwave_minus = wfsolution->Element( (radius-delta), -cosine_rz, 
                                     p_pi.P(), kInNatural);
      uwaveDr->set(i, j, ((uwave_plus - uwave_minus) / (2.0 * delta)) );
      
      // Calculate derivative of wavefunction in the angle space
      delta = 0.0001;
      if( (cosine_rz - delta) < -1.0 )  delta = cosine_rz + 1.0;
      else if( (cosine_rz + delta) > 1.0 )  delta = 1.0 - cosine_rz;
      uwave_plus = wfsolution->Element(radius, -(cosine_rz+delta), 
                                      p_pi.P(), kInNatural);
      uwave_minus = wfsolution->Element(radius, -(cosine_rz-delta), 
                                      p_pi.P(), kInNatural);
      uwaveDtheta->set( i, j, ((uwave_plus - uwave_minus) / (2.0 * delta)) );
    }
  }
  
}

cdouble AlvarezRusoCOHPiPDXSec::DeltaPropagatorInMed(LorentzVector delta_momentum)
{
  //Energy dependent in-medium Delta propagator
  double W = TMath::Abs( delta_momentum.mag() );
  double width = DeltaWidthPauliBlocked(delta_momentum, 0.0);  
  double imSigma = DeltaSelfEnergyIm(0.0);  
  double reSigma = DeltaSelfEnergyRe(0.0);
  
  cdouble denom1( (W + constant.DeltaPMass() + reSigma), 0.0);
  cdouble denom2( (W - constant.DeltaPMass() - reSigma), ((width/2.0) - imSigma) );
  
  cdouble result = 1.0 / (denom1 * denom2);
  return result;
}

double AlvarezRusoCOHPiPDXSec::DeltaWidthPauliBlocked(LorentzVector delta_momentum, double density)
{
   //In-medium Delta width including Pauli blocking  
  
  double width;
  double free_width = DeltaWidthFree(delta_momentum);
  
  if(free_width == 0.0)
  {
    width = 0.0;
  }
  else
  {
    double f = 0.0;
    double p_f = TMath::Power( ((3.0/2.0)*constants::kPi2*density), (1.0/3.0) );  // Fermi-momentum
    double p_cm = PionMomentumCM(delta_momentum);  // nucleon (and pion) momentum in CoM
    
    if(formfactors == kNieves)
    {
      // Use the approximation from Nieves et al. NPA 554(93)554
      double r = p_cm / p_f;
      if(r > 1.0)
      {
        double r2 = r*r;
        f = 1.0 + ( (-2.0)/(5.0*r2) + (9.0)/(35.0*r2*r2) - (2.0)/(21.0*r2*r2*r2) );
      }
      else
      {
        f = (34.0/35.0)*r - (22.0/105.0)*r*r*r;
      }
    }
    else if(formfactors == kGarcia)
    {
      //Use the approximation from Garcia-Recio, NPA 526(91)685
      // Delta inv. mass
      double wd = TMath::Abs( delta_momentum.M());
      // Fermi-energy
      double E_f = TMath::Sqrt( constant.NucleonMassSq() + p_f*p_f );
       // Delta 3-momentum in Lab.
      double kd = delta_momentum.R();
      // Nucleon energy in CoM
      double E_n = TMath::Sqrt( constant.NucleonMassSq() + p_cm*p_cm );
      f = (kd*p_cm + delta_momentum.E()*E_n - E_f*wd) / (2.0*kd*p_cm);
      
      if(f < 0.0)      f = 0;
      else if(f > 1.0)  f = 1.0;
    }
    else
    {
      //LOG("AlvarezRusoCohPi",pERROR) << "[ERROR] Choice of form-factor approximation not properly made";
      exit(1);
    }
    
    width = free_width*f;
  }
  return width;
}

double AlvarezRusoCOHPiPDXSec::DeltaWidthFree(LorentzVector delta_momentum)
{  
   // Free Delta width  
  double pre_factor_1 = 1.0 / (6.0 * kPi );
  double pre_factor_2 = constant.DeltaNCoupling_over_MassPiP_Squared();
  double qcm = PionMomentumCM(delta_momentum);
  double qcm3 = qcm*qcm*qcm;
  // Luis' code has the next pre-factor equal to
  // double pre_factor_3 = constant.NucleonMass() / TMath::Sqrt(TMath::Abs( delta_momentum.mag2() ));
  // but the paper uses the Delta mass?
  double pre_factor_3 = constant.NucleonMass() / TMath::Sqrt(TMath::Abs( delta_momentum.mag2() ));
  
  double delta_width = pre_factor_1 * pre_factor_2 * pre_factor_3 * qcm3;
  
  return delta_width;
}

// 1.1.1.1
/*
 * Calculate the three-momentum of a pion coming from the decay of a
 * Delta resonance.
 */
double AlvarezRusoCOHPiPDXSec::PionMomentumCM(LorentzVector delta_momentum) // qcm
{
  double m_pi2 = m_pi*m_pi;
  double m_n2 = constant.NucleonMassSq();
  double s = TMath::Abs( delta_momentum.mag2() );
  
  double p_pi_cm = ((s-m_pi2-m_n2)*(s-m_pi2-m_n2)) - 4.0*m_pi2*m_n2;
  
  if(p_pi_cm <= 0.0)
  {
  //  LOG("AlvarezRusoCohPi",pERROR) << "[ERROR] Delta's s is insufficient for decay";
    return 0.0;
  }
  
  p_pi_cm = TMath::Sqrt(p_pi_cm/ s) / 2.0;
  return p_pi_cm;
}

double AlvarezRusoCOHPiPDXSec::PNVertexFactor(LorentzVector momentum, double mass)
{
  // s-channel form factor for the piNDelta/piNN vertex as in M. Post thesis (pg 35)
  // Calculated for a NUCLEON AT REST
  
  double lambda = 1.0 / constant.HBar();
  double lambda4 = lambda*lambda*lambda*lambda;
  
  double mass2 = mass*mass;
  
  LorentzVector decaying_momentum(momentum.x(), momentum.y(),momentum.z(), (momentum.E()+constant.NucleonMass()));
  
  double factor = decaying_momentum.mag2() - mass2;
  factor *= factor;
  
  double result = lambda4 / (lambda4 + factor);
  
  return result;
}

double AlvarezRusoCOHPiPDXSec::DeltaSelfEnergyRe(double density)
{
  double result = (0.04/constant.HBar())*(density/constant.Rho0());
  return result;
}

double AlvarezRusoCOHPiPDXSec::DeltaSelfEnergyIm(double density)
{
  // Oset, Salcedo, NPA 468(87)631
  // Using eq. (3.5) to relate the energy of the delta with the pion
  // energy used in the parametrization
  
  double E = p_pi.E();
  
  // The parameterization is valid for  85 MeV < tpi < 315
  // above which we take a contant values

  if( E >= (m_pi + (0.315/constant.HBar())) )
  {
    E = m_pi + 0.315/constant.HBar();
  }
  
  double Cq = DeltaSelfEnergyConstant(-5.19, 15.35, 2.06, E) / (1000.0 * constant.HBar());
  double Ca2 = DeltaSelfEnergyConstant(1.06, -6.64, 22.66, E) / (1000.0 * constant.HBar());

  // Ca3 extrapolated to zero at low kin. energies
  double Ca3;
  if( E <= (m_pi + (0.085/constant.HBar())) )
  {
    Ca3 = DeltaSelfEnergyConstant(-13.46, 46.17 , -20.34, (m_pi + 0.085/(1000.0*constant.HBar())));
    Ca3 /= 85.0;
    Ca3 *= (E - m_pi);
  }
  else
  {
    Ca3 = DeltaSelfEnergyConstant(-13.46, 46.17, -20.34, E) / (1000.0 * constant.HBar());
  }
  double alpha = DeltaSelfEnergyConstant(0.382, -1.322, 1.466, E);
  double beta  = DeltaSelfEnergyConstant(-0.038, 0.204, 0.613, E);
  double gamma = 2.0*beta;
  
  double ratio = density / constant.Rho0();
  
  double result = Cq * TMath::Power(ratio, alpha);
  result += Ca2 * TMath::Power(ratio, beta);
  result += Ca3 * TMath::Power(ratio, gamma);
  result *= -1.0;
  
  return result;
}

double AlvarezRusoCOHPiPDXSec::DeltaSelfEnergyConstant(double a, double b, double c, double E)
{
  double x = (E / m_pi) - 1.0;
  return (a*x*x + b*x + c);
}

void AlvarezRusoCOHPiPDXSec::NuclearCurrent(LorentzVector q, LorentzVector pdir, LorentzVector pcrs, LorentzVector ppi, cdouble *jHadCurrent)
//calculates the nuclear current
{
  CVector j1;
  CVector j2;
  CVector j3;
  CVector j4;

  //double ga = constant.GAxial();
  //double fpi = constant.PiDecayConst();
  double mn = constant.NucleonMass();
  double mn2 = mn*mn;
  double mn3 = mn*mn*mn;
  double mdel = constant.DeltaPMass();
  double mdel2 = mdel*mdel;
  double mdel3 = mdel*mdel*mdel;
  double mpi = m_pi;
  double q0 = q.E();
  double q02 = q0*q0;
  double q03 = q0*q0*q0;
  double q04 = q0*q0*q0*q0;
  double q05 = q0*q0*q0*q0*q0;
  double q1 = q.X();
  double q12 = q1*q1;
  double q13 = q1*q1*q1;
  double q14 = q1*q1*q1*q1;
  double q3 = q.Z();
  double q32 = q3*q3;
  double q33 = q3*q3*q3;
  double q34 = q3*q3*q3*q3;
  double pi = constants::kPi;
  double ppim = ppi.P();
  double ppim2 = ppim*ppim;
  double ppi1 = ppi.X();
  double ppi12 = ppi1*ppi1;
  double ppi2 = ppi.Y();
  double ppi22 = ppi2*ppi2;
  double ppi3 = ppi.Z();
  double ppi32 = ppi3*ppi3;
  double ppitr = TMath::Sqrt( ppi12 + ppi22 );
  double ppitr2 = ppitr*ppitr;
  double fs = constant.DeltaNCoupling();
  double rmax = nucleus.RadiusMax();
  
  cdouble I(0,1);
  cdouble twoI(0,2);
  
  double hbarsq = constant.HBar()*constant.HBar();

  double alp;
  if( current == kCC )
  {
    alp = 3.0;
  }
  else
  {
    alp = 1.0;
  }
  
  double mod;
  if(current == kCC)
  {
    mod = 1.0;
  }
  else
  {
    mod = constants::kSqrt2;
  }

  double t = q.mag2() * hbarsq; // in GeV
  double Ma2_Delta = constant.Ma_Delta() * constant.Ma_Delta();
  double Mv2_Delta = constant.Mv_Delta() * constant.Mv_Delta();

  double C3v = 2.05 / ( (1.0 - (t/Mv2_Delta))*(1.0 - (t/Mv2_Delta)) );
  double C4v = (-constant.NucleonMass() / constant.DeltaPMass() ) * C3v;
  double C5v = 0.0;

  if( current == kNC )
  {
    double nc_factor = constant.NCFactor();
    C3v *= nc_factor;
    C4v *= nc_factor;
    C5v *= nc_factor;
  }
  
  double C5a = 1.2 * ( 1.0 - ((constant.CA5_A() * t)/(constant.CA5_B() - t)) ) / 
                                           ( ( 1.0 - (t / Ma2_Delta))*( 1.0 - (t / Ma2_Delta)) );
  double C4a = -C5a / 4.0;
  double C6a = (C5a * mn*mn) / ((mpi*mpi) - (t / hbarsq));
  
  // QE Form Factors
  double F1;
  double F2;
  double FA;
  double FP;
  {
    double mun = -1.913;
    double mup = 2.793;
     
    double MNucleon = constant.NucleonMass() * constant.HBar();
    double MPion = m_pi * constant.HBar();
    double Mv2_Nucleon = constant.Mv_Nucleon()*constant.Mv_Nucleon();
    double Ma2_Nucleon = constant.Ma_Nucleon()*constant.Ma_Nucleon();

    double Q_s = -t;
    double Q_s2 = Q_s* Q_s;
    double Q_s3 = Q_s2*Q_s;
    double Q_s4 = Q_s3*Q_s;
    double Q_s5 = Q_s4*Q_s;
    double Q_s6 = Q_s5*Q_s;
    
    double tau = Q_s / (4.0 * MNucleon*MNucleon);
    
    // parametrization by Budd, Bodek, Arrington (hep-ex/0308005) - BBA2003 formfactors
    // valid up to t = 6 GeV**2
    
    double GEp = 1.0 / (1.0 + 3.253*Q_s + 1.422*Q_s2 + 0.08582*Q_s3 + 0.3318*Q_s4 - 0.09371*Q_s5 + 0.01076*Q_s6);
    double GMp = mup / (1.0 + 3.104*Q_s + 1.428*Q_s2 + 0.1112*Q_s3 - 0.006981*Q_s4 + 0.0003705*Q_s5 - 7.063e-6*Q_s6);
    double GEn = ((-mun * 0.942 * tau) / (1.0 + 4.61*tau)) / ( (1.0 + Q_s/Mv2_Nucleon) * (1.0 + Q_s/Mv2_Nucleon) );
    
    // parametrization of Krutov (hep-ph/0202183)

    double GMn = mun / (1.+3.043*Q_s + 0.8548*Q_s2 + 0.6806*Q_s3 - 0.1287*Q_s4 + 0.008912*Q_s5);
    F1 = ((GEp - GEn) + tau*(GMp - GMn)) / (1.0 + tau);
    F2 = ((GMp - GMn) - (GEp - GEn)) / (1.0 + tau);
  
          FA = 1.0 + ( Q_s / Ma2_Nucleon );
          FA *= FA;
    FA = constant.GAxial() / FA;
    
    FP = (2.0 * MNucleon*MNucleon);
    FP /= ( MPion*MPion + Q_s );
    FP *= FA;
    FP /= constant.NucleonMass();
    
    if( current == kNC )
    {
      double nc_factor = constant.NCFactor();
      F1 *= nc_factor;
      F2 *= nc_factor;
    }
  }
  
  // Get q momentum component perpendicular to pion momentum
  double qper2 = q.P2();
  double tot = ppi.P2();
  double dot = q.Vect().Dot(ppi.Vect());
  if (tot > 0.) qper2 -=dot*dot/tot;
  qper2 = TMath::Max(qper2,0.);
  
  double qper = TMath::Sqrt(qper2);
  double qpar = dot / ppim;
  
  int n = this->sampling * 2;
  
  std::vector<cdouble > empty_row(n, cdouble(0.0,0.0));
  std::vector< std::vector<cdouble > > ordez(4, empty_row);
  std::vector< std::vector<cdouble > > ordez1(4, empty_row);
  std::vector< std::vector<cdouble > > ordez2(4, empty_row);
  std::vector< std::vector<cdouble > > ordez3(4, empty_row);
  std::vector< std::vector<cdouble > > ordez4(4, empty_row);
  
  std::vector< std::vector<cdouble > > ordeb2(4, empty_row);
  // IMPORTANT !!!
  // ORDEB HAS ITS INDICES REVERSED W.R.T THE FORTRAN
  std::vector<cdouble > empty_row_backwards(4, cdouble(0.0,0.0));
  std::vector< std::vector<cdouble > > ordeb(n, empty_row_backwards);
  
  cdouble ppi1d, ppi2d, ppi3d;
  
  std::vector<cdouble > jnuclear(4);
  
  
  for(int i = 0; i != n; ++i)
  {
    double be = nucleus.SamplePoint1(i);
    double bej0 = TMath::BesselJ0( qper*be );
    double bej1 = TMath::BesselJ1( qper*be );
    
    for(int l = 0; l != n; ++l)
    {
      double za = nucleus.SamplePoint2(l);
      double r = TMath::Sqrt(za*za + be*be);
      double r2 = r*r;
      
      double dens = nucleus.Density(i,l);
      
      cdouble exp_i_qpar_za = exp(I*qpar*za);
      
      const cdouble & uwavefunc = (*uwave)[i][l];
      cdouble A = I * ( (*uwaveDr)[i][l] - ((za/r2) * (*uwaveDtheta)[i][l]) ) * (be/r);
      cdouble B = I * ( ((*uwaveDr)[i][l]*za) + ((be*be/r2)*(*uwaveDtheta)[i][l]) ) / r;
      

      // Calculate distorted pion momentum components
      if( (qper == 0.0) && ppitr != 0.0)
      {
        ppi1d = ( q1 * ((ppi12*ppi32)+(ppi22*ppim2)) - q3*ppi1*ppi3*ppitr2 ) /
                  (ppim2*ppitr2)*A*(I*be/2.0) + ppi1/ppim *B*bej0;
        ppi2d = -ppi2*(ppi1*q1+ppi3*q3)/ppim2*A*(I*be/2.)+ppi2/ppim*B*bej0;
        ppi3d = -(q1*ppi1*ppi3-q3*ppitr2)/ppim2*A*(I*be/2.)+ppi3/ppim*B*bej0 ;
      }
      else if( (qper != 0.0) && (ppitr == 0.0) )
      {
        ppi1d = (q1/qper)*A*(I*bej1);
        ppi2d = 0.0;
        ppi3d = B*bej0;
      }
      else if( (qper == 0.0) && (ppitr == 0.0) )
      {
        ppi1d = q1*A*(I*be/2.0);
        ppi2d = 0.0;
        ppi3d = B*bej0;
      }
      else
      {
        ppi1d=(q1*((ppi12*ppi32)+(ppi22*ppim2))-q3*ppi1*ppi3*ppitr2)/(ppim2*ppitr2)/
                                                                                qper*A*(I*bej1)+ppi1/ppim*B*bej0;
        ppi2d = -ppi2 * (ppi1*q1 + ppi3*q3) / ppim2/qper*A*(I*bej1)+ppi2/ppim*B*bej0;
        ppi3d=-(q1*ppi1*ppi3-q3*ppitr2)/ppim2/qper*A*(I*bej1)+ppi3/ppim*B*bej0;
      }
      
      // Calculate the current for four different processes (See Fig 1 of Alvarez-Ruso et al,
      // "Neutral current coherent pion production", arXiv:0707.2172
      // j1 : current for direct delta production
      // j2 : current for Crossed delta production
      // j3 : current for direct nucleon production
      // j4 : current for crossed nucleon production
      j1[0] = -4*(mdel + mn + q0)*((C5a*mdel2*mn2*q0 + C6a*mdel2*q03 + C5a*mn2*(-mn - q0)*q0*(mn + q0) -
         C4a*mdel2*q0*q12 - C4a*mdel2*q0*q32 - C6a*q02*(mn + q0)*(q0*(mn + q0) - q12 - q32))*bej0*uwavefunc +
         ppi1d*(-(C5a*mn2*(-mn - q0)*q1) - C6a*mdel2*q0*q1 + C4a*mdel2*(mn + q0)*q1 +
               C6a*q0*q1*(q0*(mn + q0) - q12 - q32)) + ppi3d*(-(C5a*mn2*(-mn - q0)*q3) - 
           C6a*mdel2*q0*q3 + C4a*mdel2*(mn + q0)*q3 + C6a*q0*q3*(q0*(mn + q0) - q12 - q32)));

      j1[1] = (-4*C6a*mdel3*q02*q1 - 4*C6a*mdel2*mn*q02*q1 + 4*C6a*mdel*mn2*q02*q1 + 4*C6a*mn3*q02*q1 - 
        4*C6a*mdel2*q03*q1 + 8*C6a*mdel*mn*q03*q1 + 12*C6a*mn2*q03*q1 + 4*C6a*mdel*q04*q1 + 
        12*C6a*mn*q04*q1 + 4*C6a*q05*q1 + 4*C4a*mdel2*q02*(mdel + mn + q0)*q1 + 
        4*C5a*mn2*q0*(mn + q0)*(mdel + mn + q0)*q1 - 4*C6a*mdel*mn*q0*q13 - 4*C6a*mn2*q0*q13 - 
        4*C6a*mdel*q02*q13 - 8*C6a*mn*q02*q13 - 4*C6a*q03*q13 - 
        4*C6a*q0*(mn + q0)*(mdel + mn + q0)*q1*q32)*bej0*uwavefunc + 
        ppi1d*(-4*C4a*mdel2*q0*(mn + q0)*(mdel + mn + q0) + 4*C6a*mdel3*q12 + 4*C6a*mdel2*mn*q12 + 
          4*C6a*mdel2*q0*q12 - 4*C6a*mdel*mn*q0*q12 - 4*C6a*mn2*q0*q12 - 4*C6a*mdel*q02*q12 - 
        8*C6a*mn*q02*q12 - 4*C6a*q03*q12 + 4*C6a*mdel*q14 + 4*C6a*mn*q14 + 4*C6a*q0*q14 - 
        4*C5a*mn2*(mdel + mn + q0)*(mdel2 + q12) + 4*C4a*mdel2*(mdel + mn + q0)*q32 + 
        4*C6a*(mdel + mn + q0)*q12*q32) + ppi2d*(twoI*C4v*mdel2*(q0*(mn + q0) - q12)*q3 + 
        twoI*C3v*mdel*mn*(2*mdel2 + 2*mdel*mn - q0*(mn + q0) + q12)*q3 - 
           twoI*mdel*(C4v*mdel - C3v*mn)*q33) + ppi3d*(-4*C4a*mdel2*(mdel + mn + q0)*q1*q3 - 
        4*C5a*mn2*(mdel + mn + q0)*q1*q3 + 4*C6a*(mdel + mn + q0)*q1*(mdel2 - q0*(mn + q0) + q12)*q3 + 
        4*C6a*(mdel + mn + q0)*q1*q33) + ppi2d*twoI*C5v*mdel2*mn*q0*q3;
    
      j1[2] = -2.*I*(-2.*I*C5a*mdel*mn2*ppi2d*(mdel + mn + q0) - 
        twoI*C4a*mdel*ppi2d*(mdel + mn + q0)*(q0*(mn + q0) - q12 - q32) - 
        (ppi3d*q1 - ppi1d*q3)*(C4v*mdel*(q0*(mn + q0) - q12 - q32) + 
        C3v*mn*(2*mdel2 + 2*mdel*mn - q0*(mn + q0) + q12 + q32)))*mdel + 
        twoI*C5v*mdel*mn*q0*(ppi3d*q1 - ppi1d*q3)*mdel;

      j1[3] = (4*C4a*mdel2*q02*(mdel + mn + q0)*q3 - 4*C6a*mdel2*q02*(mdel + mn + q0)*q3 + 
        4*C5a*mn2*q0*(mn + q0)*(mdel + mn + q0)*q3 + 4*C6a*q0*(mn + q0)*(mdel + mn + q0)*
        (q0*(mn + q0) - q12)*q3 - 4*C6a*q0*(mn + q0)*(mdel + mn + q0)*q33)*bej0*uwavefunc + 
        ppi2d*(-twoI*mdel*q1*(C4v*mdel*(q0*(mn + q0) - q12) + 
        C3v*mn*(2*mdel2 + 2*mdel*mn - q0*(mn + q0) + q12)) + twoI*mdel*
        (C4v*mdel - C3v*mn)*q1*q32) + ppi1d*(-4*C4a*mdel2*(mdel + mn + q0)*q1*q3 + 
        4*C6a*mdel2*(mdel + mn + q0)*q1*q3 - 4*C5a*mn2*(mdel + mn + q0)*q1*q3 - 
        4*C6a*(mdel + mn + q0)*q1*(q0*(mn + q0) - q12)*q3 + 4*C6a*(mdel + mn + q0)*q1*q33) + 
        ppi3d*(-4*C4a*mdel2*mn*q0*(mdel + mn + q0) - 4*C4a*mdel2*(mdel + mn + q0)*(q0 - q1)*(q0 + q1) + 
        4*C6a*(mdel + mn + q0)*(mdel2 - q0*(mn + q0) + q12)*q32 + 4*C6a*(mdel + mn + q0)*q34 - 
        4*C5a*mn2*(mdel + mn + q0)*(mdel2 + q32)) + -twoI*C5v*mdel2*mn*ppi2d*q0*q1;
      
      // Crossed Delta
      
      j2[0]=-4*(mdel + mn - q0)*((C5a*mdel2*mn2*q0 - C5a*mn2*(mn - q0)*(mn - q0)*q0 + C6a*mdel2*q03 - 
        C4a*mdel2*q0*q12 - C4a*mdel2*q0*q32 - C6a*(mn - q0)*q02*((mn - q0)*q0 + q12 + q32))*bej0*uwavefunc + 
        ppi1d*(-(C4a*mdel2*mn*q1) - C5a*mn2*(mn - q0)*q1 + C4a*mdel2*q0*q1 - 
        C6a*mdel2*q0*q1 - C6a*q0*q1*((mn - q0)*q0 + q12 + q32)) + 
        ppi3d*(-(C4a*mdel2*mn*q3) - C5a*mn2*(mn - q0)*q3 + C4a*mdel2*q0*q3 - 
        C6a*mdel2*q0*q3 - C6a*q0*q3*((mn - q0)*q0 + q12 + q32)));

      j2[1]=(-4*C5a*mn2*(mn - q0)*(mdel + mn - q0)*q0*q1 - 4*C6a*mdel3*q02*q1 - 
        4*C6a*mdel2*mn*q02*q1 + 4*C6a*mdel*mn2*q02*q1 + 4*C6a*mn3*q02*q1 + 
        4*C4a*mdel2*(mdel + mn - q0)*q02*q1 + 4*C6a*mdel2*q03*q1 - 
        8*C6a*mdel*mn*q03*q1 - 12*C6a*mn2*q03*q1 + 4*C6a*mdel*q04*q1 + 
        12*C6a*mn*q04*q1 - 4*C6a*q05*q1 + 4*C6a*mdel*mn*q0*q13 + 
        4*C6a*mn2*q0*q13 - 4*C6a*mdel*q02*q13 - 8*C6a*mn*q02*q13 + 
        4*C6a*q03*q13 + 4*C6a*(mn - q0)*(mdel + mn - q0)*q0*q1*q32)*bej0*uwavefunc + 
        ppi1d*(-4*C5a*mdel2*mn2*(mdel + mn - q0) + 4*C4a*mdel2*mn*(mdel + mn - q0)*q0 - 
        4*C4a*mdel2*(mdel + mn - q0)*q02 + 4*C6a*mdel3*q12 + 
        4*C6a*mdel2*mn*q12 - 4*C5a*mn2*(mdel + mn - q0)*q12 - 
        4*C6a*mdel2*q0*q12 + 4*C6a*mdel*mn*q0*q12 + 4*C6a*mn2*q0*q12 - 
        4*C6a*mdel*q02*q12 - 8*C6a*mn*q02*q12 + 4*C6a*q03*q12 + 
        4*C6a*mdel*q14 + 4*C6a*mn*q14 - 4*C6a*q0*q14 + 
        4*C4a*mdel2*(mdel + mn - q0)*q32 + 4*C6a*(mdel + mn - q0)*q12*q32) + 
        ppi2d*(twoI*mdel*(mdel*q0*(-((C4v + C5v)*mn) + C4v*q0) + 
        C3v*mn*(2*mdel*(mdel + mn) + mn*q0 - q02))*q3 + 
        twoI*mdel*(-(C4v*mdel) + C3v*mn)*q12*q3 - twoI*mdel*(C4v*mdel - C3v*mn)*q33) 
        + ppi3d*(-4*C4a*mdel2*(mdel + mn - q0)*q1*q3 - 
        4*C5a*mn2*(mdel + mn - q0)*q1*q3 + 4*C6a*(mdel + mn - q0)*(mdel2 + (mn - q0)*q0)*q1*q3 + 
        4*C6a*(mdel + mn - q0)*q13*q3 + 4*C6a*(mdel + mn - q0)*q1*q33);

      j2[2]=-(twoI*ppi2d*(-twoI*C5a*mdel*mn2*(mdel + mn - q0) +
        twoI*C4a*mdel*(mdel + mn - q0)*((mn - q0)*q0 + q12 + q32)) - 
        twoI*ppi3d*q1*(C3v*mn*(2*mdel*(mdel + mn) + mn*q0 - q02 + q12 + q32) - 
        mdel*(C5v*mn*q0 + C4v*((mn - q0)*q0 + q12 + q32))) + 
        twoI*ppi1d*q3*(C3v*mn*(2*mdel*(mdel + mn) + mn*q0 - q02 + q12 + q32) - 
        mdel*(C5v*mn*q0 + C4v*((mn - q0)*q0 + q12 + q32))))*mdel;

      j2[3] = (-4*C5a*mn2*(mn - q0)*(mdel + mn - q0)*q0*q3 +
        4*C4a*mdel2*(mdel + mn - q0)*q02*q3 - 4*C6a*mdel2*(mdel + mn - q0)*q02*q3 + 
        4*C6a*(mn - q0)*(mdel + mn - q0)*q0*((mn - q0)*q0 + q12)*q3 + 
        4*C6a*(mn - q0)*(mdel + mn - q0)*q0*q33)*bej0*uwavefunc + 
        ppi2d*(-twoI*mdel*q1*(C3v*mn*(2*mdel*(mdel + mn) + mn*q0 - q02 + q12) - 
        mdel*(q0*((C4v + C5v)*mn - C4v*q0) + C4v*q12)) + 
        twoI*mdel*(C4v*mdel - C3v*mn)*q1*q32) + 
        ppi1d*(-4*C4a*mdel2*(mdel + mn - q0)*q1*q3 + 4*C6a*mdel2*(mdel + mn - q0)*q1*q3 - 
        4*C5a*mn2*(mdel + mn - q0)*q1*q3 + 
        4*C6a*(mdel + mn - q0)*q1*((mn - q0)*q0 + q12)*q3 + 
        4*C6a*(mdel + mn - q0)*q1*q33) + 
        ppi3d*(-4*C5a*mdel2*mn2*(mdel + mn - q0) + 
        4*C4a*mdel2*(mdel + mn - q0)*((mn - q0)*q0 + q12) - 4*C5a*mn2*(mdel + mn - q0)*q32 + 
        4*C6a*(mdel + mn - q0)*(mdel2 + (mn - q0)*q0 + q12)*q32 + 
        4*C6a*(mdel + mn - q0)*q34);
      
      //
      // Direct Nucleon
      
      j3[0]=-2.0*(FA - FP*q0)*(q02*bej0*uwavefunc - ppi1d*q1 - ppi3d*q3);
      
      j3[1] = 2.0*(-(FA*q0*q1) + FP*q02*q1) * bej0 * 
        uwavefunc + 2.0*ppi1d*(2.0*FA*mn + FA*q0 - FP*q12) -
        twoI*(F1 + F2)*ppi2d*q3 - 2*FP*ppi3d*q1*q3;

      j3[2]=twoI*(-I*FA*ppi2d*(2.0*mn + q0) - (F1 + F2)*(ppi3d*q1 - ppi1d*q3));

      j3[3]=twoI*(F1 + F2)*ppi2d*q1 - 2.0*FP*ppi1d*q1*q3 + 
        2.0*(-(FA*q0*q3) + FP*q02*q3)*bej0*uwavefunc + 2.0*ppi3d*(2.0*FA*mn + FA*q0 - FP*q32);
      
      //
      // Crossed Nucleon
      
      j4[0] = 2.0*(FA + FP*q0)*(q02  *bej0*uwavefunc - ppi1d*q1 - ppi3d*q3);

      j4[1] = -2.0*(-(FA*q0*q1) - FP*q02*q1)*bej0*uwavefunc - 2.0*ppi1d*(-2.0*FA*mn + FA*q0 + FP*q12) - 
        twoI*(F1 + F2)*ppi2d*q3 - 2.0*FP*ppi3d*q1*q3;

      j4[2] = twoI*(-I*FA*ppi2d*(2.0*mn - q0) - (F1 + F2)*(ppi3d*q1 - ppi1d*q3));

      j4[3] = twoI*(F1 + F2)*ppi2d*q1 - 2.0*FP*ppi1d*q1*q3 + 2.0*(FA*q0*q3 + FP*q02*q3)*bej0*uwavefunc + 
        2.0*ppi3d*(2.0*FA*mn - FA*q0 - FP*q32);

      cdouble pre_factor_1 = mod * I *
        (fs/mpi) / constants::kSqrt3 *
        exp_i_qpar_za *
        (alp*dens+dens) *
        DeltaCouplingInMed(pdir,ppi,dens,dens) *
        pi/(3.0*mn*mn*mdel*mdel)*F_direct_delta;

      cdouble pre_factor_2 = mod * I * (fs/mpi) / constants::kSqrt3 *exp_i_qpar_za *
        (dens + alp*dens)*pi/(3.*mn2*mdel2)*F_cross_delta *
        DeltaCouplingInMed(pcrs,ppi,dens,dens);

      double PreFacMult = (constant.GAxial()/constants::kSqrt2/constant.PiDecayConst());
      cdouble pre_factor_3 = 1./mod*(-I)*PreFacMult*exp_i_qpar_za*
                             dens*NucleonPropagator(pdir)*pi*F_direct_nucleon;
      cdouble pre_factor_4 =  1./mod*(-I)*PreFacMult*exp_i_qpar_za*
                             dens*NucleonPropagator(pcrs)*pi*F_cross_nucleon;
                      
      for(int m = 0; m != 4; ++m)
      {
        ordez1[m][l] = pre_factor_1 * j1[m];
        ordez2[m][l] = pre_factor_2 * j2[m];
        ordez3[m][l] = pre_factor_3 * j3[m];
        ordez4[m][l] = pre_factor_4 * j4[m];
        
        ordez[m][l] = ordez1[m][l] + ordez2[m][l] + ordez3[m][l] + ordez4[m][l];

      }
    }
    
    // IMPORTANT !!!
    // ORDEB HAS ITS INDICES REVERSED W.R.T THE FORTRAN
    
    std::vector<cdouble >result(4);
    integrationtools::RG202D(-rmax, rmax, 2, 0, 3, ordez, sampling, ordeb[i]);


    
    for(unsigned int z = 0; z != 4; ++z)
    {
      ordeb[i][z] *= be;
    }
  }
  
  for(int z = 0; z != n; ++z)
  {
    for(unsigned int y = 0; y != 4; ++y)
    {
      ordeb2[y][z] = ordeb[z][y];
    }
  }
  
  integrationtools::RG202D(0.0, rmax, 2, 0, 3, ordeb2, sampling, jnuclear);
  for (int i=0; i<4; i++) {
    cdouble & jn = jnuclear[i];
    *(jHadCurrent+i) = cdouble(jn);
  }
}


cdouble AlvarezRusoCOHPiPDXSec::DeltaCouplingInMed(LorentzVector delta_momentum, LorentzVector pion_momentum, double density_p, double density_n) 
{
  cdouble gdmed;
  cdouble s_delta (delta_momentum.mag2(),0);
  cdouble I(0,1);
  if( real(s_delta) < ((constant.NucleonMass() + m_pi)*(constant.NucleonMass() + m_pi)) )
  {
    gdmed = 1.0 / ( s_delta - constant.DeltaPMass()*constant.DeltaPMass() ) ;
  }
  else if( delta_momentum.E() < 0.0 )
  {
    gdmed = 0.0;
  }
  else
  {
    cdouble sqrt_delta = sqrt(s_delta);
    double density = density_p + density_n;
    double gamdpb = DeltaWidthPauliBlocked(delta_momentum, density);
    double real = DeltaSelfEnergyRe(density);
    double imaginary = DeltaSelfEnergyIm(density);
    double ofshel = PiDecayVertex(pion_momentum, constant.DeltaPMass());

    cdouble part_1 = sqrt_delta - constant.DeltaPMass() + (ofshel*ofshel*(I*gamdpb)/2.0) - real - I*imaginary;
    cdouble part_2 = sqrt_delta + constant.DeltaPMass();
    
    gdmed = 1.0 / (part_1 * part_2);
  }
  
  return gdmed;
}

cdouble AlvarezRusoCOHPiPDXSec::NucleonPropagator(LorentzVector nucleon_momentum)
{
  // relativistic nucleon propagator (its denominator)
  
  cdouble gn( nucleon_momentum.mag2() - constant.NucleonMassSq() ,
           ( constant.NucleonMass()*10.0 / (1000.0*constant.HBar()) ) );
  gn = 1.0 / gn;
  
  return gn;
}

} // namespace alvarezruso
} // namespace genie
