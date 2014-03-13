//____________________________________________________________________________
/*!

\class    genie::AtharSingleKaonPXSec

\brief    Differential cross section for single kaon production.

\ref      Physical Review D82 (2010) 033001 (arXiv:1004.5484 [hep-ph])

\author   Chris Marshall and Martti Nirkko

\created  2014-02-14

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _ATHAR_SINGLE_KAON_PXSEC_H_
#define _ATHAR_SINGLE_KAON_PXSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class XSecIntegratorI;

class AtharSingleKaonPXSec : public XSecAlgorithmI {

public:
  AtharSingleKaonPXSec();
  AtharSingleKaonPXSec(string config);
  virtual ~AtharSingleKaonPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;
  bool   ValidKinematics (const Interaction * i) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  void LoadConfig(void);

  const XSecIntegratorI * fXSecIntegrator;  ///< cross section integrator

};

}       // genie namespace
#endif  // _ATHAR_SINGLE_KAON_PXSEC_H_
