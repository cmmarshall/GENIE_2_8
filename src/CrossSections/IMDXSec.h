//____________________________________________________________________________
/*!

\class    genie::IMDXSec

\brief    Computes the Inverse Muon Decay cross section.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Fabruary 14, 2005

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _IMD_XSEC_H_
#define _IMD_XSEC_H_

#include "Base/XSecIntegratorI.h"

namespace genie {

class IMDXSec : public XSecIntegratorI {

public:
  IMDXSec();
  IMDXSec(string config);
  virtual ~IMDXSec();

  //! XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

  //! overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig (void);
};

}       // genie namespace
#endif  // _IMD_XSEC_H_
