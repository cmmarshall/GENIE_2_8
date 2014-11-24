//____________________________________________________________________________
/*!

\class    genie::ASKKinematicsGenerator

\brief    Generates values for the kinematic variables describing coherent 
          neutrino-nucleus pion production events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 03, 2004

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _ASK_KINEMATICS_GENERATOR_H_
#define _ASK_KINEMATICS_GENERATOR_H_

#include "EVGModules/KineGeneratorWithCache.h"
#include "Utils/Range1.h"

class TF3;

namespace genie {

class ASKKinematicsGenerator : public KineGeneratorWithCache {

public :
  ASKKinematicsGenerator();
  ASKKinematicsGenerator(string config);
 ~ASKKinematicsGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

public:
  // methods to load sub-algorithms and config data from the Registry
  void LoadConfig (void);

  // different kinematics calculators for different models
  void   CalculateKin_AtharSingleKaon(GHepRecord * event_rec) const;

  //double MaxXSec (const GHepRecord * ev) const;
  //double MaxXSec (const Interaction * in) const;
  double ComputeMaxXSec (const Interaction * in) const;

  // overload KineGeneratorWithCache method to get energy
  double Energy         (const Interaction * in) const;

//  mutable TF3 * fEnvelope; ///< 3-D envelope used for importance sampling
};

}      // genie namespace
#endif // _ASK_KINEMATICS_GENERATOR_H_
