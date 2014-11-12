//____________________________________________________________________________
/*!

\class    genie::DISInteractionListGenerator

\brief    Concrete implementations of the InteractionListGeneratorI interface.
          Generate a list of all the Interaction (= event summary) objects that
          can be generated by the DIS EventGenerator.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 13, 2005

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DIS_INTERACTION_LIST_GENERATOR_H_
#define _DIS_INTERACTION_LIST_GENERATOR_H_

#include <map>

#include "EVGCore/InteractionListGeneratorI.h"

using std::multimap;

namespace genie {

class Interaction;

class DISInteractionListGenerator : public InteractionListGeneratorI {

public :
  DISInteractionListGenerator();
  DISInteractionListGenerator(string config);
 ~DISInteractionListGenerator();

  // implement the InteractionListGeneratorI interface
  InteractionList * CreateInteractionList(const InitialState & init) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfigData(void);

  multimap<int,bool> GetHitQuarks(const Interaction * interaction) const;

  bool fIsCC;
  bool fIsNC;
  bool fIsEM;
  bool fSetHitQuark;
  bool fIsCharm;
};

}      // genie namespace

#endif // _DIS_INTERACTION_LIST_GENERATOR_H_
