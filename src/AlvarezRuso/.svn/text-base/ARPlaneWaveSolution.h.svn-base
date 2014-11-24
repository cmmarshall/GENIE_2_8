//____________________________________________________________________________
/*!

\class    genie::alvarezruso::ARPlaneWaveSolution

\brief    Plane wave wavefunction solution for Alvarez-Ruso Coherent Pion Production xsec

\ref      

\author   Steve Dennis
          University of Warwick, Rutherford Appleton Laboratory

\created  05/12/2013

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _AR_PLANEWAVE_SOLUTION_H_
#define _AR_PLANEWAVE_SOLUTION_H_

#include <complex>

#include "AlvarezRuso/ARWFSolution.h"

namespace genie
{
namespace alvarezruso
{

class ARPlanewaveSolution: public ARWFSolution
{
  public:
    
    ARPlanewaveSolution(bool debug = false);
    virtual ~ARPlanewaveSolution();
    
    virtual std::complex<double>  Element(const double radius, const double cosine_rz, const double p_pion,const pUnits_t piUnit);
    
    virtual void Solve();
};

} //namespace alvarezruso
} //namespace genie

#endif
