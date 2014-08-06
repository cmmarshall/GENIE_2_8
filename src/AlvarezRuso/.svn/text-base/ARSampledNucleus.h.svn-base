//____________________________________________________________________________
/*!

\class    genie::alvarezruso::ARSampledNucleus

\brief    Nucleus class for Alvarez-Ruso Coherent Pion Production xsec

\ref      

\author   Steve Dennis
          University of Warwick, Rutherford Appleton Laboratory

\created  05/12/2013

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________
#ifndef _AR_NUCLEUS_H_
#define _AR_NUCLEUS_H_

namespace genie
{
namespace alvarezruso
{

class ARSampledNucleus
{
  public:
    
    ARSampledNucleus(const unsigned int& ZNumber, const unsigned int& ANumber, const bool debug = false);
    ARSampledNucleus();
    
    ~ARSampledNucleus();
    
    unsigned int A() const  {  return a;  }
    
    unsigned int Z() const  {  return z;  }
    
    unsigned int N() const  {  return (a-z);  }
    
    double Density(const int i, const int j) const;
    double Radius (const int i, const int j) const;
    
    double RadiusMax() const
    {
      return r_max;
    }
    double SamplePoint1(const unsigned int i) const // absib in original fortran
    {
      return sample_points_1[i];
    }
    double SamplePoint2(const unsigned int i) const // absiz in original fortran
    {
      return sample_points_2[i];
    }
    double SampleWeight1(const unsigned int i) const
    {
      return sample_weights_1[i];
    }
    double SampleWeight2(const unsigned int i) const
    {
      return sample_weights_2[i];
    }
    
    void SetSampling(const unsigned int & nsamp);
    
  private:
  
    void Fill();
    void FillSamplePoints();
    void FillDensities();
    // Members
    
    unsigned int sampling;
    bool debug_;
    unsigned int n_densities;
    
    unsigned int a;
    unsigned int z;
    
    double r_max;
    double** radii;
    double** densities;
    double* sample_points_1;  // absib: 0 < r < r_max
    double* sample_points_2;  // absiz: -r_max < r < r_max
    double* sample_weights_1;
    double* sample_weights_2;

};

} //namespace alvarezruso
} //namespace genie

#endif
