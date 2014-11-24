//____________________________________________________________________________
/*!

\program gCCohARTester

\brief   Sets up a Alvarez-Ruso object, for testing

         Syntax :
           gCCohARTester -t target_list

         Options :
           -t A list of comma-delimited targets to calculate the wavefunction for.
          Valid targets are carbon, oxygen, aluminium, calcium and iron.

     Notes :
       If no targets are specified, lookup tables will be made for all valid
       target lists.

         Examples :

           1) shell% gcohartest 

          will generate lookup tables for all nuclear target types.

           2) shell% gcohartest -t C,Fe

          will generated a lookup table for C and Fe only

\author  Steve Dennis <stephen.dennis \at warwick.ac.uk>
         University of Warwick, RAL

\created August 06, 2013

\cpright Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include <TSystem.h>
#include <TMath.h>
#include <TRandom3.h>

#include "Base/XSecAlgorithmI.h"
#include "Conventions/XmlParserStatus.h"
#include "Messenger/Messenger.h"
#include "Utils/XSecSplineList.h"
#include "Utils/StringUtils.h"
#include "Utils/SystemUtils.h"
#include "Utils/CmdLnArgParser.h"
#include "AlvarezRuso/AlvarezRusoCOHPiPDXSec.h"

using std::string;
using std::vector;
using std::ostringstream;

using namespace genie;
using namespace genie::alvarezruso;

void           GetCommandLineArgs (int argc, char ** argv);
void           PrintSyntax        (void);

//User-specified options:
string         gOutFile;   ///< output XML file
vector<string> gTargets;  ///< list of all input files

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs(argc,argv);

  int nTargets = gTargets.size();
  std::vector<std::string> NucleusVector;
  if (gTargets.size() < 1) 
  { 
    LOG("gcohartest",pINFO) << "No target defined. We will calculate them all"; 
    nTargets = 5;
    NucleusVector.push_back("C12");
    NucleusVector.push_back("O16");
    NucleusVector.push_back("Al27");
    NucleusVector.push_back("Ca40");
    NucleusVector.push_back("Fe56");
  }
  else {
    std::vector<string>::iterator it;
    for (it = gTargets.begin(); it != gTargets.end(); ++it) {
        std::cout << (*it) << std::endl;
        std::string target = genie::utils::str::ToUpper((*it));
        if ( target == "C" || target == "C12" ) { 
            NucleusVector.push_back("C12");
        } else if ( target == "O" || target == "O16" ) {
                    NucleusVector.push_back("O16");
        } else if ( target == "AL" || target == "AL27" ) {
                    NucleusVector.push_back("Al27");
        } else if ( target == "CA" || target == "CA40" ) {
                    NucleusVector.push_back("Ca40");
        } else if ( target == "FE" || target == "FE56" ) {
                    NucleusVector.push_back("Fe56");
        } else {
           LOG("gcohartest",pINFO) << "Target " << *it << " is not currently defined in the Alvarez-Ruso model";
        }
    }
  }
  if (NucleusVector.size() < 1) {
    LOG("gcohartest",pINFO) << " No valid target type requested. Exiting ";
    return 0;
  }
  
  flavour_t flavours[] = {kMu,kE};
  nutype_t types[]     = {kNu,kAntiNu};
  current_t currents[] = {kCC,kNC};
  
  TRandom3 tr;// Fixed seed when initiated without constant, always the same
  
  int Z,A;
  double Enu = 0.6;
  unsigned int n_iter = 50;
  
  for (unsigned int iflav = 0 ; iflav < 2 ; ++iflav) {
    for (unsigned int itype = 0 ; itype < 2 ; ++itype) {
      for (unsigned int icurr = 0 ; icurr < 2 ; ++icurr) {
        for (unsigned int i=0; i<NucleusVector.size(); ++i) {
          if      (NucleusVector[i].compare("C12" )==0) {Z= 6; A=12;}
          else if (NucleusVector[i].compare("O16" )==0) {Z= 8; A=16;}
          else if (NucleusVector[i].compare("Al27")==0) {Z=13; A=27;}
          else if (NucleusVector[i].compare("Ca40")==0) {Z=20; A=40;}
          else if (NucleusVector[i].compare("Fe56")==0) {Z=26; A=56;}
          else { std::cout << "Unknown Nucleus. Exiting" << std::endl; exit(2); }
          
          std::cout<<NucleusVector[i]<<std::endl;
          AlvarezRusoCOHPiPDXSec coh(Z, A,currents[icurr], flavours[iflav],types[itype]);
          for (unsigned int j=0 ; j < n_iter ; j++) {
            double Elep = tr.Rndm() * Enu;
            double theta_lep = tr.Rndm() * TMath::Pi();
            double phi_lep   = tr.Rndm() * TMath::Pi() * 2;
            double theta_pi  = tr.Rndm() * TMath::Pi();
            double phi_pi    = tr.Rndm() * TMath::Pi() * 2;
            double dxsec = coh.DXSec(Enu,Elep,theta_lep, phi_lep, theta_pi, phi_pi);
            std::cout<<dxsec<<"\n";
          }
        }// loop over nuclear target list
      }
    }
  }
 
  return 0;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gcohartest", pNOTICE) << "Parsing command line arguments";
  
  CmdLnArgParser parser(argc,argv);
  std::cout << argv << std::endl;
  if( parser.OptionExists('t') ) {
    LOG("gcohartest", pINFO) << "Reading input files";
    string targets = parser.ArgAsString('t');
    if(targets.find(",") != string::npos) {
      // split the comma separated list
      gTargets = utils::str::Split(targets, ",");
    } else {
      // there is just one file
      gTargets.push_back(targets);
    }
  }
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gcohartest", pNOTICE)
   << "\n\n" << "Syntax:" << "\n"
   << "   gcohartest  [-t target_list] [--table] \n";
}
//____________________________________________________________________________
