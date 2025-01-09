//-------------------------------------------------------------------------//
//                                                                         //
//                                CAnalyse.h                               //
//                               Version 1.0                               //
//                        Matthieu Lebois November 2021                    //
//                                                                         //
//  This file contain the definition of the function for a root file       //
//  and the declaration to extract some physical results from the data     //
//  contained in the root file                                             //
//                                                                         //
//-------------------------------------------------------------------------//

#ifndef CEVENTBUILDER_H
#define CEVENTBUILDER_H 1

// Call for the principal parts of the functions
#include "Analyse_Main.h"
#include "CExperiment.h"
#include "CDetectors.h"
#include "CDetecteurType.h"
#include "colormod.h"
#include "CAnalyse.h"
#include "CHit.h"
#include "CHitCollection.h"

// Call for standard C/C++ Libraries
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <map>
#include <numeric>
#include <thread>
#include <algorithm>

// Call for ROOT Libraries
#include "TSpectrum.h"
#include "TApplication.h"
#include "TFile.h"
#include "TCutG.h"
#include "TMath.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH3.h"
#include "TH3F.h"
#include "THStack.h"
#include "TObject.h"
#include "TTree.h"
#include "TBasket.h"
#include "TString.h"
#include "TBranch.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TMath.h"

// Definition of the standard types
using label_type = UShort_t;
using tm_Caltype = Double_t;
using tm_Uncaltype = ULong64_t;
using nrj_Ucaltype = Int_t;
using nrj_Caltype = Double_t;
using pu_type = Bool_t;

// I call the other classes from the analyse program
class CExperiment;
class CDetectors;
//template class CHit<tm_Uncaltype,nrj_Caltype>;

// Definition of the I/O function
using namespace std;

// Function for event building

// Test to reconstruct Ionization Chamber events in the context of JRC Tests
int ICBuilder(CExperiment *experiment, Double_t deltaT);

// Event builder for fold >=2 gamma coincidences
int GammaCoinc(CExperiment *experiment, Double_t deltaT, Char_t gammafoldmin);


#endif //TANALYSEF_H
