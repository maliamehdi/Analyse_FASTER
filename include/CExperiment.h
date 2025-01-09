//-------------------------------------------------------------------------//
//                                                                         //
//                               CExperiment.h                             //
//                               Version 1.0                               //
//                        Matthieu Lebois November 2021                    //
//                                                                         //
//  This file contain the definition of the function to define meta        //
//  Parameter to explicitely create everything needed to analyse an        //
//  experiment                                                             //
//                                                                         //
//-------------------------------------------------------------------------//

#ifndef CEXPERIMENT_H
#define CEXPERIMENT_H 1

#include "TString.h"
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <string>
#include "CDetectors.h"
#include "CDetecteurType.h"
#include "colormod.h"



// Definition of the I/O function
using namespace std;

class CExperiment{
public:
  CExperiment( TString name );
  virtual ~CExperiment();

public:
  int nthreads = 1;
  TString Experiencename;
  // Definition of the file name where to find detectors list
  TString  PIDfilename;
  TString  PIDlistname;
  TString  Typelistname;
  TString  Runlist_filename;
  TString  Parameter_folder;
  TString  NRJCalibration_filename;
  TString  TimeCalibration_filename;
  TString  Timewindows_filename;
  TString  PARISangles_filename;
  TString  FileDirectory_IN;
  TString  FileDirectory_OUT;
  std::string  DataTreeName;
  std::vector<std::string>::iterator it_Outputfilenames;
  std::vector<std::string> Outputfilenames;
  std::vector<std::string>::iterator it_Datafilenames;
  std::vector<std::string> Datafilenames;
  Int_t    Nbrofdetectors = 0;
  Int_t    NbrofADC = 0;
  Int_t    NbrofQDC = 0;
  Int_t    NbrofRF = 0;
  Int_t    NbrofSAMPLER = 0;
  Int_t    NbrofModules = 0;
  Int_t    NbrofTChains = 0;
  Bool_t   isMultiChain = kFALSE;
  Int_t    maxLabel = 0;
  Bool_t   isQDC2 = kFALSE;
  Bool_t   isQDC3 = kFALSE;
  Bool_t   isQDC4 = kFALSE;
  Bool_t   isPARISin = kTRUE;

  // Mapping between FASTER label and spectrum number
  std::map<TString, int> DetType2nbr;std::map<TString, int> DetName2nbr;
  std::vector<int> it_Label2Detnbr;
  std::vector<int> Label2Detnbr;
  std::vector<CDetectors*>::iterator it_Detectors;
  std::vector<CDetectors*> tab_Detectors;
  CDetectors *Ref_detector = nullptr;

  // Data I want to manipulate
  std::vector<TChain *>::iterator it_chained_oak;
  std::vector<TChain *> chained_oak;

public:
  void SetThreadsNbr(const int &n){nthreads = n;};
  void SetExperienceName(const TString &name ){Experiencename = name;};
  void SetPIDfilename(const TString &name ){PIDfilename = name;};
  void SetPIDlistname(const TString &name ){PIDlistname = name;};
  void SetTypelistname(const TString &name ){Typelistname = name;};
  void SetRUNlistname(const TString &name );
  void SetParameterFolder(const TString &name ){Parameter_folder=name;};
  void SetNRJCalibration_filename(const TString &name ){NRJCalibration_filename = name;};
  void SetTimeCalibration_filename(const TString &name ){TimeCalibration_filename = name;};
  void SetTimewindows_filename(const TString &name ){Timewindows_filename = name;};
  void SetPARISAngle_filename(const TString &name ){PARISangles_filename = name;};
  void SetFileDirectory_IN(const TString &name ){FileDirectory_IN = name;};
  void SetFileDirectory_OUT(const TString &name ){FileDirectory_OUT = name;};
  void SetDataTreeName(const std::string &name ){DataTreeName = name;};
  Int_t SetNbrofdetectors(const TString name );
  void SetNbrofdetectors(const Int_t &nbr){Nbrofdetectors=nbr;};
  void SetNbrofADC(const Int_t &nbr ){NbrofADC=nbr;};
  void SetNbrofQDC(const Int_t &nbr ){NbrofQDC=nbr;};
  void SetNbrofRF(const Int_t &nbr ){NbrofRF=nbr;};
  void SetNbrofModules(const Int_t &nbr ){NbrofModules=nbr;};
  void SetNbrofTChains(const Int_t &nbr ){NbrofTChains=nbr;};
  void SetQDCNbr(const int &nbr);
  void SetThereIsPARIS(const Bool_t &isPARIS) {isPARISin = isPARIS;};

// Experiment variables Get functions
public:
  int      GetThreadsNbr() const {return nthreads;};
  TString  GetExperienceName() const {return Experiencename;};
  TString  GetPIDfilename() const {return PIDfilename;};
  TString  GetPIDlistname() const {return PIDlistname;};
  TString  GetTypelistname() const {return Typelistname;};
  TString  GetRUNlistname() const {return Runlist_filename;};
  TString  GetParameterFolder() const {return Parameter_folder;};
  TString  GetNRJCalibration_filename() const {return NRJCalibration_filename;};
  TString  GetTimeCalibration_filename() const {return TimeCalibration_filename;};
  TString  GetPARISAngle_filename() const {return PARISangles_filename;};
  TString  GetTimewindows_filename() const {return Timewindows_filename;};
  TString  GetFileDirectory_IN() const {return FileDirectory_IN;};
  TString  GetFileDirectory_OUT() const {return FileDirectory_OUT;};
  std::string  GetDataTreeName() const {return DataTreeName;};
  std::vector<std::string> GetOutputFileNames() const {return Outputfilenames;};
  std::vector<std::string> GetDataFileNames() const {return Datafilenames;};
  Int_t    GetNbrofdetectors() const {return Nbrofdetectors;};
  Int_t    GetNbrofADC() const {return NbrofADC;};
  Int_t    GetNbrofQDC() const {return NbrofQDC;};
  Int_t    GetNbrofRF() const {return NbrofRF;};
  Int_t    GetNbrofModules() const {return NbrofModules;};
  Int_t    GetNbrofTChains() const {return NbrofTChains;};
  Bool_t   GetisQDC2() const {return isQDC2;}
  Bool_t   GetisQDC3() const {return isQDC3;}
  Bool_t   GetisQDC4() const {return isQDC4;}
  Bool_t   IsTherePARIS() const {return isPARISin;};

// Detectors variables
public:
  int SetDetectors( const TString inputfilename, const TString outputfilename, const TString Typefilename );
  std::vector<CDetectors*> GetDetectors() const {return tab_Detectors;};
  void PrintallDetectorsInfo();
  void PrintLabel2Detnbrs();
  std::map<TString, int> GetDetTypeValues() const {return DetType2nbr;};
  std::map<TString, int> GetDetName2Detnbrs() const {return DetName2nbr;};
  //std::vector<int>     GetLabel2Detnbrs() const {return Label2Detnbr;};
  int    GetLabel2Detnbrs(const int & i) const {return Label2Detnbr.at(i);};
  void SetReferenceDetector(const Int_t &label);
  CDetectors* GetReferenceDetector() const {return Ref_detector;};
  CDetectors* GetDetector(Int_t i) const {return tab_Detectors.at(i);};
  Int_t LoadCalibration(TString CalibrationFileName);
  Int_t LoadTimeCalibration(TString CalibrationFileName);


// Function to set the data I need to analyse
public:
  std::vector<TChain *> GettheTChain() const {return chained_oak;};
  Int_t ChainBuilder(Bool_t isMultiChain, Bool_t per_Run);
  Int_t ChainCounter(Bool_t per_Run);

protected:
  void SetDetTypeValues();
  void SetDetName2Detnbrs();
  void SetLabel2Detnbrs();

  //ClassDef(CExperiment,1);
};

#endif
