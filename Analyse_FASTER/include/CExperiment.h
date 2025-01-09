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
  int nthreads;
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
  Int_t    Nbrofdetectors;
  Int_t    NbrofADC;
  Int_t    NbrofQDC;
  Int_t    NbrofRF;
  Int_t    NbrofSAMPLER;
  Int_t    NbrofModules;
  Int_t    NbrofTChains;
  Bool_t   isMultiChain;
  Int_t    maxLabel;
  Bool_t   isQDC2;
  Bool_t   isQDC3;
  Bool_t   isQDC4;
  Bool_t   isPARISin;

  // Mapping between FASTER label and spectrum number
  std::map<TString, int> DetType2nbr;std::map<TString, int> DetName2nbr;
  std::vector<int> it_Label2Detnbr;
  std::vector<int> Label2Detnbr;
  std::vector<CDetectors*>::iterator it_Detectors;
  std::vector<CDetectors*> tab_Detectors;
  CDetectors *Ref_detector;

  // Data I want to manipulate
  std::vector<TChain *>::iterator it_chained_oak;
  std::vector<TChain *> chained_oak;

public:
  void SetThreadsNbr(int n){nthreads = n;};
  void SetExperienceName( TString name ){Experiencename = name;};
  void SetPIDfilename( TString name ){PIDfilename = name;};
  void SetPIDlistname( TString name ){PIDlistname = name;};
  void SetTypelistname( TString name ){Typelistname = name;};
  void SetRUNlistname( TString name );
  void SetParameterFolder( TString name ){Parameter_folder=name;};
  void SetNRJCalibration_filename( TString name ){NRJCalibration_filename = name;};
  void SetTimeCalibration_filename( TString name ){TimeCalibration_filename = name;};
  void SetTimewindows_filename( TString name ){Timewindows_filename = name;};
  void SetPARISAngle_filename( TString name ){PARISangles_filename = name;};
  void SetFileDirectory_IN( TString name ){FileDirectory_IN = name;};
  void SetFileDirectory_OUT( TString name ){FileDirectory_OUT = name;};
  void SetDataTreeName( std::string name ){DataTreeName = name;};
  Int_t SetNbrofdetectors( TString name );
  void SetNbrofdetectors(Int_t nbr){Nbrofdetectors=nbr;};
  void SetNbrofADC( Int_t nbr ){NbrofADC=nbr;};
  void SetNbrofQDC( Int_t nbr ){NbrofQDC=nbr;};
  void SetNbrofRF( Int_t nbr ){NbrofRF=nbr;};
  void SetNbrofModules( Int_t nbr ){NbrofModules=nbr;};
  void SetNbrofTChains( Int_t nbr ){NbrofTChains=nbr;};
  void SetQDCNbr(int nbr);
  void SetThereIsPARIS(Bool_t isPARIS) {isPARISin = isPARIS;};

// Experiment variables Get functions
public:
  int      GetThreadsNbr(){return nthreads;};
  TString  GetExperienceName(){return Experiencename;};
  TString  GetPIDfilename(){return PIDfilename;};
  TString  GetPIDlistname(){return PIDlistname;};
  TString  GetTypelistname(){return Typelistname;};
  TString  GetRUNlistname(){return Runlist_filename;};
  TString  GetParameterFolder(){return Parameter_folder;};
  TString  GetNRJCalibration_filename(){return NRJCalibration_filename;};
  TString  GetTimeCalibration_filename(){return TimeCalibration_filename;};
  TString  GetPARISAngle_filename(){return PARISangles_filename;};
  TString  GetTimewindows_filename(){return Timewindows_filename;};
  TString  GetFileDirectory_IN(){return FileDirectory_IN;};
  TString  GetFileDirectory_OUT(){return FileDirectory_OUT;};
  std::string  GetDataTreeName(){return DataTreeName;};
  std::vector<std::string> GetOutputFileNames(){return Outputfilenames;};
  std::vector<std::string> GetDataFileNames(){return Datafilenames;};
  Int_t    GetNbrofdetectors(){return Nbrofdetectors;};
  Int_t    GetNbrofADC(){return NbrofADC;};
  Int_t    GetNbrofQDC(){return NbrofQDC;};
  Int_t    GetNbrofRF(){return NbrofRF;};
  Int_t    GetNbrofModules(){return NbrofModules;};
  Int_t    GetNbrofTChains(){return NbrofTChains;};
  Bool_t   GetisQDC2(){return isQDC2;}
  Bool_t   GetisQDC3(){return isQDC3;}
  Bool_t   GetisQDC4(){return isQDC4;}
  Bool_t   IsTherePARIS(){return isPARISin;};

// Detectors variables
public:
  int SetDetectors(  TString inputfilename, TString outputfilename, TString Typefilename );
  std::vector<CDetectors*> GetDetectors(){return tab_Detectors;};
  void PrintallDetectorsInfo();
  void PrintLabel2Detnbrs();
  std::map<TString, int> GetDetTypeValues(){return DetType2nbr;};
  std::map<TString, int> GetDetName2Detnbrs(){return DetName2nbr;};
  //std::vector<int>     GetLabel2Detnbrs(){return Label2Detnbr;};
  int    GetLabel2Detnbrs(const int & i){return Label2Detnbr.at(i);};
  void SetReferenceDetector(Int_t label);
  CDetectors* GetReferenceDetector(){return Ref_detector;};
  CDetectors* GetDetector(Int_t i){return tab_Detectors.at(i);};
  Int_t LoadCalibration(TString CalibrationFileName);
  Int_t LoadTimeCalibration(TString CalibrationFileName);


// Function to set the data I need to analyse
public:
  std::vector<TChain *> GettheTChain(){return chained_oak;};
  Int_t ChainBuilder(Bool_t isMultiChain, Bool_t per_Run);
  Int_t ChainCounter(Bool_t per_Run);

protected:
  void SetDetTypeValues();
  void SetDetName2Detnbrs();
  void SetLabel2Detnbrs();

  //ClassDef(CExperiment,1);
};

#endif
