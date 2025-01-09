//-------------------------------------------------------------------------//
//                                                                         //
//                               CDetectors.h                              //
//                               Version 1.0                               //
//                        Matthieu Lebois November 2021                    //
//                                                                         //
//  This file contain the definition of the function to define meta        //
//  Parameter to explicitely create everything needed to analyse an        //
//  experiment                                                             //
//                                                                         //
//-------------------------------------------------------------------------//

#ifndef CDETECTORS_H
#define CDETECTORS_H 1

#include "TString.h"
#include "TMath.h"
#include <map>
#include <string>
#include <iostream>
#include "colormod.h"
#include "TCutG.h"

// Definition of the I/O function
using namespace std;

class CDetectors{
public:
  CDetectors( TString name , Int_t label );
  CDetectors();
  virtual ~CDetectors();

public:
  TString Detname;
  Int_t   Detlabel;
  TString DetCARDType;
  TString DetType;
  Int_t ring;
  Int_t alveole;
  Int_t PID;
  Float_t x3,x2,x1,x0;
  Float_t x3prime,x2prime,x1prime,x0prime;
  Int_t FASTERMaxChannelnbr;
  Int_t maxch;
  Int_t nbrchannels;
  Int_t chargenbr;
  Long64_t tm_shift;
  Float_t period;
  Bool_t isReferenceDetector;
  Bool_t isbgo;
  Bool_t isprompt;
  TCutG *mycut;
  Float_t angle;
  Float_t tan_angle;
  Float_t labr_discri_pos;
  Float_t labr_discri_sigma;
  Bool_t  ispurelabr;
  Float_t nai_discri_pos;
  Float_t nai_discri_sigma;
  Bool_t  ispurenai;

public:
  void SetDetectorName ( TString name ){Detname = name;};
  void SetDetectorCARDType ( TString name ){DetCARDType = name;};
  void SetDetectorType ( TString name ){DetType = name;};
  void SetDetectorlabel( Int_t nbr)    {Detlabel=nbr;};
  void SetIsReferenceDetector(Bool_t isReference) {isReferenceDetector = isReference;};
  void PrintDetInfo();
  void SetCalib( Float_t a1, Float_t b1, Float_t c1, Float_t d1){x3=a1;x2=b1;x1=c1;x0=d1;};
  void SetSecondaryCalib( Float_t a1, Float_t b1, Float_t c1, Float_t d1){x3prime=a1;x2prime=b1;x1prime=c1;x0prime=d1;};
  void SetTimeShift( Double_t time){tm_shift = time;};
  void SetRing(int number){ring = number;};
  void SetAlveole(int number){alveole = number;};
  void SetPID(int number){PID = number;};
  Bool_t IsPrompt(TCutG, Double_t , Double_t);

public:
  TString GetDetectorName (){return Detname;};
  TString GetDetectorCARDType (){return DetType;};
  TString GetDetectorType (){return DetType;};
  Int_t GetDetectorlabel  (){return Detlabel;};
  Bool_t IsReferenceDetector(){return isReferenceDetector;};
  Int_t  GetRing() {return ring;};
  Int_t  GetAlveole() {return alveole;};
  Int_t  GetPID() {return PID;};
  Float_t GetCalibx3() {return x3;};
  Float_t GetCalibx2() {return x2;};
  Float_t GetCalibx1() {return x1;};
  Float_t GetCalibx0() {return x0;};
  Float_t GetSecondCalibx3() {return x3prime;};
  Float_t GetSecondCalibx2() {return x2prime;};
  Float_t GetSecondCalibx1() {return x1prime;};
  Float_t GetSecondCalibx0() {return x0prime;};
  Double_t GetTimeShift() {return tm_shift;};
  Double_t GetEnergy(Double_t NRJ);
  Double_t GetSecondaryEnergy(Double_t NRJ);

// Virtual functions to access daughter class informations
public:
  virtual Int_t    GetMaxchNumber() {return maxch;};
  virtual Int_t    GetNbrChannels() {return nbrchannels;};
  virtual Int_t    GetFASTERNbrChannels() {return FASTERMaxChannelnbr;}
  virtual Int_t    GetChargeNbr() {return chargenbr;};
  virtual Float_t  GetRotationAngle(){return angle;};
  virtual Float_t  GetRotationAngleTan(){return tan_angle;};
  virtual Float_t  GetLaBrDiscriPosition(){return labr_discri_pos;};
  virtual Float_t  GetLaBrDiscriSigma(){return labr_discri_sigma;};
  virtual Float_t  GetNaIDiscriPosition(){return nai_discri_pos;};
  virtual Float_t  GetNaIDiscriSigma(){return nai_discri_sigma;};
  virtual Bool_t   IsBGO() {return isbgo;};
  virtual Bool_t   IsPureLaBr3(Double_t PSD);
  virtual Bool_t   IsPureNaI(Double_t PSD);
  virtual Bool_t   IsBeyondLaBr3andNaI(Double_t PSD);
  virtual void     SetChargeNbr( Int_t charge) {chargenbr = charge;};
  virtual void     SetIsBGO( Bool_t isBGO) {isbgo = isBGO;} ;
  virtual void     SetPeriod(Double_t Period) {period=Period;};
  virtual void     SetRotationAngle(Float_t angle1){angle=angle1;};
  virtual void     SetRotationAngleTan(Float_t angle1){tan_angle=angle1;};
  virtual void     SetLaBrDiscriPosition(Float_t pos){labr_discri_pos=pos;};
  virtual void     SetLaBrDiscriSigma(Float_t sigma){labr_discri_sigma=sigma;};
  virtual void     SetNaIDiscriPosition(Float_t pos){nai_discri_pos=pos;};
  virtual void     SetNaIDiscriSigma(Float_t sigma){nai_discri_sigma=sigma;};
  virtual Double_t GetPeriod() {return period;};
};

#endif
