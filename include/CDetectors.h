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

protected:
  TString Detname;
  Int_t   Detlabel;
  TString DetCARDType;
  TString DetType;
  Int_t ring;
  Int_t alveole;
  Int_t PID;
  Float_t x3,x2,x1,x0;
  Float_t x3prime,x2prime,x1prime,x0prime;
  Double_t res_cst,res_amp,res_pow;
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
  Float_t resA;
  Float_t respower;
  Float_t caliba;
  Float_t calibb;
  Bool_t  ispurenai;

public:
  void SetDetectorName (const TString& name ){Detname = name;};
  void SetDetectorCARDType ( const TString&  name ){DetCARDType = name;};
  void SetDetectorType ( const TString&  name ){DetType = name;};
  void SetDetectorlabel( Int_t nbr)    {Detlabel=nbr;};
  void SetIsReferenceDetector(Bool_t isReference) {isReferenceDetector = isReference;};
  void PrintDetInfo();
  void SetCalib( Float_t a1, Float_t b1, Float_t c1, Float_t d1){x3=a1;x2=b1;x1=c1;x0=d1;};
  void SetSecondaryCalib( Float_t a1, Float_t b1, Float_t c1, Float_t d1){x3prime=a1;x2prime=b1;x1prime=c1;x0prime=d1;};
  void SetResolution(Float_t a1, Float_t b1, Float_t c1){res_cst=a1;res_amp=b1;res_pow=c1;};
  void SetTimeShift( Double_t time){tm_shift = time;};
  void SetRing(int number){ring = number;};
  void SetAlveole(int number){alveole = number;};
  void SetPID(int number){PID = number;};
  Bool_t IsPrompt(TCutG, Double_t , Double_t);  

public:
  const TString& GetDetectorName () const {return Detname;};
  const TString&  GetDetectorCARDType () const {return DetType;};
  const TString&  GetDetectorType () const {return DetType;};
  Int_t GetDetectorlabel  () const {return Detlabel;};
  Bool_t IsReferenceDetector() const {return isReferenceDetector;};
  Int_t  GetRing() const {return ring;};
  Int_t  GetAlveole() const {return alveole;};
  Int_t  GetPID() const {return PID;};
  Float_t GetCalibx3() const {return x3;};
  Float_t GetCalibx2() const {return x2;};
  Float_t GetCalibx1() const {return x1;};
  Float_t GetCalibx0() const {return x0;};
  Float_t GetSecondCalibx3() const {return x3prime;};
  Float_t GetSecondCalibx2() const {return x2prime;};
  Float_t GetSecondCalibx1() const {return x1prime;};
  Float_t GetSecondCalibx0() const {return x0prime;};
  Float_t GetResolution_cst() const {return res_cst;};
  Float_t GetResolution_amp() const {return res_amp;};
  Float_t GetResolution_pow() const {return res_pow;};
  Float_t GetResA() const {return resA;};
  Float_t GetRespower() const {return respower;};
  Float_t GetCaliba() const {return caliba;};
  Float_t GetCalibb() const {return calibb;};
  Double_t GetTimeShift() const {return tm_shift;};
  Double_t GetEnergy(Double_t NRJ);
  Double_t GetSecondaryEnergy(Double_t NRJ);
  Double_t GetDetectorResolution(Double_t NRJ);

// Virtual functions to access daughter class informations
public:
  virtual Int_t    GetMaxchNumber() const {return maxch;};
  virtual Int_t    GetNbrChannels() const {return nbrchannels;};
  virtual Int_t    GetFASTERNbrChannels() const {return FASTERMaxChannelnbr;}
  virtual Int_t    GetChargeNbr() const {return chargenbr;};
  virtual Float_t  GetRotationAngle() const {return angle;};
  virtual Float_t  GetRotationAngleTan() const {return tan_angle;};
  virtual Float_t  GetLaBrDiscriPosition() const {return labr_discri_pos;};
  virtual Float_t  GetLaBrDiscriSigma() const {return labr_discri_sigma;};
  virtual Float_t  GetNaIDiscriPosition() const {return nai_discri_pos;};
  virtual Float_t  GetNaIDiscriSigma() const {return nai_discri_sigma;};
  virtual Bool_t   IsBGO(){return isbgo;};
  //virtual Bool_t   IsPureLaBr3(Double_t PSD);
  virtual Bool_t   IsPureLaBr3(const Double_t PSD, const Double_t Qs,const Double_t QL);
  //virtual Bool_t   IsPureNaI(Double_t PSD);
  virtual Bool_t   IsPureNaI(const Double_t PSD, const Double_t Qs,const Double_t QL);
  //virtual Bool_t   IsBeyondLaBr3andNaI(Double_t PSD);
  virtual Bool_t   IsBeyondLaBr3andNaI(Double_t PSD, const Double_t Qs,const Double_t QL);
  virtual void     SetChargeNbr( Int_t charge) {chargenbr = charge;};
  virtual void     SetIsBGO( Bool_t isBGO) {isbgo = isBGO;} ;
  virtual void     SetPeriod(Double_t Period) {period=Period;};
  virtual void     SetRotationAngle(Float_t angle1){angle=angle1;};
  virtual void     SetRotationAngleTan(Float_t angle1){tan_angle=angle1;};
  virtual void     SetLaBrDiscriPosition(Float_t pos){labr_discri_pos=pos;};
  virtual void     SetLaBrDiscriSigma(Float_t sigma){labr_discri_sigma=sigma;};
  virtual void     SetNaIDiscriPosition(Float_t pos){nai_discri_pos=pos;};
  virtual void     SetNaIDiscriSigma(Float_t sigma){nai_discri_sigma=sigma;};
  virtual Double_t GetPeriod() const {return period;};
};

#endif
