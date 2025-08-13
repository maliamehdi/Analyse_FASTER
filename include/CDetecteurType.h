#ifndef CDETECTEURTYPE_H
#define CDETECTEURTYPE_H 1

//-------------------------------------------------------------------------//
//                                                                         //
//                             TDetecteurType.h                            //
//                               Version 1.0                               //
//                        Matthieu Lebois November 2021                    //
//                                                                         //
//  This file contain the definition of the various types of detectors     //                                        //
//                                                                         //
//-------------------------------------------------------------------------//
#include "CDetectors.h"
#include "colormod.h"

// Call for libraries needed locally
#include "TApplication.h"
#include "TFile.h"
#include "TCutG.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TObject.h"
#include "TTree.h"
#include "TBasket.h"
#include "TString.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TChain.h"


class CExperiment;

class ADC: public CDetectors
{

public:
  // Constructors
  ADC(TString name, Int_t label);
  ADC();
  virtual ~ADC();

public:
  Int_t adcindex;
  Int_t FASTERMaxChannelnbr;
  Int_t maxch = 2500000;
  Int_t nbrchannels = static_cast<Int_t>(TMath::Power(2,22));//65536;
  Int_t detectorposition;
  Bool_t isbgo=false;

public:
  void SetADCIndex( Int_t );
  void SetADCPosition( Int_t );
  void SetIsBGO( Bool_t isBGO);

public:
  Int_t    GetADCIndex() const {return adcindex;}
  Int_t    GetADCPosition() const {return detectorposition;}
  using    CDetectors::GetMaxchNumber;
  Int_t    GetMaxchNumber() const {return maxch;}
  using    CDetectors::GetNbrChannels;
  Int_t    GetNbrChannels() const {return nbrchannels;}
  using    CDetectors::GetFASTERNbrChannels;
  Int_t    GetFASTERNbrChannels() const {return FASTERMaxChannelnbr;}
  using    CDetectors::IsBGO;
  Bool_t   IsBGO() const {return isbgo;}
};

class QDC: public CDetectors
{
public:
  QDC(TString name, Int_t label);
  QDC();
  ~QDC();

public:
  Int_t qdcindex;
  Int_t detectorposition;
  Int_t FASTERMaxChannelnbr;
  Int_t maxch = 2390;
  Int_t nbrchannels = static_cast<Int_t>(TMath::Power(2,20));//16384;
  Int_t chargenbr= 4;
  Float_t theta_angle = 0.;
  Float_t theta_angle_tan = 0.;
  Float_t discriLaBr_pos=0.;
  Float_t discriLaBr_sigma=0.;
  Float_t discriNaI_pos=0.;
  Float_t discriNaI_sigma=0.;

public:
  // The Set Functions
  void SetQDCIndex( Int_t );
  void SetQDCPosition( Int_t );
  void SetChargeNbr( Int_t charge);
  void SetRotationAngle(Float_t angle);
  void SetRotationAngleTan(Float_t angle);
  void SetLaBrDiscriPosition(Float_t pos);
  void SetLaBrDiscriSigma(Float_t sigma);
  void SetNaIDiscriPosition(Float_t pos);
  void SetNaIDiscriSigma(Float_t sigma);
  void SetResA(Float_t resA);
  void SetRespower(Float_t respower);
  void SetCaliba(Float_t caliba);
  void SetCalibb(Float_t calibb);
  void SetCaliba2(Float_t caliba2);
  void Print();

  // The Get Functions
  Int_t   GetQDCIndex() const {return qdcindex;}
  Int_t   GetQDCPosition() const {return detectorposition;}
  using   CDetectors::GetMaxchNumber;
  Int_t   GetMaxchNumber() const {return maxch;}
  using   CDetectors::GetNbrChannels;
  Int_t   GetNbrChannels() const {return nbrchannels;}
  using   CDetectors::GetFASTERNbrChannels;
  Int_t   GetFASTERNbrChannels() const {return FASTERMaxChannelnbr;}
  using   CDetectors::GetChargeNbr;
  Int_t   GetChargeNbr() const {return chargenbr;}
  using   CDetectors::GetRotationAngle;
  Float_t GetRotationAngle() const {return theta_angle;};
  using   CDetectors::GetRotationAngleTan;
  Float_t GetRotationAngleTan() const {return theta_angle_tan;};
  using   CDetectors::GetLaBrDiscriPosition;
  Float_t GetLaBrDiscriPosition() const {return discriLaBr_pos;};
  using   CDetectors::GetLaBrDiscriSigma;
  Float_t GetLaBrDiscriSigma() const {return discriLaBr_sigma;};
  using   CDetectors::GetNaIDiscriPosition;
  Float_t GetNaIDiscriPosition() const {return discriNaI_pos;};
  using   CDetectors::GetNaIDiscriSigma;
  Float_t GetNaIDiscriSigma() const {return discriNaI_sigma;};
  using   CDetectors::GetResA;
  Float_t GetResA() const {return resA;};
  using   CDetectors::GetRespower;
  Float_t GetRespower() const {return respower;};
  using CDetectors::GetCaliba;
  Float_t GetCaliba() const {return caliba;};
  using CDetectors::GetCalibb;
  Float_t GetCalibb() const {return calibb;};
  using CDetectors::GetCaliba2;
  Float_t GetCaliba2() const {return caliba2;};


public:
  using   CDetectors::IsPureLaBr3;
  Bool_t IsPureLaBr3(const Double_t PSD, const Double_t Qs,const Double_t QL);
  using   CDetectors::IsPureNaI;
  Bool_t IsPureNaI(const Double_t PSD, const Double_t Qs,const Double_t QL);
  using   CDetectors::IsBeyondLaBr3andNaI;
  Bool_t IsBeyondLaBr3andNaI(const Double_t PSD, const Double_t Qs,const Double_t QL);
  Double_t SimpleResolution(const double E);
  Double_t SimpleResolutionNaI(const double E);
  Double_t CeBrqdcResolution(const double Qs);

};

class RF: public CDetectors
{
public:
  RF(TString name, Int_t label);
  RF();
  ~RF();

public:
  Double_t period=400.; //in ns

public:
  using  CDetectors::SetPeriod;
  void   SetPeriod( Double_t);

public:
  using    CDetectors::GetPeriod;
  Double_t GetPeriod() const {return period;}

};

class SAMPLER: public CDetectors
{
public:
  SAMPLER(TString name, Int_t label);
  SAMPLER();
  ~SAMPLER();

public:
  Int_t samplerindex;
  Int_t detectorposition;
  Int_t FASTERMaxChannelnbr;
  Int_t maxch = 2390;
  Int_t nbrchannels = 16384;

public:
  void SetSAMPLERIndex( Int_t );
  void SetSAMPLERPosition( Int_t );

  Int_t   GetSAMPLERIndex() const {return samplerindex;}
  Int_t   GetSAMPLERPosition() const {return detectorposition;}
  using   CDetectors::GetMaxchNumber;
  Int_t   GetMaxchNumber() const {return maxch;}
  using   CDetectors::GetNbrChannels;
  Int_t   GetNbrChannels() const {return nbrchannels;}
  using   CDetectors::GetFASTERNbrChannels;
  Int_t   GetFASTERNbrChannels() const {return FASTERMaxChannelnbr;}
};

#endif //TDETECTEUR_H
