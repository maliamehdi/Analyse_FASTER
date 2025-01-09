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
  Int_t nbrchannels = 65536;
  Int_t detectorposition;
  Bool_t isbgo=false;

public:
  void SetADCIndex( Int_t );
  void SetADCPosition( Int_t );
  void SetIsBGO( Bool_t isBGO);

public:
  Int_t    GetADCIndex(){return adcindex;}
  Int_t    GetADCPosition(){return detectorposition;}
  using    CDetectors::GetMaxchNumber;
  Int_t    GetMaxchNumber() {return maxch;}
  using    CDetectors::GetNbrChannels;
  Int_t    GetNbrChannels() {return nbrchannels;}
  using    CDetectors::GetFASTERNbrChannels;
  Int_t    GetFASTERNbrChannels() {return FASTERMaxChannelnbr;}
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
  Int_t nbrchannels = 16384;
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

  // The Get Functions
  Int_t   GetQDCIndex(){return qdcindex;}
  Int_t   GetQDCPosition(){return detectorposition;}
  using   CDetectors::GetMaxchNumber;
  Int_t   GetMaxchNumber() {return maxch;}
  using   CDetectors::GetNbrChannels;
  Int_t   GetNbrChannels() {return nbrchannels;}
  using   CDetectors::GetFASTERNbrChannels;
  Int_t   GetFASTERNbrChannels() {return FASTERMaxChannelnbr;}
  using   CDetectors::GetChargeNbr;
  Int_t   GetChargeNbr(){return chargenbr;}
  using   CDetectors::GetRotationAngle;
  Float_t GetRotationAngle(){return theta_angle;};
  using   CDetectors::GetRotationAngleTan;
  Float_t GetRotationAngleTan(){return theta_angle_tan;};
  using   CDetectors::GetLaBrDiscriPosition;
  Float_t GetLaBrDiscriPosition(){return discriLaBr_pos;};
  using   CDetectors::GetLaBrDiscriSigma;
  Float_t GetLaBrDiscriSigma(){return discriLaBr_sigma;};
  using   CDetectors::GetNaIDiscriPosition;
  Float_t GetNaIDiscriPosition(){return discriNaI_pos;};
  using   CDetectors::GetNaIDiscriSigma;
  Float_t GetNaIDiscriSigma(){return discriNaI_sigma;};

public:
  using   CDetectors::IsPureLaBr3;
  bool IsPureLaBr3(Double_t PSD);
  using   CDetectors::IsPureNaI;
  bool IsPureNaI(Double_t PSD);
  using   CDetectors::IsBeyondLaBr3andNaI;
  bool IsBeyondLaBr3andNaI(Double_t PSD);

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
  Double_t GetPeriod(){return period;}

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

  Int_t   GetSAMPLERIndex(){return samplerindex;}
  Int_t   GetSAMPLERPosition(){return detectorposition;}
  using   CDetectors::GetMaxchNumber;
  Int_t   GetMaxchNumber() {return maxch;}
  using   CDetectors::GetNbrChannels;
  Int_t   GetNbrChannels() {return nbrchannels;}
  using   CDetectors::GetFASTERNbrChannels;
  Int_t   GetFASTERNbrChannels() {return FASTERMaxChannelnbr;}
};

#endif //TDETECTEUR_H
