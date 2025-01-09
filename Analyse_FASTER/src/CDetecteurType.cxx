//-------------------------------------------------------------------------//
//                                                                         //
//                               TAnalyse.cxx                              //
//                               Version 1.2                               //
//                        Matthieu Lebois January 2006                     //
//                            with I. Hrivnacova                           //
//                                                                         //
//  This file contain the definition of the function for a root file       //
//  and the declaration to extract some physical results from the data     //
//  contained in the root file                                             //
//                                                                         //
//-------------------------------------------------------------------------//

// Call for the principal parts of the functions
#include "include/CDetecteurType.h"
#include "include/CDetectors.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ADC::ADC(TString name, Int_t label): CDetectors(name,label),
adcindex(0),FASTERMaxChannelnbr(2097152),maxch(560000),nbrchannels(56000), detectorposition(0),
isbgo(kFALSE)
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ADC::ADC(): CDetectors(),
adcindex(0),FASTERMaxChannelnbr(2097152),maxch(560000),nbrchannels(56000), detectorposition(0),
isbgo(kFALSE)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ADC::~ADC(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ADC::SetADCIndex( Int_t index){adcindex = index;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ADC::SetADCPosition( Int_t channeltag){detectorposition = channeltag;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ADC::SetIsBGO( Bool_t isBGO){isbgo = isBGO;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

QDC::QDC(TString name, Int_t label):CDetectors(name,label),
qdcindex(0),detectorposition(0),FASTERMaxChannelnbr(1073741824),maxch(2100000),nbrchannels(4096),chargenbr(1),
theta_angle(0.), discriLaBr_pos(0.), discriLaBr_sigma(0.), discriNaI_pos(0.), discriNaI_sigma(0.)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

QDC::QDC():CDetectors(),
qdcindex(0),detectorposition(0),FASTERMaxChannelnbr(1073741824),maxch(2100000),nbrchannels(4096),chargenbr(1),
theta_angle(0.), discriLaBr_pos(0.), discriLaBr_sigma(0.), discriNaI_pos(0.), discriNaI_sigma(0.)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

QDC::~QDC(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void QDC::SetQDCIndex( Int_t index){qdcindex = index;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void QDC::SetQDCPosition( Int_t channeltag){detectorposition = channeltag;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void QDC::SetChargeNbr( Int_t charge)
{chargenbr = charge;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void QDC::SetRotationAngle( Float_t angle)
{theta_angle = angle;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void QDC::SetRotationAngleTan( Float_t angle)
{theta_angle_tan = angle;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void QDC::SetLaBrDiscriPosition( Float_t pos)
{discriLaBr_pos = pos;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void QDC::SetLaBrDiscriSigma( Float_t sigma)
{discriLaBr_sigma = sigma;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void QDC::SetNaIDiscriPosition( Float_t pos)
{discriNaI_pos = pos;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void QDC::SetNaIDiscriSigma( Float_t sigma)
{discriNaI_sigma = sigma;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool QDC::IsPureLaBr3(Double_t PSD)
{
  Double_t REF = TMath::Abs(PSD-discriLaBr_pos);
  if(REF <= 2*2.35*discriLaBr_sigma)  return true;
  else return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool QDC::IsPureNaI(Double_t PSD)
{
  Double_t REF = TMath::Abs(PSD-discriNaI_pos);
  if(REF <= 2*2.35*discriNaI_sigma)  return true;
  else return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool QDC::IsBeyondLaBr3andNaI(Double_t PSD)
{
  Double_t REFLaBr = PSD-discriLaBr_pos;
  Double_t REFNaI  = PSD-discriNaI_pos;
  Bool_t isNaI = IsPureNaI(PSD);
  Bool_t isLaBr = IsPureLaBr3(PSD);
  Bool_t outofLaBr = kFALSE;
  Bool_t outofNaI = kFALSE;
  if(REFLaBr < 0 && !isLaBr) outofLaBr = kTRUE;
  if(REFNaI > 0 && !isNaI) outofNaI = kTRUE;
  if(outofLaBr || outofNaI)  return true;
  else return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SAMPLER::SAMPLER(TString name, Int_t label):CDetectors(name,label),
samplerindex(0),detectorposition(0),FASTERMaxChannelnbr(1073741824),maxch(4800000),nbrchannels(4096)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SAMPLER::SAMPLER():CDetectors(),
samplerindex(0),detectorposition(0),FASTERMaxChannelnbr(1073741824),maxch(4800000),nbrchannels(4096)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SAMPLER::~SAMPLER(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SAMPLER::SetSAMPLERIndex( Int_t index){samplerindex = index;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SAMPLER::SetSAMPLERPosition( Int_t channeltag){detectorposition = channeltag;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RF::RF(TString name, Int_t label): CDetectors(name,label),
period(400)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RF::RF(): CDetectors(),
period(400)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RF::~RF(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RF::SetPeriod( Double_t a_in){period=a_in;};
