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
qdcindex(0),detectorposition(0),FASTERMaxChannelnbr(1073741824),maxch(2100000),nbrchannels(8192),chargenbr(1),
theta_angle(0.), discriLaBr_pos(0.), discriLaBr_sigma(0.), discriNaI_pos(0.), discriNaI_sigma(0.)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

QDC::QDC():CDetectors(),
qdcindex(0),detectorposition(0),FASTERMaxChannelnbr(1073741824),maxch(2100000),nbrchannels(8192),chargenbr(1),
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
void QDC::SetResA(Float_t resA)
{this->resA = resA;}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void QDC::SetRespower(Float_t respower)
{this->respower = respower;}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void QDC::SetRawResA(Float_t RawresA)
{this->RawresA = RawresA;}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void QDC::SetRawRespower(Float_t Rawrespower)
{this->Rawrespower = Rawrespower;}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void QDC::SetCaliba(Float_t caliba)
{this->caliba = caliba;}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void QDC::SetCalibb(Float_t calibb)
{this->calibb = calibb;}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void QDC::SetCaliba2(Float_t caliba2)
{this->caliba2 = caliba2;}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void QDC::Print( )
{
  cout << "Print QDC informations : " << endl;
  cout << "qdcindex : " << qdcindex << endl;;
  cout << " detectorposition : " << detectorposition << endl;
  cout << " FASTERMaxChannelnbr : " << FASTERMaxChannelnbr << endl;
  cout << " maxch : " << maxch << endl;
  cout << " nbrchannels  : " << nbrchannels << endl;
  cout << " chargenbr : " << chargenbr << endl;
  cout << " theta_angle : " << theta_angle << endl;
  cout << " theta_angle_tan : " << theta_angle_tan << endl;
  cout << " discriLaBr_pos : " << discriLaBr_pos << endl;
  cout << " discriLaBr_sigma : " << discriLaBr_sigma << endl;
  cout << " discriNaI_pos : " << discriNaI_pos << endl;
  cout << " discriNaI_sigma : " << discriNaI_sigma << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Double_t QDC::CeBrqdcResolution(const double Qs)
{
  double power = -0.5;
  double slope = 13.477;
  double resolution = slope*TMath::Power(Qs,power);
  
  return resolution;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Double_t QDC::SimpleResolution(const double E)
{
  double power = -0.49;
  double slope = 4.68;
  double resolution = slope*TMath::Power(E,power);
  
  return resolution;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Double_t QDC::SimpleResolutionNaI(const double E)
{
  double power = -0.49;
  double slope = 2*4.68;
  double resolution = slope*TMath::Power(E,power);
  
  return resolution;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Bool_t QDC::IsPureLaBr3(const Double_t PSD, const Double_t Qs,const Double_t QL)
{
  
  // First I check the angle
  Bool_t istheta_OK = false;
  Bool_t isQs_inf_O_sup = false;
  Bool_t isQs_sup_O_inf = false;
  double labr_discri_pos = GetLaBrDiscriPosition();
  double labr_discri_sigma = GetLaBrDiscriSigma();
  double QSref = QL*tan(labr_discri_pos);
  double resA = GetRawResA();
  double respower = GetRawRespower();
  Double_t REF = TMath::Abs(PSD-labr_discri_pos);
  //if(REF <= 2.35*labr_discri_sigma)  istheta_OK = true;
  
  // I calculate the resolution for this charge
  double resolution = RawresA*TMath::Power(Qs,Rawrespower);
  //This is only valid for PARIS305
  //double resolution = QDC::CeBrqdcResolution(Qs);
  
  // I define coordinate at the origin for LaBr upper and lower cut lines
  double ordinate_sup = +3000.;
  double ordinate_inf = -1500.;
  
  // I need to calculate the slope of the upper cut line
  double Qsmax = QSref+resolution*QSref;//+ordinate_sup; //QL*tan(labr_discri_pos+2.35*labr_discri_sigma);//QSref+resolution*QSref+ordinate_sup;
  double slopemax = (Qsmax)/QL;
  double thetamax = TMath::ATan(Qsmax/QL);
  
  // I can check that my Qs is below the limit corresponding to its QL
  //if(Qs <= (slopemax*QL+ordinate_sup) ) isQs_inf_O_sup = true;
  //if(Qs <= (Qsmax+ordinate_sup) ) isQs_inf_O_sup = true;
  
  // I do the same for the lower line
  double Qsmin = (QSref-resolution*QSref);//QL*tan(labr_discri_pos-2.35*labr_discri_sigma);+ordinate_inf
  double slopemin = (Qsmin)/QL;
  double thetamin = TMath::ATan(Qsmin/QL);
  //if(Qs >= (slopemin*QL+ordinate_inf) ) isQs_sup_O_inf = true;
  //if(Qs >= (Qsmin+ordinate_inf) ) isQs_sup_O_inf = true;

  if(PSD <= thetamax && PSD >= thetamin && Qs>=Qsmin && Qs<=Qsmax)  istheta_OK = true;
  
  //if (Qs < 10000) {isQs_inf_O_sup = true;isQs_sup_O_inf = true;}
  // Then I check the short charge variation (within the LaBr3 normal resolution 3%)
  // Double_t localQs = QL*TMath::Tan(labr_discri_pos);
  // Bool_t isQs_OK = false;
  // if(TMath::Abs(localQs-Qs)/Qs < 3.e-2) isQs_OK=true;
  // cout << endl << endl;
  // cout << "Detector # = " << qdcindex << endl;
  // cout << "PSD = " << PSD << "; labr_discri_pos = " << labr_discri_pos << " and labr_discri_sigma = " << labr_discri_sigma << endl;
  // cout << "Qs = " << Qs << "; Ql = " << QL << endl;
  // cout << "istheta_OK = " << (int) istheta_OK << " isQs_OK = " << (int) isQs_OK << endl;
  
  //if(isQs_inf_O_sup && isQs_sup_O_inf && istheta_OK) return true;
  if(istheta_OK) return true;
  else return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


Bool_t QDC::IsPureNaI(const Double_t PSD, const Double_t Qs,const Double_t QL)
{
  
  // First I check the angle
  Bool_t istheta_OK = false;
  // Bool_t isQs_inf_O_sup = false;
  // Bool_t isQs_sup_O_inf = false;
  double nai_discri_pos = GetNaIDiscriPosition();
  double nai_discri_sigma = GetNaIDiscriSigma();
  double QSref = QL*tan(nai_discri_pos);
  
  //cout<<"Qsref NaI : "<<QSref<<endl;

  Double_t REF = TMath::Abs(PSD-nai_discri_pos);
  //if(REF <= 2.35*labr_discri_sigma)  istheta_OK = true;
  
  // I calculate the resolution for this charge
  double resolution = GetDetectorResolution(Qs)*2;
  
  // I define coordinate at the origin for NaI upper and lower cut lines
  double ordinate_sup = 2000.;
  double ordinate_inf = -200.;
  
  // I need to calculate the slope of the upper cut line
  double Qsmax = QL*tan(nai_discri_pos + 2.35 * nai_discri_sigma)+ordinate_sup;//QSref+resolution*QSref+ordinate_sup;

  //cout<<"Qsmax NaI : "<<Qsmax<<endl;
  //double slopemax = (Qsmax)/QL;
  double thetamax = TMath::ATan(Qsmax/QL);
  
  // I can check that my Qs is below the limit corresponding to its QL
  //if(Qs <= (slopemax*QL+ordinate_sup) ) isQs_inf_O_sup = true;
  //if(Qs <= (Qsmax+ordinate_sup) ) isQs_inf_O_sup = true;
  
  // I do the same for the lower line
  double Qsmin = QL*tan(nai_discri_pos - 2.35 * nai_discri_sigma)+ordinate_inf; //(QSref-resolution*QSref)+ordinate_inf;
  //cout<<"Qsmin NaI : "<<Qsmin<<endl;
  //double slopemin = (Qsmin)/QL;
  double thetamin = TMath::ATan(Qsmin/QL);
  //if(Qs >= (slopemin*QL+ordinate_inf) ) isQs_sup_O_inf = true;
  //if(Qs >= (Qsmin+ordinate_inf) ) isQs_sup_O_inf = true;
  if(PSD <= thetamax && PSD >= thetamin)  istheta_OK = true;
  //if (Qs < 10000) {isQs_inf_O_sup = true;isQs_sup_O_inf = true;}
  // Then I check the short charge variation (within the LaBr3 normal resolution 3%)
  // Double_t localQs = QL*TMath::Tan(labr_discri_pos);
  // Bool_t isQs_OK = false;
  // if(TMath::Abs(localQs-Qs)/Qs < 3.e-2) isQs_OK=true;
  // cout << endl << endl;
  // cout << "Detector # = " << qdcindex << endl;
  // cout << "PSD = " << PSD << "; labr_discri_pos = " << labr_discri_pos << " and labr_discri_sigma = " << labr_discri_sigma << endl;
  // cout << "Qs = " << Qs << "; Ql = " << QL << endl;
  // cout << "istheta_OK = " << (int) istheta_OK << " isQs_OK = " << (int) isQs_OK << endl;
  
  //if(isQs_inf_O_sup && isQs_sup_O_inf && istheta_OK) return true;
  if(istheta_OK) return true;
  else return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Bool_t QDC::IsBeyondLaBr3andNaI(Double_t PSD, const Double_t Qs,const Double_t QL)
{
  double labr_discri_pos = GetLaBrDiscriPosition();
  double nai_discri_pos = GetNaIDiscriPosition();
  //Double_t REFLaBr = PSD-labr_discri_pos;
  //Double_t REFNaI  = PSD-nai_discri_pos;
  Bool_t isNaI = IsPureNaI(PSD,Qs,QL);
  Bool_t isLaBr = IsPureLaBr3(PSD,Qs,QL);
  Bool_t outofLaBr = kFALSE;
  Bool_t outofNaI = kFALSE;
  if(PSD > labr_discri_pos && !isLaBr) outofLaBr = kTRUE;
  if(PSD < nai_discri_pos && !isNaI) outofNaI = kTRUE;
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
