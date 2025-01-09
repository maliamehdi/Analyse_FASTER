// Call for the principal parts of the functions
#include "include/CDetectors.h"


CDetectors::CDetectors(TString name , Int_t  label):
Detname("UNKNW"),
Detlabel(0), DetCARDType("UNKNW"), DetType("UNKNW"),
ring(0),alveole(0),PID(0),x3(0.),x2(0.),x1(1.),x0(0.),x3prime(0.),x2prime(0.),x1prime(1.),x0prime(0.),
FASTERMaxChannelnbr(1),maxch(1),nbrchannels(1),chargenbr(1),tm_shift(0.),
isReferenceDetector(kFALSE),isbgo(kFALSE),isprompt(kFALSE),
angle(0.),tan_angle(0.),labr_discri_pos(0.),labr_discri_sigma(0.),ispurelabr(kFALSE),
nai_discri_pos(0.),nai_discri_sigma(0.),ispurenai(kFALSE)
{
  Detname = name;
  Detlabel = label;
  //cout << "Creating Detector " << Detname << " with the DAQ label " << Detlabel << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CDetectors::CDetectors():
Detname("UNKNW"),
Detlabel(0), DetCARDType("UNKNW"), DetType("UNKNW"),
ring(0),alveole(0),PID(0),x3(0.),x2(0.),x1(1.),x0(0.),x3prime(0.),x2prime(0.),x1prime(1.),x0prime(0.),
FASTERMaxChannelnbr(1),maxch(1),nbrchannels(1),chargenbr(1),tm_shift(0.),
isReferenceDetector(kFALSE),isbgo(kFALSE),isprompt(kFALSE),
angle(0.),tan_angle(0.),labr_discri_pos(0.),labr_discri_sigma(0.),ispurelabr(kFALSE),
nai_discri_pos(0.),nai_discri_sigma(0.),ispurenai(kFALSE)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CDetectors::~CDetectors(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CDetectors::PrintDetInfo()
{
  cout << "Detector " << Detname << " labeled " << Detlabel << endl;
  cout << "Is Coded with a " << DetCARDType << endl;
  cout << "And it is a " << DetType << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t CDetectors::GetEnergy(Double_t NRJ)
{

  Double_t e3 = TMath::Power((Double_t)NRJ,3);
  Double_t e2 = TMath::Power((Double_t)NRJ,2);
  Double_t e1 = TMath::Power((Double_t)NRJ,1);
  Double_t nrj = x3*e3+x2*e2+x1*e1+x0;

  //cout << "Energy before = " << NRJ << " & after " << nrj << endl;
  //cout << "Parameters were : " << x3 << "\t" << x2<< "\t" << x1<< "\t" << x0 << endl;

  return nrj;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t CDetectors::GetSecondaryEnergy(Double_t NRJ)
{
  Double_t e3 = TMath::Power((Double_t)NRJ,3);
  Double_t e2 = TMath::Power((Double_t)NRJ,2);
  Double_t e1 = TMath::Power((Double_t)NRJ,1);
  Double_t nrj = x3prime*e3+x2prime*e2+x1prime*e1+x0prime;

  return nrj;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool CDetectors::IsPureLaBr3(Double_t PSD)
{
  Double_t REF = TMath::Abs(PSD-labr_discri_pos);
  if(REF <= 2*2.35*labr_discri_sigma)  return true;
  else return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


bool CDetectors::IsPureNaI(Double_t PSD)
{
  Double_t REF = TMath::Abs(PSD-nai_discri_pos);
  if(REF <= 2*2.35*nai_discri_sigma)  return true;
  else return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool CDetectors::IsBeyondLaBr3andNaI(Double_t PSD)
{
  Double_t REFLaBr = PSD-labr_discri_pos;
  Double_t REFNaI  = PSD-nai_discri_pos;
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
