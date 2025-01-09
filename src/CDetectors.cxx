// Call for the principal parts of the functions
#include "include/CDetectors.h"


CDetectors::CDetectors(TString name , Int_t  label):
Detname("UNKNW"),
Detlabel(0), DetCARDType("UNKNW"), DetType("UNKNW"),
ring(0),alveole(0),PID(0),x3(0.),x2(0.),x1(1.),x0(0.),x3prime(0.),x2prime(0.),x1prime(1.),x0prime(0.),res_cst(0.),res_amp(1),res_pow(0.),
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
ring(0),alveole(0),PID(0),x3(0.),x2(0.),x1(1.),x0(0.),x3prime(0.),x2prime(0.),x1prime(1.),x0prime(0.),res_cst(0.),res_amp(1),res_pow(0.),
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
  // cout << " Detname : " << Detname << endl;
  // cout << " Detlabel : " << Detlabel << endl;
  // cout << " DetCARDType: " << DetCARDType << endl;
  // cout << " DetType: " << DetType << endl;
  // cout << " ring: " << ring << endl;
  // cout << " alveole: " << alveole << endl;
  // cout << " PID: " << PID << endl;
  // cout << " x3,x2,x1,x0: " << x3<<","<<x2<<","<<x1<<","<< x0 << endl;
  // cout << " x3prime,x2prime,x1prime,x0prime: " << x3prime <<","<<x2prime<<","<<x1prime<<","<< x0prime << endl;
  // cout << " FASTERMaxChannelnbr: " << FASTERMaxChannelnbr << endl;
  // cout << " maxch: " << maxch << endl;
  // cout << " nbrchannels: " << nbrchannels << endl;
  // cout << " chargenbr: " << chargenbr << endl;
  // cout << " tm_shift: " << tm_shift << endl;
  // cout << " period: " << period << endl;
  // cout << " isReferenceDetector: " << isReferenceDetector << endl;
  // cout << " isbgo: " << isbgo << endl;
  // cout << " isprompt: " << isprompt << endl;
  // cout << " angle: " << angle << endl;
  // cout << " tan_angle: " << tan_angle << endl;
  // cout << " labr_discri_pos: " << labr_discri_pos << endl;
  // cout << " labr_discri_sigma: " << labr_discri_sigma << endl;
  // cout << "  ispurelabr: " << Detlabel << endl;
  // cout << " nai_discri_pos: " << nai_discri_pos << endl;
  // cout << " nai_discri_sigma: " << nai_discri_sigma << endl;
  // cout << "  ispurenai: " << ispurenai << endl;
  // cout << endl << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t CDetectors::GetEnergy(Double_t NRJ)
{

  //Double_t e3 = TMath::Power((Double_t)NRJ,3);
  Double_t e2 = TMath::Power((Double_t)NRJ,2);
  Double_t e1 = TMath::Power((Double_t)NRJ,1);
  Double_t nrj = x2*e2+x1*e1+x0;//x3*e3+x2*e2+x1*e1+x0;

  // cout << "Energy before = " << NRJ << " & after " << nrj << endl;
  // cout << "Parameters were : " << x3 << "\t" << x2<< "\t" << x1<< "\t" << x0 << endl;

  return nrj;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t CDetectors::GetSecondaryEnergy(Double_t NRJ)
{
  //Double_t e3 = TMath::Power((Double_t)NRJ,3);
  Double_t e2 = TMath::Power((Double_t)NRJ,2);
  Double_t e1 = TMath::Power((Double_t)NRJ,1);
  Double_t nrj = x2prime*e2+x1prime*e1+x0prime;//x3prime*e3+x2prime*e2+x1prime*e1+x0prime;

  return nrj;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t CDetectors::GetDetectorResolution(Double_t NRJ)
{
  Double_t reso = res_cst+res_amp*TMath::Power(NRJ,res_pow);
  return reso;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Bool_t CDetectors::IsPureLaBr3(const Double_t PSD, const Double_t Qs,const Double_t QL)
{
  // First I check the angle
  Bool_t istheta_OK = false;
  double labr_discri_pos = GetLaBrDiscriPosition();
  double labr_discri_sigma = GetLaBrDiscriSigma();
  double resA = GetResA();
  double respower = GetRespower();
  Double_t REF = TMath::Abs(PSD-labr_discri_pos);
  if(REF <= 2*labr_discri_sigma)  istheta_OK = true;
  
  // Then I check the short charge variation (within the LaBr3 normal resolution 3%)
  Double_t localQs = QL*TMath::Tan(PSD);
  Bool_t isQs_OK = false;
  if(TMath::Abs(localQs-Qs)/Qs < 3.e-2) isQs_OK=true;
  // cout << endl << endl;
  // cout << "PSD = " << PSD << "; labr_discri_pos = " << labr_discri_pos << " and labr_discri_sigma = " << labr_discri_sigma << endl;
  // cout << "Qs = " << Qs << "; Ql = " << QL << endl;
  // cout << "istheta_OK = " << (int) istheta_OK << " isQs_OK = " << (int) isQs_OK << endl;
  
  if(isQs_OK && istheta_OK) return true;
  else return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


Bool_t CDetectors::IsPureNaI(const Double_t PSD, const Double_t Qs,const Double_t QL)
{
// First I check the angle
  Bool_t istheta_OK = false;
  double nai_discri_pos = GetNaIDiscriPosition();
  double nai_discri_sigma = GetNaIDiscriSigma();
  Double_t REF = TMath::Abs(PSD-nai_discri_pos);
  if(REF <= 2*nai_discri_sigma)  istheta_OK = true;
  
  // Then I check the short charge variation (within the LaBr3 normal resolution 3%)
  Double_t localQs = QL*TMath::Tan(PSD);
  Bool_t isQs_OK = false;
  if(TMath::Abs(localQs-Qs)/Qs < 7.e-2) isQs_OK=true;
  
  if(isQs_OK && istheta_OK) return true;
  else return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Bool_t CDetectors::IsBeyondLaBr3andNaI(Double_t PSD, const Double_t Qs,const Double_t QL)
{
  double labr_discri_pos = GetLaBrDiscriPosition();
  double nai_discri_pos = GetNaIDiscriPosition();
  Double_t REFLaBr = PSD-labr_discri_pos;
  Double_t REFNaI  = PSD-nai_discri_pos;
  Bool_t isNaI = IsPureNaI(PSD,Qs,QL);
  Bool_t isLaBr = IsPureLaBr3(PSD,Qs,QL);
  Bool_t outofLaBr = kFALSE;
  Bool_t outofNaI = kFALSE;
  if(REFLaBr < 0 && !isLaBr) outofLaBr = kTRUE;
  if(REFNaI > 0 && !isNaI) outofNaI = kTRUE;
  if(outofLaBr || outofNaI)  return true;
  else return false;
}

