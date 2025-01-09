// Call for the principal parts of the functions
#include "CHit.h"
#include <stdlib.h>


// Definition of the I/O function
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
CHit::CHit(ULong64_t hitnumber)
{
  // Defining the name of the experiment
  NHit = hitnumber;
  Label = 0; Time = 0; NRJ1=0; NRJ2=0; NRJ3=0; NRJ4=0; PileUp = kFALSE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
CHit::~CHit(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void CHit::Clear()
{
  Label = 0; Time = 0; NRJ1=0; NRJ2=0; NRJ3=0; NRJ4=0; PileUp = kFALSE;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void CHit::PrintHit()
{
  cout << RESETTEXT << FOREMAG << "############################################" << endl;
  cout << "# Printing Hit #" << NHit <<  endl;
  cout << "############################################" << RESETTEXT << endl;
  cout << "# Detector Label    " << Label << endl;
  cout << "# Hit time          " << (ULong64_t)Time << endl;
  cout << "# Hit First Energy  " << NRJ1 << endl;
  cout << "# Hit Second Energy " << NRJ2 << endl;
  cout << "# Hit Third Energy  " << NRJ3 << endl;
  cout << "# Hit Fourth Energy " << NRJ4 << endl;
  cout << "# Hit Pile-Up       " << PileUp << endl;
  cout << "############################################" << RESETTEXT << endl <<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Double_t CHit::PerformPSD()
{
  // First I identify long from short gate:
  // if(NRJ1 < NRJ2) // NRJ2 is the long gate
  // {
  //   return (NRJ2-NRJ1)/NRJ2;
  // }
  // else
  // {
  //   return (NRJ1-NRJ2)/NRJ1;
  // }
  return (NRJ2-NRJ1)/NRJ2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Double_t CHit::PerformSimplePSD()
{
  return (NRJ1/NRJ2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Double_t CHit::PerformPARISPSD()
{
  return TMath::ATan(static_cast<double>(NRJ1)/static_cast<double>(NRJ2));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......