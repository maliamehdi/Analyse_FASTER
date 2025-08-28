// Call for the principal parts of the functions
#include "CHitCollection.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <vector>
#include <math.h>

// Definition of the I/O function
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CHitCollection::CHitCollection()
{
  collectionsize = 0;
  collectionTimesize = 100.; // Time windows size set by default to 100 ns
  averageenergy = NAN;
  totalenergy = NAN;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CHitCollection::~CHitCollection()
{
  Clear(); // 🔥 libère tout à la destruction
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CHitCollection::Clear()
{
  for (CHit* hit : hitcollection) {
    delete hit;
  }
  hitcollection.clear();
  hitcollection.shrink_to_fit();

  collectionsize = 0;
  averageenergy = NAN;
  totalenergy = NAN;
}

void CHitCollection::PrintHitCollection()
{
  cout << RESETTEXT << FOREMAG << BACKWHT << "############################################" << BACKBLK << endl;
  cout << BACKWHT << "# Printing Hit Collection " << BACKBLK <<  endl;
  cout << BACKWHT << "############################################" << RESETTEXT << endl;

  for (CHit* hit : hitcollection)
  {
    hit->PrintHit();
  }
}

Double_t CHitCollection::GetAverageEnergy()
{
  if (averageenergy != averageenergy)
  {
    averageenergy = 0;
    for (int i = 0; i < collectionsize; i++)
    {
      averageenergy += hitcollection.at(i)->GetHitE1();
    }
    averageenergy /= collectionsize;
  }
  return averageenergy;
}

Double_t CHitCollection::GetTotalEnergy()
{
  if (totalenergy != totalenergy)
  {
    totalenergy = 0;
    for (int i = 0; i < collectionsize; i++)
    {
      totalenergy += hitcollection.at(i)->GetHitE1();
    }
  }
  return totalenergy;
}

Bool_t CHitCollection::IsHitInside(CHit *hit)
{
  if (collectionsize == 0) return true;
  else
  {
    Double_t deltaT = (hit->GetHitTime() - hitcollection.at(0)->GetHitTime()) / 1000.;
    return (deltaT <= collectionTimesize);
  }
}

int CHitCollection::IsReferenceDetectorIn(int referencelabel)
{
  if (collectionsize == 0) return -666;

  int finalpos = -666;
  int pos = 0;
  for (CHit* hit : hitcollection)
  {
    if (hit->GetHitLabel() == referencelabel) finalpos = pos;
    pos++;
  }
  return finalpos;
}

int CHitCollection::CountLabel(int label) const
{
  int count = 0;
  for (size_t i = 0; i < hitcollection.size(); i++) {
    if (hitcollection[i]->GetHitLabel() == label) {
      count++;
    }
  }
  return count;
}

bool CHitCollection::HasLabel(int label) const
{
  return (CountLabel(label) > 0);
}

int CHitCollection::FindLabel(int label) const
{
  for (size_t i = 0; i < hitcollection.size(); i++) {
    if (hitcollection[i]->GetHitLabel() == label) {
      return static_cast<int>(i);
    }
  }
  return -1;
}