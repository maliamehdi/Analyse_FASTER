//-------------------------------------------------------------------------//
//                                                                         //
//                            CHitCollection.h                             //
//                               Version 1.0                               //
//                        Matthieu Lebois November 2021                    //
//                                                                         //
//  This file contain the definition of the function to define meta        //
//  Parameter to explicitely create everything needed to analyse an        //
//  experiment                                                             //
//                                                                         //
//-------------------------------------------------------------------------//

#ifndef CHitCollection_h
#define CHitCollection_h 1

#include "TString.h"
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include "CDetectors.h"
#include "CDetecteurType.h"
#include "colormod.h"
#include "CHit.h"
#include "CExperiment.h"

// Definition of the standard types
using label_type = UShort_t;
using tm_Caltype = Double_t;
using tm_Uncaltype = ULong64_t;
using nrj_Ucaltype = Int_t;
using nrj_Caltype = Double_t;
using pu_type = Bool_t;

// Definition of the I/O function
using namespace std;

// I call the other classes from the analyse program
class CExperiment;
class CDetectors;

class CHitCollection
{
public:
  CHitCollection();
  virtual ~CHitCollection();

  void PrintHitCollection();
  void Clear();

public:
  UShort_t collectionsize;
  Double_t collectionTimesize; // Time in ns
  Double_t averageenergy;
  Double_t totalenergy;

  typename vector<CHit*>::iterator it_hitcollection;
  vector<CHit*> hitcollection;

// the Set functions
public:
  void SetCollectionSize(int size){collectionsize = (UShort_t)size;};
  void SetCollectionTimeSize(Double_t deltaT){collectionTimesize = deltaT;};

// The Get functions
public:
  int GetCollectionSize(){return (int)collectionsize;};
  Double_t GetCollectionTimeSize(){return collectionTimesize;};
  Double_t GetAverageEnergy();
  Double_t GetTotalEnergy();
  Int_t GetModularMultiplicity(CExperiment *experiment);
  void DoComptonSuppression(CExperiment *experiment);
  void DoAddBack(CExperiment *experiment);
  void DoABandCS(CExperiment *experiment);

public:
  void AddHit(CHit *hit){hitcollection.push_back(hit);collectionsize++;};
  CHit GetHit(int i){return *(hitcollection.at(i));}
  Bool_t IsHitInside(CHit *hit);
  int IsReferenceDetectorIn(int referencelabel);
  int  CountLabel(int label) const;      // how many times does 'label' appear?
  Bool_t HasLabel(int label) const;        // at least once?
  int  FindLabel(int label) const;       // index of first occurrence or -1 if none
  
};


#endif
