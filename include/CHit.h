//-------------------------------------------------------------------------//
//                                                                         //
//                               CExperiment.h                             //
//                               Version 1.0                               //
//                        Matthieu Lebois November 2021                    //
//                                                                         //
//  This file contain the definition of the function to define meta        //
//  Parameter to explicitely create everything needed to analyse an        //
//  experiment                                                             //
//                                                                         //
//-------------------------------------------------------------------------//

#ifndef CHIT_H
#define CHIT_H 1

#include "TString.h"
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include "CDetectors.h"
#include "CDetecteurType.h"
#include "colormod.h"

// Definition of the standard types
// Definition of the standard types
using label_Rawtype   = UShort_t;
using label_type      = UShort_t;
using tm_type         = Double_t;
using tm_Rawtype      = ULong64_t;
using nrj_Rawtype     = Int_t;
using nrj_type        = Double_t;
using pu_type         = Bool_t;


// Definition of the I/O function
using namespace std;


//template <typename...> class CHit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//template <class T, class E>
class CHit
{

public:
  CHit(ULong64_t hitnumber);
  virtual ~CHit();

public:
  ULong64_t NHit;

  Int_t Label;
  //template<typename T>
  Double_t Time;
  //template<typename E>
  Double_t NRJ1;
  Double_t NRJ2;
  Double_t NRJ3;
  Double_t NRJ4;
  Bool_t PileUp;

public:
  void SetHitLabel(label_type label){Label = static_cast <Int_t>(label);};
  template <typename T> void SetHitTime(T tm){Time = static_cast <Double_t>(tm);};
  template <typename E> void SetHitE1(E nrj1){NRJ1 = static_cast <Double_t>(nrj1);};
  template <typename E> void SetHitE2(E nrj2){NRJ2 = static_cast <Double_t>(nrj2);};
  template <typename E> void SetHitE3(E nrj3){NRJ3 = static_cast <Double_t>(nrj3);};
  template <typename E> void SetHitE4(E nrj4){NRJ4 = static_cast <Double_t>(nrj4);};
  template <typename E> void SetHitPileUp(pu_type pu){PileUp = static_cast <Bool_t>(pu);};
  template <typename T, typename E> void SetHit(label_type label, T tm, E nrj1, pu_type pu)
  {
    Label = static_cast <Int_t>(label);
    Time = static_cast <Double_t>(tm);
    NRJ1 = static_cast <Double_t>(nrj1);
    NRJ2 = NAN;
    NRJ3 = NAN;
    NRJ4 = NAN;
    PileUp = static_cast <Bool_t>(pu);
  };
  template <typename T, typename E> void SetHit(label_type label, T tm, E nrj1, E nrj2, pu_type pu)
  {
    Label = static_cast <Int_t>(label);
    Time = static_cast <Double_t>(tm);
    NRJ1 = static_cast <Double_t>(nrj1);
    NRJ2 = static_cast <Double_t>(nrj2);
    NRJ3 = NAN;
    NRJ4 = NAN;
    PileUp = static_cast <Bool_t>(pu);
  };
  template <typename T, typename E> void SetHit(label_type label, T tm, E nrj1, E nrj2, E nrj3, pu_type pu)
  {
    Label  = static_cast <Int_t>(label);
    Time   = static_cast <Double_t>(tm);
    NRJ1   = static_cast <Double_t>(nrj1);
    NRJ2   = static_cast <Double_t>(nrj2);
    NRJ3   = static_cast <Double_t>(nrj3);
    NRJ4   = NAN;
    PileUp = static_cast <Bool_t>(pu);
  };
  template <typename T, typename E> void SetHit(label_type label, T tm, E nrj1, E nrj2, E nrj3, E nrj4, pu_type pu)
  {
    Label = static_cast <Int_t>(label);
    Time = static_cast <Double_t>(tm);
    NRJ1 = static_cast <Double_t>(nrj1);
    NRJ2 = static_cast <Double_t>(nrj2);
    NRJ3 = static_cast <Double_t>(nrj3);
    NRJ4 = static_cast <Double_t>(nrj4);
    PileUp = static_cast <Bool_t>(pu);
  };

public:
  Int_t GetHitLabel() const {return Label;};
  Double_t GetHitTime() const {return Time;};
  Double_t GetHitE1() const {return NRJ1;};
  Double_t GetHitE2() const {return NRJ2;};
  Double_t GetHitE3() const {return NRJ3;};
  Double_t GetHitE4() const {return NRJ4;};
  Bool_t GetHitPileUp() const {return PileUp;};

public:
  void PrintHit();
  void Clear();
  Double_t PerformPSD();
  Double_t PerformSimplePSD();
  Double_t PerformPARISPSD();

};

#endif
