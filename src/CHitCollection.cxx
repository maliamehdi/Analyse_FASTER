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
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CHitCollection::Clear()
{
  // I empty the Hit vector
  hitcollection.clear();
  hitcollection.shrink_to_fit();

  // I reinitialize the class variables
  collectionsize = 0;
  averageenergy = NAN;
  totalenergy = NAN;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CHitCollection::PrintHitCollection()
{

  cout << RESETTEXT << FOREMAG << BACKWHT << "############################################" << BACKBLK << endl;
  cout << BACKWHT << "# Printing Hit Collection " << BACKBLK <<  endl;
  cout << BACKWHT << "############################################" << RESETTEXT << endl;

  for(auto &&hit : hitcollection)
  {
    hit.PrintHit();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t CHitCollection::GetAverageEnergy()
{
  if(averageenergy !=averageenergy)
  {
    averageenergy = 0;
    for(int i=0; i < collectionsize; i++)
    {
      averageenergy += hitcollection.at(i).GetHitE1();
    }
    averageenergy = averageenergy/collectionsize;
  }
  return averageenergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t CHitCollection::GetTotalEnergy()
{
  if(totalenergy !=totalenergy)
  {
    totalenergy = 0;
    for(int i=0; i < collectionsize; i++)
    {
      totalenergy += hitcollection.at(i).GetHitE1();
    }
  }
  return totalenergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Bool_t CHitCollection::IsHitInside(CHit *hit)
{
  if(collectionsize == 0) return true;

  else
  {
    Double_t deltaT = (hit->GetHitTime()-hitcollection.at(0).GetHitTime())/1000.; // Time in ns
    if (deltaT <= collectionTimesize) return true;
    else return false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int CHitCollection::IsReferenceDetectorIn(int referencelabel)
{
  if(collectionsize == 0) return -666;

  else
  {
    int finalpos = -666;
    int pos(0);
    for(auto &&hit : hitcollection)
    {
      if(hit.GetHitLabel() == referencelabel) finalpos = pos;
      pos++;
    }
    return finalpos;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Int_t CHitCollection::GetModularMultiplicity(CExperiment *experiment)
{

  // I read the event
  int modularmultiplicity(0);
  int *tab_module = new int[experiment->GetNbrofModules()+1];
  int alveolenbr(0);int ringnbr(0);

  //channeltag = ringnbr*10000+alveolenbr*100+type*10+colornbr;

  for(int i=0; i<experiment->GetNbrofModules()+1; i++) tab_module[i]=0;

  for(int ll = 0; ll < collectionsize;ll++)
  {
    int indice(0);
    //int detindex = experiment->GetLabel2Detnbrs()[hitcollection.at(ll).GetHitLabel()];
    int detindex = experiment->GetLabel2Detnbrs(hitcollection.at(ll).GetHitLabel());
    alveolenbr = experiment->GetDetector(detindex)->GetAlveole();
    ringnbr = experiment->GetDetector(detindex)->GetRing();

    // Determining which alveole it belongs to
    indice = ringnbr*10-10+alveolenbr;
    switch(ringnbr)
    {
      case 1:
        indice =  0+alveolenbr;
        break;
      case 2:
        indice =  10+alveolenbr;
        break;
      case 3:
        indice =  20+alveolenbr;
        break;
      case 4:
        indice =  30+alveolenbr;
        break;
      default:
        break;
    }
    tab_module[indice]++;
  }
  for(int i=0; i<experiment->GetNbrofModules()+1; i++)
  {
    if(tab_module[i]!=0)modularmultiplicity++;
  }
  return modularmultiplicity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CHitCollection::DoComptonSuppression(CExperiment *experiment)
{

  // I create an array of vector that will keep in memory the position in the
  // collection that belongs to the same position
  std::vector<int> *tab_module = new vector<int>[experiment->GetNbrofModules()+1];
  //int *tab_module = new int[experiment->GetNbrofModules()+1];
  int alveolenbr(0);int ringnbr(0);
  std::vector<Int_t> tab_CS_label;
  tab_CS_label.clear();
  bool CS_cond=false;

  for(int ll = 0; ll < collectionsize;ll++)
  {
    int indice(0);
    //int detindex = experiment->GetLabel2Detnbrs()[hitcollection.at(ll).GetHitLabel()];
    int detindex = experiment->GetLabel2Detnbrs(hitcollection.at(ll).GetHitLabel());
    alveolenbr = experiment->GetDetector(detindex)->GetAlveole();
    ringnbr = experiment->GetDetector(detindex)->GetRing();

    // Determining which alveole it belongs to
    switch(ringnbr)
    {
      case 1:
        indice =  0+alveolenbr;
        break;
      case 2:
        indice =  10+alveolenbr;
        break;
      case 3:
        indice =  20+alveolenbr;
        break;
      case 4:
        indice =  30+alveolenbr;
        break;
      default:
        break;
    }
    tab_module[indice].push_back(ll);
  }

  // Normally I just have a set of position in the hit collection
  // that might belong to the same position
  for(int indice=0; indice < experiment->GetNbrofModules()+1; indice++)
  {
    // I check I found two detector (or more) fired at the same position
    if(tab_module[indice].size()>1)
    {
      Bool_t isBGO = kFALSE; int BGOpos(666);
      Bool_t isGe = kFALSE; int Gepos(666);
      for(auto jj = 0; jj<(int)tab_module[indice].size(); jj++)
      {
        //int detindex = experiment->GetLabel2Detnbrs()[hitcollection.at(tab_module[indice].at(jj)).GetHitLabel()];
        int detindex = experiment->GetLabel2Detnbrs(hitcollection.at(tab_module[indice].at(jj)).GetHitLabel());
        Double_t nrj = hitcollection.at(tab_module[indice].at(jj)).GetHitE1();
        if(experiment->GetDetector(detindex)->GetDetectorType()=="Ge") {isGe = kTRUE;Gepos=jj;}
        if(experiment->GetDetector(detindex)->GetDetectorType()=="BGO" && nrj > 4000.) {isBGO = kTRUE;BGOpos=jj;}
      }
      if(isBGO && isGe)
      {
        // I have to check the timing for Compton suppression
        Double_t timeGe = hitcollection.at(tab_module[indice].at(Gepos)).GetHitTime();
        Double_t timeBGO = hitcollection.at(tab_module[indice].at(BGOpos)).GetHitTime();
        Double_t deltaT = (timeGe-timeBGO);deltaT = deltaT/1000.;
        if(deltaT > -50. && deltaT < 50.)
        {
          CS_cond=true;
          for(auto jj = 0; jj<(int)tab_module[indice].size(); jj++) tab_CS_label.push_back(tab_module[indice].at(jj));
        }
      }
    }
  }

  if(CS_cond == true){
    // I need to sort the CS_label to erase the last term first
    std::sort(tab_CS_label.begin(),tab_CS_label.end());

    //for (int ll = tab_CS_label.size()-1; ll >=0;ll--) std::cout << tab_CS_label.at(ll) << std::endl;

    // Now I loop for the last to first entry to remove for CS the Ge and  its BGO
    int memorypos(666);
    for(int ll = tab_CS_label.size()-1; ll >= 0;ll--)
    {
      int pos = tab_CS_label.at(ll);
      //std::cout << "pos = " << pos << std::endl;
      if(pos != memorypos){
        hitcollection.erase(hitcollection.begin()+pos);  // Memorisation of the detector number to check for coincidences
      }
      memorypos = pos;
    }
  }


  collectionsize=hitcollection.size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CHitCollection::DoAddBack(CExperiment *experiment)
{
  // I create an array of vector that will keep in memory the position in the
  // collection that belongs to the same position
  std::vector<int> *tab_module = new vector<int>[experiment->GetNbrofModules()+1];
  //int *tab_module = new int[experiment->GetNbrofModules()+1];
  int alveolenbr(0);int ringnbr(0);
  std::vector<Int_t> tab_CS_label;
  tab_CS_label.clear();
  bool CS_cond=false;

  for(int ll = 0; ll < collectionsize;ll++)
  {
    int indice(0);
    int detindex = experiment->GetLabel2Detnbrs(hitcollection.at(ll).GetHitLabel());
    alveolenbr = experiment->GetDetector(detindex)->GetAlveole();
    ringnbr = experiment->GetDetector(detindex)->GetRing();

    // Determining which alveole it belongs to
    //indice = ringnbr*10-10+alveolenbr;
    switch(ringnbr)
    {
      case 1:
        indice =  0+alveolenbr;
        break;
      case 2:
        indice =  10+alveolenbr;
        break;
      case 3:
        indice =  20+alveolenbr;
        break;
      case 4:
        indice =  30+alveolenbr;
        break;
      default:
        break;
    }
    tab_module[indice].push_back(ll);
  }

  // Normally I just have a set of position in the hit collection
  // that might belong to the same position
  for(int indice=0; indice < experiment->GetNbrofModules()+1; indice++)
  {
    // I check I found two detector (or more) fired at the same position
    if(tab_module[indice].size()>1)
    {
      int nGe=0; /// I'm going to count the number of crystal fired
      std::vector<int> Gepos;
      Bool_t isGe = kFALSE;
      for(auto jj = 0; jj<(int)tab_module[indice].size(); jj++)
      {
        int detindex = experiment->GetLabel2Detnbrs(hitcollection.at(tab_module[indice].at(jj)).GetHitLabel());
        if(experiment->GetDetector(detindex)->GetDetectorType()=="Ge") {nGe++;Gepos.push_back(jj);}
      }

      if(nGe > 1) // No matter what I will the first two neighboring crystals
      {
        for(int ll=0; ll < (int)Gepos.size(); ll++)
        {
          int detindex = experiment->GetLabel2Detnbrs(hitcollection.at(tab_module[indice].at(Gepos.at(ll))).GetHitLabel());
          int HPGEID1 = experiment->GetDetector(detindex)->GetPID();

          // Now I scan the others if this one has not been "suppressed" earlier
          if(hitcollection.at(tab_module[indice].at(Gepos.at(ll))).GetHitE1()>0.)
          {
            for(int kk = 0; kk < (int)Gepos.size(); kk++ )
            {
              if(kk!=ll)
              {
                detindex = experiment->GetLabel2Detnbrs(hitcollection.at(tab_module[indice].at(Gepos.at(kk))).GetHitLabel());
                int HPGEID2 = experiment->GetDetector(detindex)->GetPID();

                int PIDdiff = TMath::Abs(HPGEID1-HPGEID2);
                if(PIDdiff < 3 && PIDdiff%2 == 0)
                {
                  // Now I check the DeltaT to make sure
                  Double_t deltaT = (Double_t)(hitcollection.at(tab_module[indice].at(Gepos.at(kk))).GetHitTime()-hitcollection.at(tab_module[indice].at(Gepos.at(ll))).GetHitTime()); // to get in ns
                  deltaT = deltaT/1000.;
                  if(deltaT > -50. && deltaT < 50.)
                  {
                    // I've got to add it to the previous crystal
                    hitcollection.at(tab_module[indice].at(Gepos.at(ll))).SetHitE1(hitcollection.at(tab_module[indice].at(Gepos.at(ll))).GetHitE1()+hitcollection.at(tab_module[indice].at(Gepos.at(kk))).GetHitE1());

                    // I set the other crystal energy to zero not to sum it later
                    hitcollection.at(tab_module[indice].at(Gepos.at(kk))).SetHitE1(-999.);
                    tab_CS_label.push_back(tab_module[indice].at(Gepos.at(kk)));
                    CS_cond=true;
                  }
                }
              }
            }
          }

        }
      }
    }
  }

  if(CS_cond == true){
    // I need to sort the CS_label to erase the last term first
    std::sort(tab_CS_label.begin(),tab_CS_label.end());

    //for (int ll = tab_CS_label.size()-1; ll >=0;ll--) std::cout << tab_CS_label.at(ll) << std::endl;

    // Now I loop for the last to first entry to remove for CS the Ge and  its BGO
    int memorypos(666);
    for(int ll = tab_CS_label.size()-1; ll >= 0;ll--)
    {
      int pos = tab_CS_label.at(ll);
      //std::cout << "pos = " << pos << std::endl;
      if(pos != memorypos){
        if(hitcollection.at(pos).GetHitE1() < 0) hitcollection.erase(hitcollection.begin()+pos);  // Memorisation of the detector number to check for coincidences
      }
      memorypos = pos;
    }
  }


  collectionsize=hitcollection.size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CHitCollection::DoABandCS(CExperiment *experiment)
{
  // I create an array of vector that will keep in memory the position in the
  // collection that belongs to the same position
  std::vector<int> *tab_module = new vector<int>[experiment->GetNbrofModules()+1];
  //int *tab_module = new int[experiment->GetNbrofModules()+1];
  int alveolenbr(0);int ringnbr(0);
  std::vector<Int_t> tab_CS_label;
  tab_CS_label.clear();
  bool CS_cond=false;

  for(int ll = 0; ll < collectionsize;ll++)
  {
    int indice(0);
    int detindex = experiment->GetLabel2Detnbrs(hitcollection.at(ll).GetHitLabel());
    alveolenbr = experiment->GetDetector(detindex)->GetAlveole();
    ringnbr = experiment->GetDetector(detindex)->GetRing();

    // Determining which alveole it belongs to
    //indice = ringnbr*10-10+alveolenbr;
    switch(ringnbr)
    {
      case 1:
        indice =  0+alveolenbr;
        break;
      case 2:
        indice =  10+alveolenbr;
        break;
      case 3:
        indice =  20+alveolenbr;
        break;
      case 4:
        indice =  30+alveolenbr;
        break;
      default:
        break;
    }
    tab_module[indice].push_back(ll);
  }

  // Normally I just have a set of position in the hit collection
  // that might belong to the same position
  for(int indice=0; indice < experiment->GetNbrofModules()+1; indice++)
  {
    // I check I found two detector (or more) fired at the same position
    if(tab_module[indice].size()>1)
    {
      int nGe=0; /// I'm going to count the number of crystal fired
      Bool_t isBGO = kFALSE; int BGOpos(666);
      std::vector<int> Gepos;
      Bool_t isGe = kFALSE;
      for(auto jj = 0; jj<(int)tab_module[indice].size(); jj++)
      {
        int detindex = experiment->GetLabel2Detnbrs(hitcollection.at(tab_module[indice].at(jj)).GetHitLabel());
        Double_t nrj = hitcollection.at(tab_module[indice].at(jj)).GetHitE1();
        if(experiment->GetDetector(detindex)->GetDetectorType()=="Ge") {isGe = kTRUE;nGe++;Gepos.push_back(jj);}
        if(experiment->GetDetector(detindex)->GetDetectorType()=="BGO" && nrj > 4000.) {isBGO = kTRUE;BGOpos=jj;}
      }
      if(isBGO && isGe)
      {
        // I have to check the timing for Compton suppression
        for(int kk = 0; kk < Gepos.size(); kk++)
        {
          Double_t timeGe = hitcollection.at(tab_module[indice].at(Gepos.at(kk))).GetHitTime();
          Double_t timeBGO = hitcollection.at(tab_module[indice].at(BGOpos)).GetHitTime();
          Double_t deltaT = (timeGe-timeBGO);deltaT = deltaT/1000.;
          if(deltaT > -50. && deltaT < 50.)
          {
            CS_cond=true;
            for(auto jj = 0; jj<(int)tab_module[indice].size(); jj++) tab_CS_label.push_back(tab_module[indice].at(jj));
          }
        }

      }
      if(nGe > 1) // No matter what I will the first two neighboring crystals
      {
        for(int ll=0; ll < (int)Gepos.size(); ll++)
        {
          int detindex = experiment->GetLabel2Detnbrs(hitcollection.at(tab_module[indice].at(Gepos.at(ll))).GetHitLabel());
          int HPGEID1 = experiment->GetDetector(detindex)->GetPID();

          // Now I scan the others if this one has not been "suppressed" earlier
          if(hitcollection.at(tab_module[indice].at(Gepos.at(ll))).GetHitE1()>0.)
          {
            for(int kk = 0; kk < (int)Gepos.size(); kk++ )
            {
              if(kk!=ll)
              {
                detindex = experiment->GetLabel2Detnbrs(hitcollection.at(tab_module[indice].at(Gepos.at(kk))).GetHitLabel());
                int HPGEID2 = experiment->GetDetector(detindex)->GetPID();

                int PIDdiff = TMath::Abs(HPGEID1-HPGEID2);
                if(PIDdiff < 3 && PIDdiff%2 == 0)
                {
                  // Now I check the DeltaT to make sure
                  Double_t deltaT = (Double_t)(hitcollection.at(tab_module[indice].at(Gepos.at(kk))).GetHitTime()-hitcollection.at(tab_module[indice].at(Gepos.at(ll))).GetHitTime()); // to get in ns
                  deltaT = deltaT/1000.;
                  if(deltaT > -50. && deltaT < 50.)
                  {
                    // I've got to add it to the previous crystal
                    hitcollection.at(tab_module[indice].at(Gepos.at(ll))).SetHitE1(hitcollection.at(tab_module[indice].at(Gepos.at(ll))).GetHitE1()+hitcollection.at(tab_module[indice].at(Gepos.at(kk))).GetHitE1());

                    // I set the other crystal energy to zero not to sum it later
                    hitcollection.at(tab_module[indice].at(Gepos.at(kk))).SetHitE1(-999.);
                    tab_CS_label.push_back(tab_module[indice].at(Gepos.at(kk)));
                    CS_cond=true;
                  }
                }
              }
            }
          }

        }
      }
    }
  }

  if(CS_cond == true){
    // I need to sort the CS_label to erase the last term first
    std::sort(tab_CS_label.begin(),tab_CS_label.end());
    tab_CS_label.erase( unique( tab_CS_label.begin(), tab_CS_label.end() ), tab_CS_label.end() );

    //for (int ll = tab_CS_label.size()-1; ll >=0;ll--) std::cout << tab_CS_label.at(ll) << std::endl;

    // Now I loop for the last to first entry to remove for CS the Ge and  its BGO
    int memorypos(666);
    for(int ll = tab_CS_label.size()-1; ll >= 0;ll--)
    {
      int pos = tab_CS_label.at(ll);
      //std::cout << "pos = " << pos << std::endl;
      if(pos != memorypos){
        hitcollection.erase(hitcollection.begin()+pos);  // Memorisation of the detector number to check for coincidences
      }
      memorypos = pos;
    }
  }


  collectionsize=hitcollection.size();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
