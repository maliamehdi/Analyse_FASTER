//-------------------------------------------------------------------------//
//                                                                         //
//                                CAnalyse.h                               //
//                               Version 1.0                               //
//                        Matthieu Lebois November 2021                    //
//                                                                         //
//  This file contain the definition of the function for a root file       //
//  and the declaration to extract some physical results from the data     //
//  contained in the root file                                             //
//                                                                         //
//-------------------------------------------------------------------------//

#ifndef CANALYSE_H
#define CANALYSE_H 1

// Call for the principal parts of the functions
#include "Analyse_Main.h"
#include "CExperiment.h"
#include "CDetectors.h"
#include "CDetecteurType.h"
#include "colormod.h"
#include "CMTTHxx.h"
#include "CHit.h"
#include "CHitCollection.h"
#include "CTools.h"
#include "ProgressBar.h"

// Call for standard C/C++ Libraries
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <map>
#include <numeric>
#include <thread>
#include <algorithm>
#include <type_traits>
#include <exception>
#include <sstream>
#include <functional>
#include <chrono>
#include <tuple>

// Call for ROOT Libraries
#include <TROOT.h>
#include <ROOT/TSeq.hxx>
#include <ROOT/TProcessExecutor.hxx>
#include <ROOT/TThreadedObject.hxx>
#include <ROOT/TTreeProcessorMT.hxx>
#include "TThread.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TSpectrum.h"
#include "TApplication.h"
#include <TObject.h>
#include <TFriendElement.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TEntryList.h>
#include <TEnv.h>
#include "TFile.h"
#include "TCutG.h"
#include "TMath.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH3.h"
#include "TH3F.h"
#include "THStack.h"
#include "TObject.h"
#include "TTree.h"
#include "TBasket.h"
#include "TString.h"
#include "TBranch.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TSpectrum2.h"
#include "TF2.h"

// I call the other classes from the analyse program
class CExperiment;
class CDetectors;

// Definition of the I/O function
using namespace std;
using namespace indicators;

// General Functions
//int CheckCoincidenceWindow(std::vector<tm_Uncaltype> tab_memory_tm, Double_t deltaTfin);

// Functions about Energy
int CalculateCalibrationCoefficients(CExperiment *experiment, TString sourcename, Bool_t Calib_BGO, Bool_t Calib_Ge, Bool_t Calib_LaBr);
std::vector<TH1F*> DrawAllEnergyUncalibratedSpectra(CExperiment *experiment);
int DrawAllCalibrationSpectra(CExperiment *experiment);
std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>> SetCalibrationSource(TString sourcename, Bool_t Calib_BGO, Bool_t Calib_Ge, Bool_t Calib_LaBr);
std::tuple<int, TGraphErrors*, TF1*, TF1*, TGraph*, Double_t, Double_t , Double_t ,Double_t>
OneSpectrumPeakSearchandAnalysis(CExperiment *experiment,
                                 Bool_t Calib_BGO, Bool_t Calib_Ge, Bool_t Calib_LaBr,
                                 TH1F *spectrum, int sindex,
                                 std::vector<Double_t> Ref_nrj, std::vector<Double_t> Ref_nrj_err, std::vector<Double_t> Ref_nrj_intensity);
int DrawAllParisUncalibratedSpectra(CExperiment *experiment);
std::tuple<std::vector<Double_t>,std::vector<Double_t>>  PSDSpectrumAnalyzer(TH1F *spectrum);
Double_t PSDMatrixAnalyzer(TH2F *matrix);

// Functions about time
std::vector<TH1F*>  DrawTimeShifts_ECal(CExperiment *experiment, Double_t deltaTinit, Double_t deltaTfin);
std::vector<TH1F*>  DrawTimeShifts(CExperiment *experiment, Double_t deltaTinit, Double_t deltaTfin,  TString filename);
int CalculateTimealignementShifts(CExperiment *experiment, Bool_t isCalibrated);
std::vector<Double_t> DeltaTmeasurer(TH1F *timespectrum, bool isqdc);
int CheckCoincidenceWindow(std::vector<tm_Rawtype> tab_memory_tm, Double_t deltaTfin);

// Functions to manipulate data (calibration, time alignement, sorting...)
int EnergyCalibrator(CExperiment *experiment);
int NRJCalibrator(CExperiment *experiment,TChain *chained_oak, DynamicProgress<ProgressBar> &BAR, int localchainnbr, mutex *stop);

/*
int CalculateResidues(CExperiment *experiment, TString sourcename, Bool_t Calib_BGO, Bool_t Calib_Ge, Bool_t Calib_LaBr);
std::vector<TH1F*> DrawAllEnergyUncalibratedSpectra(CExperiment *experiment);
std::vector<TH1F*> DrawAllEnergyCalibratedSpectra(CExperiment *experiment, Double_t Emin, Double_t Emax, Int_t binning);
TH1F* DrawOneEnergyUncalibratedSpectra(CExperiment *experiment, Int_t sindex);// Not DAQ label
Double_t Faster2bitsNRJConverter(CExperiment *experiment, Int_t detectorindex);


// Functions about time
std::vector<TH1F*>  DrawTimeShifts(CExperiment *experiment, Double_t deltaTinit, Double_t deltaTfin,  TString filename);
std::vector<TH1F*>  DrawTimeShifts_ECal(CExperiment *experiment, Double_t deltaTinit, Double_t deltaTfin);
int CheckCoincidenceWindow(std::vector<tm_Uncaltype> tab_memory_tm, Double_t deltaTfin);
int CalculateTimealignementShifts(CExperiment *experiment, Bool_t isCalibrated);
int CheckTimealignementShifts(CExperiment *experiment, Bool_t isCalibrated);
std::vector<Double_t> DeltaTmeasurer(TH1F *timespectrum, bool isqdc);
int TimeAlignator(CExperiment *experiment);
int TMCalibrator(CExperiment *experiment, TChain *chained_oak);
int TimeSorter(CExperiment *experiment);
int TMSorter(CExperiment *experiment, TChain *chained_oak);
*/

#endif //TANALYSEF_H
