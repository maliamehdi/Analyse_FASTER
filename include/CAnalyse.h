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
#include <unordered_map>
#include <array>
#include <numeric>
#include <thread>
#include <algorithm>
#include <type_traits>
#include <exception>
#include <sstream>
#include <functional>
#include <chrono>
#include <tuple>
#include <random>
#include <filesystem>  // C++17


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
#include "TGraph.h"
#include "TSpline.h"
#include "TMultiGraph.h"
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
#include "TTreeIndex.h"
#include <regex>
#include <array>
// I call the other classes from the analyse program
class CExperiment;
class CDetectors;

// Definition of the I/O function
using namespace std;
using namespace indicators;
namespace fs = std::filesystem;
using AlignMap = std::unordered_map<int, std::array<double,3>>;
AlignMap loadAlignFile(const std::string& filename);
using CalibMap = std::unordered_map<std::string, std::pair<double,double>>;
CalibMap loadCalibFile(const std::string& filename);

// General Functions
//int CheckCoincidenceWindow(std::vector<tm_Uncaltype> tab_memory_tm, Double_t deltaTfin);

// Functions about Energy
int CalculateCalibrationCoefficients(const CExperiment &experiment, TString sourcename, Bool_t Calib_BGO, Bool_t Calib_Ge, Bool_t Calib_LaBr);
std::vector<TH1F*> DrawAllEnergyUncalibratedSpectra(const CExperiment &experiment);
int DrawAllEnergyCalibratedSpectra(const CExperiment &experiment, Double_t Emin, Double_t Emax, Int_t binning);
int DrawOneDetectorTypeEnergyCalibratedSpectra(const CExperiment &experiment, const TString DetectorType, Double_t Emin, Double_t Emax, Int_t binning);
int DrawAllCalibrationSpectra(const CExperiment &experiment);
int DrawAllCorrectedCalibrationSpectra(const CExperiment &experiment);
std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>> SetCalibrationSource(TString sourcename, Bool_t Calib_BGO, Bool_t Calib_Ge, Bool_t Calib_LaBr);
std::tuple<int, TGraphErrors*, TF1*, TF1*, TGraph*, Double_t, Double_t , Double_t ,Double_t>
OneSpectrumPeakSearchandAnalysis(const CExperiment &experiment,
                                 Bool_t Calib_BGO, Bool_t Calib_Ge, Bool_t Calib_LaBr,
                                 TH1F *spectrum, int sindex,
                                 std::vector<Double_t> Ref_nrj, std::vector<Double_t> Ref_nrj_err, std::vector<Double_t> Ref_nrj_intensity);
int DrawAllParisUncalibratedSpectra(const CExperiment &experiment);
int DrawAllParisUncalibratedSpectra_with_rotation(const CExperiment &experiment);
int DrawAllParisCalibratedSpectra(const CExperiment &experiment);


std::tuple<std::vector<Double_t>,std::vector<Double_t>>  PSDSpectrumAnalyzer(TH1F *spectrum);
Double_t PSDMatrixAnalyzer(TH2F *matrix);

// Functions about time
std::vector<TH1F*>  DrawTimeShifts_NOTECal(const CExperiment &experiment, Double_t deltaTinit, Double_t deltaTfin);
std::vector<TH1F*>  DrawTimeShifts(const CExperiment &experiment, Double_t deltaTinit, Double_t deltaTfin,  TString filename);
std::vector<TH1F*> DrawTimeShifts_fissionevents(const CExperiment &experiment, Double_t deltaTinit, Double_t deltaTfin);
std::vector<TH1F*> DrawTimeShifts_fissionevents_Calibrated(const CExperiment &experiment, Double_t deltaTinit, Double_t deltaTfin);
std::vector<TH1F*>  FissionEventReconstruction(const CExperiment &experiment, Double_t deltaTinit, Double_t deltaTfin, Double_t neutronwindow);
int BISFissionEventReconstruction(const CExperiment &experiment, Double_t deltaTinit, Double_t deltaTfin, Double_t neutronwindow);
int CalculateTimealignementShifts(const CExperiment &experiment, Bool_t isCalibrated);
std::vector<Double_t> DeltaTmeasurer(TH1F *timespectrum, bool isqdc);
int CheckCoincidenceWindow(std::vector<tm_Rawtype> tab_memory_tm, Double_t deltaTfin);
std::vector<TH1F*>  CheckTimeShifts(const CExperiment &experiment, Double_t deltaTinit, Double_t deltaTfin);

// Functions to manipulate data (calibration, time alignement, sorting...)
int EnergyCalibrator(const CExperiment &experiment);
int NRJCalibrator(const CExperiment &experiment,TChain *chained_oak, DynamicProgress<ProgressBar> &BAR, int localchainnbr, mutex *stop);
const std::string ApplyMyEnergyCalibration(const CExperiment &experiment);
const std::string ApplyMyCorrection(const CExperiment &experiment);
const std::string TimeAlignator(const CExperiment &experiment, Bool_t isCalibrated);
double alignCalib(const AlignMap& data, int detectorId, double x);
CHitCollection* BuildCenteredWindow(TChain* chain, Double_t centerTime_ns, Double_t window_ns, ULong64_t startIndex);

UInt_t CompressFASTERValue(const nrj_Rawtype &enrj, const int &originalBits, const int &targetBits);

// void functions
// Prototype de la fonction :
void SortTreeByTime(const std::string &treename, const std::string &filename);

void checktimeorder(const std::string &filename, const std::string &treename);

void generate_dat_files_ECal_TShift(const std::string& ecal_file_path);
void generate_dat_files_Ecal(const std::string& ecal_file_path);
void generate_dat_files_CORR(const std::string& CORR_file_path);
void LoadCorrectionParameters(const std::string& filename);
// Structure pour stocker les coefficients
struct CorrectionParams {
    double a_opt = 1.0;
    double b_opt = 0.0;
};

// Map globale
extern std::map<std::pair<int, std::string>, CorrectionParams> correctionMap;
/*
int CalculateResidues(const CExperiment &experiment, TString sourcename, Bool_t Calib_BGO, Bool_t Calib_Ge, Bool_t Calib_LaBr);
std::vector<TH1F*> DrawAllEnergyUncalibratedSpectra(const CExperiment &experiment);
std::vector<TH1F*> DrawAllEnergyCalibratedSpectra(const CExperiment &experiment, Double_t Emin, Double_t Emax, Int_t binning);
TH1F* DrawOneEnergyUncalibratedSpectra(const CExperiment &experiment, Int_t sindex);// Not DAQ label
Double_t Faster2bitsNRJConverter(const CExperiment &experiment, Int_t detectorindex);


// Functions about time
std::vector<TH1F*>  DrawTimeShifts(const CExperiment &experiment, Double_t deltaTinit, Double_t deltaTfin,  TString filename);
std::vector<TH1F*>  DrawTimeShifts_ECal(const CExperiment &experiment, Double_t deltaTinit, Double_t deltaTfin);
int CheckCoincidenceWindow(std::vector<tm_Uncaltype> tab_memory_tm, Double_t deltaTfin);
int CalculateTimealignementShifts(const CExperiment &experiment, Bool_t isCalibrated);
int CheckTimealignementShifts(const CExperiment &experiment, Bool_t isCalibrated);
std::vector<Double_t> DeltaTmeasurer(TH1F *timespectrum, bool isqdc);
int TimeAlignator(const CExperiment &experiment);
int TMCalibrator(const CExperiment &experiment, TChain *chained_oak);
int TimeSorter(const CExperiment &experiment);
int TMSorter(const CExperiment &experiment, TChain *chained_oak);
*/

#endif //TANALYSEF_H
