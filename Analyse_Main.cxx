//-------------------------------------------------------------------------//
//                                                                         //
//                             Analyse_Main.cxx                            //
//                               Version 1.0                               //
//                       Matthieu Lebois November 2021                     //
//                                                                         //
//   This program make some useful graphe from the root file created by    //
//                               Convert_FASTER                            //
//                                                                         //
//-------------------------------------------------------------------------//

// Call of the project general library
#include "Analyse_Main.h"

// Definition of the I/O function
using namespace std;

//-------------------------------------------------------------//
//                                                             //
//                    Principal part of the program            //
//                                                             //
//-------------------------------------------------------------//
// Main part of the program
int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <runlist>" << std::endl;
    return 1;
  }

  
  // Creation of a time counter
  TStopwatch Configurationtimer;
  TStopwatch Analysistimer;

  // Chronometer to benchmark the graphes constitution process
  Configurationtimer.Reset();
  Configurationtimer.Start();

  // Creation of an experiment.
  TString expname = "frozen";
  CExperiment *experiment = new CExperiment(expname);

  // Defining the number of threads available on the machine
  const auto max_workers = std::thread::hardware_concurrency();
  int nworkers = static_cast<int>(2*max_workers/3.);
  experiment->SetThreadsNbr(nworkers);
  ROOT::EnableImplicitMT(experiment->GetThreadsNbr());

  // Name of the folder to find various parametrization Parameter_Files
  TString folder_name = "FROZEN/frozen/";//"FROZEN/Test/";//"nuball2/Calibrations/";//"FROZEN/PoC/";
  const char* runlist = argv[1];
  TString Runlist_filename = runlist;//"rootlist_premanip_60Co.dat";
  // Names of the calibration file list
  //"listerun_calib_NRJ_152Eu_endofexperiment.dat";
  //"listerun_AmBe_endofexperiment.dat";
  //"listerun_CoandCs_endexperiment.dat";
  //"listerun_calib_NRJ_152Eu.dat";
  //"listerun_calib_NRJ_152Eu_Cal.dat";
  //"listerun_PARIS_Calibration.dat";//"listerun_CoandCs_endexperiment.dat";//"listerun_nearline.dat";//"listerun_calib_Time.dat";//"listerun_calib_NRJ.dat";

  // Some Boolean to know what I'm doing with Data
  Bool_t isPlottingSpectra = kTRUE;
  Bool_t isCalculating_E_Calibration = kFALSE;
  Bool_t isApplying_E_Calibration = kFALSE;
  Bool_t isChecking_E_Calib = kFALSE;
  Bool_t isCalculating_T_Calibration = kFALSE;
  Bool_t isApplying_T_Calibration = kFALSE;
  Bool_t isChecking_T_Calib = kFALSE;
  Bool_t isEventBuilding = kFALSE;
  Bool_t isEventAnalysing = kFALSE;

  // Multichain is for multithreading the data analysis
  // kFALSE : one TChain
  // kTRUE : Several TChains sent to several threads. (useful to treat each file independantly)
  Bool_t isMultiChain = kFALSE;
  //if(isApplying_E_Calibration || isApplying_T_Calibration || isEventBuilding) isMultiChain = kTRUE;
  Bool_t Run_process = kFALSE; // In case of multichain is the chain built per Run (kTRUE) or per File (kFALSE)
  if(isEventBuilding) Run_process = kTRUE;
  //Run_process = kFALSE; isMultiChain = kFALSE;// Ponctuellement pour les Tests


  // Defining PID files for detector definition
  TString PIDfilename    = "./PID/";PIDfilename += folder_name; PIDfilename += "sample.pid";
  TString PIDlistname = "./PID/";PIDlistname += folder_name; PIDlistname += "ID.dat";
  TString Typelistname = "./PID/";Typelistname += folder_name; Typelistname += "Type.dat";
  // Definition of the name of file containing the run list
  TString listefilename = "./Parameter_Files/";listefilename += folder_name;
  listefilename += Runlist_filename;
  experiment->SetRUNlistname(listefilename);
  // Open the runlist file.
  std::ifstream runlist_file(listefilename.Data());
  if (!runlist_file.is_open()) {
    std::cerr << "Error opening runlist file: " << listefilename.Data() << std::endl;
    return 1;
  }

  std::vector<std::string> fileEntries;
  std::string line;
  // Read the file line by line.
  while (std::getline(runlist_file, line)) {
    if (!line.empty()) {
      fileEntries.push_back(line);
    }
  }
  runlist_file.close();

  // Check that there is at least a third entry.
  if (fileEntries.size() < 4) {
    std::cerr << "The runlist file does not contain enough entries." << std::endl;
    return 1;
  }

  // The third entry (index 2) is taken as the PARIS_Angles file path.
  TString PARIS_Angles(fileEntries[3].c_str());
  experiment->SetPARISAngle_filename(PARIS_Angles);

  std::cout << "PARIS_Angles file set to: " << PARIS_Angles.Data() << std::endl;
  //TString PARIS_Angles = "./Parameter_Files/"; PARIS_Angles+=folder_name;PARIS_Angles+="PARIS_Angles.txt";
  experiment->SetPIDfilename(PIDfilename);
  experiment->SetPIDlistname(PIDlistname);
  experiment->SetTypelistname(Typelistname);
  //experiment->SetPARISAngle_filename(PARIS_Angles);

  // Defining the detectors
  int check = experiment->SetNbrofdetectors(experiment->GetPIDfilename());
  if(check == 0)
  {
    cout << SetBOLD << SetForeRED << endl;
    cout << " Problem to get the right number of detectors" << endl;
    cout << " Check file " << experiment->GetPIDfilename() << endl,
    cout << RESETTEXT << endl;
    return 0;
  }


  // Now I know how many of each type of detectors I'll need.
  // I create them
  // For the moment are supported detector types:
  // Ge      associated to a internal label (not related to FASTER)   1;
  // BGO     associated to a internal label (not related to FASTER)   2;
  // LaBr    associated to a internal label (not related to FASTER)   3;
  // Cathode associated to a internal label (not related to FASTER)   4;
  // Anode   associated to a internal label (not related to FASTER)   5;
  // TETRA   associated to a internal label (not related to FASTER)   6;
  // OPSA    associated to a internal label (not related to FASTER)   7;
  // IC      associated to a internal label (not related to FASTER)   9;
  // PARIS   associated to a internal label (not related to FASTER)  10;
  // RF      associated to a internal label (not related to FASTER) 250;
  
  // For building Detectors I set a PARIS variable that will require PSD analysis
  experiment->SetThereIsPARIS(kTRUE);
  
  // I build the detectors within the experiment
  // Normally this load the detector parameters (angles, calibrations, PSD informations, ...)
  int checkbuilding = experiment->SetDetectors(experiment->GetPIDfilename(),experiment->GetPIDlistname(),experiment->GetTypelistname());
  if(checkbuilding == 1) {
    cout << SetBOLD << SetForeGRN << endl << "---------------------------------------" << endl;
    cout << "Detectors are built: READY for analysis" << endl;
    cout << "---------------------------------------" << RESETTEXT << endl << endl;
  }
  if(checkbuilding == 0) {
    cout << SetBOLD << SetForeRED << endl << "---------------------------------------" << endl;
    cout << "Detectors are NOT BUILT Properly: STOP analysis" << endl;
    cout << "---------------------------------------" << RESETTEXT << endl << endl;
    return(0);
  }

  // For multiplicity counting later I define a number of modules that can be fired
  // In the context of nu-ball2 campaign:
  experiment->SetNbrofModules(96);

  // Uncomment to check all the detectors properties
  //experiment->PrintallDetectorsInfo();
  //experiment->PrintLabel2Detnbrs();

  // Now I'm ready to look for the data
  
  // Now the experiment is set.
  // I should start the analysis
  // First I need to load the run I want to analyse.
  // I define the name of the TTree I'll have in the Files
  experiment->SetDataTreeName("DataTree");

  // Multichain is for multithreading the data analysis
  check = experiment->ChainBuilder(isMultiChain,Run_process);
  if(check == 0)
  {
    cout << SetBOLD << SetForeRED << endl;
    cout << " Problem to get the Chain built properly" << endl;
    cout << " Check data file(s) " << endl,
    cout << RESETTEXT << endl;
    return 0;
  }

  // Printing of chronometer measurement
  Configurationtimer.Stop();
  Double_t rtime = Configurationtimer.RealTime();
  Double_t ctime = Configurationtimer.CpuTime();

  cout << endl;
  cout << "To configure the experiment it took : " << endl;
  cout << "# RealTime=" << rtime << " seconds, CpuTime="<< ctime << " seconds" <<endl;
  cout << endl;

  // Now I have the data prepared let's do something with it!
  // Chronometer to benchmark the graphes constitution process
  Analysistimer.Reset();
  Analysistimer.Start();
  cout << endl << "----------------------------------------------" << endl;
  cout <<         "Ready to process the data of the experiment >>" << endl;


  // First let's get the timing Right (I need eneryg calibrated spectrum for 60Co coinc)
  // I need to determine the time shift regarding the reference detector
  // Just to make sure I define the reference detector before any time manipulation in data
  experiment->SetReferenceDetector(22);
  if(isCalculating_T_Calibration)
  {
    std::cout<<FOREYEL<<SetBOLD<<"isCalculating_T_Calibration"<<endl;
    experiment->SetReferenceDetector(22);
    
    Bool_t isCalibrated = kTRUE;
    //DrawTimeShifts_ECal(*experiment,-200.,200.); // Begin/End of Time Window in ns
    //DrawTimeShifts(*experiment,-20.,20., "filename"); // Begin/End of Time Window in ns
    //DrawTimeShifts_fissionevents(*experiment,-800.,800.);

    check = CalculateTimealignementShifts(*experiment,isCalibrated);
    /*if(check == 0)
    {
      cout << SetBOLD << SetForeRED << endl;
      cout << " Problem to read the Time calibration run properly" << endl;
      cout << " Check data file(s) " << endl,
      cout << RESETTEXT << endl;
      return 0;
    }*/
  }

  // Then let's apply the time shifts (extracted from calibrated spectrum for 60Co coinc)
  if(isApplying_T_Calibration)
  {
    std::cout<<FOREYEL<<SetBOLD<<"isApplying_T_Calibration"<<endl;
    Bool_t isCalibrated = kTRUE;
    TString CalibrationFileName = experiment->GetTimeCalibration_filename();
    cout << "Using file " << CalibrationFileName << " to time align detectors" << endl;
    experiment->LoadTimeCalibration(CalibrationFileName);
    check = TimeAlignator(*experiment, isCalibrated);
    if(check == 0)
    {
      cout << SetBOLD << SetForeRED << endl;
      cout << " Problem to read the Time calibration run properly" << endl;
      cout << " Check data file(s) " << endl,
      cout << RESETTEXT << endl;
      return 0;
    }
  }

  if(isChecking_T_Calib)
  {
    std::cout<<FOREYEL<<SetBOLD<<"isChecking_T_Calib"<<endl;
    experiment->SetReferenceDetector(22);
    //experiment->SetReferenceDetector(1);
    CheckTimeShifts(*experiment,-20.,20.); // Begin/End of Time Window in ns
    Bool_t isCalibrated = kFALSE;
    //check = CheckTimealignementShifts(*experiment,isCalibrated);
    // if(check == 0)
    // {
    //   cout << SetBOLD << SetForeRED << endl;
    //   cout << " Problem to read the Time calibration run properly" << endl;
    //   cout << " Check data file(s) " << endl,
    //   cout << RESETTEXT << endl;
    //   return 0;
    // }
  }

  // Generally I need to calibrate the data
  // First Let's plot all the calibration spectra
  if(isPlottingSpectra)
  {
    std::cout<<FOREYEL<<SetBOLD<<"isPlottingSpectra"<<endl;
    // 1st step to get the PARIS Angles and "PSD" parameters
    // Make sure to run on 60Co source for this
    //check = DrawAllParisUncalibratedSpectra(*experiment);
    
    // If the function before has been used to determine the PARIS rotation angles then you can run the following on your calibration data
    //check = DrawAllParisUncalibratedSpectra_with_rotation(*experiment);
    //TO apply calibration as well
    //TString CalibrationFileName = experiment->GetNRJCalibration_filename();
    //cout << "Using file " << CalibrationFileName << " to calibrate detectors" << endl;
    //CalibrationFileName += "Calibrations/";
    //CalibrationFileName += "Calibrations_Calib_Ge.data";
    //experiment->LoadCalibration(CalibrationFileName);


    check = DrawAllParisCalibratedSpectra(*experiment);

    //check =  DrawAllCalibrationSpectra(*experiment);
    // After Getting PARIS info, I can plot all the calibration spectra
    //Create a Ecal tree

    //check = ApplyMyEnergyCalibration(*experiment);
    
    // Ideally a 60Co, 22Na, Cs, and AmBe+Ni run should be good for PARIS
    // 152Eu, 60Co, AmBe+Ni should be good for HPGe as it covers from 121 keV to ~ 9MeV
    //check = DrawAllCalibrationSpectra(*experiment);//DrawAllParisUncalibratedSpectra(*experiment);
    if(check == 0)
    {
      cout << SetBOLD << SetForeRED << endl;
      cout << " Problem to plot all the PARIS Calibration spectra properly" << endl;
      cout << " Check data file(s) " << endl,
      cout << RESETTEXT << endl;
      return 0;
    }

  }
  // Second let's get an energy calibration
  if(isCalculating_E_Calibration)
  {
    std::cout<<SetBOLD<<FOREYEL<<"isCalculating_E_Calibration"<<endl;
    check = CalculateCalibrationCoefficients(*experiment,
                                   "Eu",    // Describe the source
                                   kFALSE,  // Calibrate BGO?
                                   kTRUE,   // Calibrate Ge?
                                   kFALSE); // Calibrate LaBr?

    //check = DrawAllParisUncalibratedSpectra(*experiment);
    if(check == 0)
    {
      cout << SetBOLD << SetForeRED << endl;
      cout << " Problem to calculate the Calibration coefficients properly" << endl;
      cout << " Check data file(s) " << endl,
      cout << RESETTEXT << endl;
      return 0;
    }
  }

  // // Now I should have all the coefficients. Let's apply them

  if(isApplying_E_Calibration)
  {
    std::cout<<FOREYEL<<SetBOLD<<"isApplying_E_Calibration"<<endl;
    //TString CalibrationFileName = experiment->GetNRJCalibration_filename();
    //cout << "Using file " << CalibrationFileName << " to calibrate detectors" << endl;
    //CalibrationFileName += "Calibrations/";
    //CalibrationFileName += "Calibrations_Calib_Ge.data";
    //check = experiment->LoadCalibration(CalibrationFileName);
    //if(check == 0)
    //{
      //cout << SetBOLD << SetForeRED << endl;
      //cout << " Problem to get the Calibration loaded properly" << endl;
      //cout << " Check data file(s) " << endl,
      //cout << RESETTEXT << endl;
      //return 0;
    //}

    // I Apply the Calibrations to the data
    //check = EnergyCalibrator(*experiment);
    check = ApplyMyEnergyCalibration(*experiment);
    if(check == 0)
    {
      cout << SetBOLD << SetForeRED << endl;
      cout << " Problem to apply the Calibration properly" << endl;
      cout << " Check data file(s) " << endl,
      cout << RESETTEXT << endl;
      return 0;
    }
  }
  
  // Verifying Calibration
  if(isChecking_E_Calib)
  {
    std::cout<<FOREYEL<<SetBOLD<<"isChecking_E_Calib"<<endl;
    check = DrawAllCalibrationSpectra(*experiment);
    //check = DrawAllParisCalibratedSpectra(*experiment);
    //check = DrawAllEnergyCalibratedSpectra(*experiment,0,2000,4000);
    //check = DrawOneDetectorTypeEnergyCalibratedSpectra(*experiment,"Ge",0,2000,4000);
    //check = DrawOneDetectorTypeEnergyCalibratedSpectra(*experiment,"PARIS",0,2000,4000);
    // // Now let's check the quality on Calibration
    // check = CalculateResidues(*experiment,
    //                         "Co",    // Describe the source
    //                         kFALSE,  // Calibrate BGO?
    //                         kTRUE,   // Calibrate Ge?
    //                         kFALSE);  // Calibrate LaBr?
    if(check == 0)
    {
      cout << SetBOLD << SetForeRED << endl;
      cout << " Problem to read the Calibrated calibration run properly" << endl;
      cout << " Check data file(s) " << endl,
      cout << RESETTEXT << endl;
      return 0;
    }
  }


/*
  // Now I reconstruct physics events
  if(isEventBuilding)
  {
    // Now let's check the quality on Calibration
    //check = ICBuilder(*experiment,650.); // time is given in ns
    check = GammaCoinc(*experiment,60.,2);
    if(check == 0)
    {
      cout << SetBOLD << SetForeRED << endl;
      cout << " Problem to do the event building" << endl;
      cout << " Check data file(s) and or reconstruction " << endl,
      cout << RESETTEXT << endl;
      return 0;
    }
  }

  // Now I read the final data
  if(isEventAnalysing)
  {
    Double_t Emin = 30;
    Double_t Emax = 2000;
    Int_t nbins = Emax*2;
    //DrawAllEnergyCalibratedSpectra(*experiment, Emin, Emax, nbins);
    cout << "MEGA WHOAA" << endl;
  }
*/
  
  // Printing of chronometer measurement
  Analysistimer.Stop();
  rtime = Analysistimer.RealTime();
  ctime = Analysistimer.CpuTime();
  cout << endl;
  cout << "To analyse the data it took : " << endl;
  cout << "# RealTime=" << rtime << " seconds, CpuTime="<< ctime << " seconds" <<endl;
  cout << endl;

  // Printing of the end of the program
  cout << "Exiting program" << endl;
  cout << endl;
  cout << endl;

  // Cleaning Memory
  delete experiment;
  
  // End of the routine
  return (0);

}
