// Call for the principal parts of the functions
#include "CExperiment.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <vector>

// Definition of the I/O function
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CExperiment::CExperiment( TString name ):
nthreads(1),
Experiencename(name),
PIDfilename("UNKNW"),PIDlistname("UNKNW"), Runlist_filename("UNKNW"),Parameter_folder("UNKNW"),
NRJCalibration_filename("UNKNW"),TimeCalibration_filename("UNKNW"),Timewindows_filename("UNKNW"),PARISangles_filename("UNKNW"),
FileDirectory_IN("UNKNW"),FileDirectory_OUT("UNKNW"),DataTreeName("UNKNW"),
Nbrofdetectors(0), NbrofADC(0), NbrofQDC(0), NbrofRF(0),NbrofTChains(0),
isMultiChain(kFALSE),maxLabel(0),
isQDC2(kFALSE),isQDC3(kFALSE),isQDC4(kFALSE), isPARISin(kFALSE)
{
  // Defining the name of the experiment
  Experiencename = name;

  // Initializing a map between detector type and an int for easier switch
  SetDetTypeValues();

  cout << "Analysing Data from " << Experiencename << " Experiment" << endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CExperiment::~CExperiment(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CExperiment::SetDetTypeValues()
{
  // Defining a dictionary between the name of the type of detectors and an Int label
  DetType2nbr.insert(pair<TString, int>("Ge", 1));
  DetType2nbr.insert(pair<TString, int>("BGO", 2));
  DetType2nbr.insert(pair<TString, int>("LaBr", 3));
  DetType2nbr.insert(pair<TString, int>("Cathode", 4));
  DetType2nbr.insert(pair<TString, int>("Anode", 5));
  DetType2nbr.insert(pair<TString, int>("TETRA", 6));
  DetType2nbr.insert(pair<TString, int>("OPSA", 7));
  DetType2nbr.insert(pair<TString, int>("IC", 9));
  DetType2nbr.insert(pair<TString, int>("PARIS", 10));
  DetType2nbr.insert(pair<TString, int>("RF", 250));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CExperiment::SetDetName2Detnbrs()
{

  // Je check que les detecteurs existent
  if(tab_Detectors.size()==0){cout << FORERED << "Error, Detectors have not been set yet" << RESETTEXT << endl;}

  else
  {
    // Je boucle sur tous les detecteurs
    for (int detnbr = 0; detnbr < (int)tab_Detectors.size(); detnbr++)
    {
      //cout << "Associating " << tab_Detectors.at(detnbr)->GetDetectorName() << " with number " << detnbr << endl;
      DetName2nbr.insert(pair<TString, int>(tab_Detectors.at(detnbr)->GetDetectorName(), detnbr));
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int CExperiment::SetNbrofdetectors(const TString name )
{
  string buffer;
  int it1;

  // Reading of the PID list
  cout << endl << "---------------------------------------" << endl;
  cout << "Reading file " << name << endl;
  ifstream detcounter(name,ios::in) ;
  if ( !detcounter.good())
  {
    cout << SetBOLD << SetForeRED << endl;
    cout << "Problem in " << name << endl ;
    cout << RESETTEXT << endl;
    //cout << BOLD(FRED("Problem in" << name.Data()))  << endl ;
    return 0 ;
  }

  // Now I loop on lines
  while(!detcounter.eof())
  {
    getline(detcounter,buffer);
    if(buffer != "") Nbrofdetectors++;
    TString inputfile = buffer.c_str() ;

    // Skipping the line starting with #
    if ( inputfile.Index("#") == -1)
    {
      // Searching for the first : that separate the FASTER label from the channel type
      it1 = inputfile.Index(":",1,inputfile.kExact);
      TString bufftype = inputfile(it1+1,inputfile.Length());
      it1 = bufftype.Index(":",1,inputfile.kExact);

      // Getting the first letter of the detector type
      TString buff = bufftype (0,it1);
      TString comp=bufftype(0);
      if(comp == "A" || comp == "T" || comp == "C") {NbrofADC++;} //---------------------------------- ADC
      if(comp == "Q") {NbrofQDC++;}                               //---------------------------------- QDC
      if(comp == "R") {NbrofRF++;}                                //---------------------------------- RF
      if(comp == "S") {NbrofSAMPLER++;}                           //---------------------------------- SAMPLER
    }
  }

  Nbrofdetectors = Nbrofdetectors/2;
  cout << "We will have " << Nbrofdetectors << " detectors to consider" << endl;
  cout << "Including " << endl;
  cout << NbrofADC << " ADC; " << NbrofQDC << " QDC; " << NbrofRF << " RF; "  << NbrofSAMPLER << " SAMPLER" <<  endl;
  cout << "---------------------------------------" << endl;
  cout << endl;

  return 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int CExperiment::SetDetectors( const TString inputfilename, const TString outputfilename, const TString Typefilename)
{
  int nbrline = 0;
  int localADCnbr = 0;
  int localQDCnbr = 0;
  int localRFnbr = 0;
  int localSAMPLERnbr = 0;
  Int_t detindex_main = 0;
  Int_t chargenumber = 0;
  int maxchargenumber = 0;
  int it1 = 0;
  string buffer_main;

  cout << "Now loading detectors informations : " << endl;

  std::vector<string> names;string name;
  Float_t pos_labr, sigma_labr, pos_nai, sigma_nai,theta, cos, sin, tan, resA, respower, caliba, calibb;
  std::vector<Float_t> vec_pos_labr; std::vector<Float_t> vec_pos_nai;
  std::vector<Float_t> vec_sigma_labr; std::vector<Float_t> vec_sigma_nai;
  std::vector<Float_t> vec_theta;
  std::vector<Float_t> vec_theta_tan;

  // Turn to false if you don't want resolution dependent binning
  Bool_t resolution_binning = true;
  std::vector<Float_t> vec_resA;
  std::vector<Float_t> vec_respower;

  //To apply energy calibration
  Bool_t calibrating = true;
  std::vector<Float_t> vec_caliba;
  std::vector<Float_t> vec_calibb;

  if(isPARISin)
  {
    // I need to store the PARIS data now
    ifstream PARISangleliste(PARISangles_filename,ios::in);
    cout << "Loading PARIS parameters from " << PARISangles_filename << endl;
    if ( !PARISangleliste.good())
    {
      cout << SetBOLD << SetForeRED << endl;
      cout << "Problem in PARIS Angle file named: " << PARISangles_filename << endl ;
      cout << RESETTEXT << endl;
      return 0;
    }

    else
    {
      // Now reading the PARIS informations
      getline(PARISangleliste,buffer_main); // Je skip la premiere ligne

      while(!PARISangleliste.eof())
      { if (resolution_binning== true && calibrating == true){
        PARISangleliste >> name >> pos_labr >> sigma_labr >> pos_nai >> sigma_nai >> theta >> cos >> sin >> tan >> resA >> respower >> caliba >> calibb;
        names.push_back(name);vec_pos_labr.push_back(pos_labr);vec_sigma_labr.push_back(sigma_labr);vec_pos_nai.push_back(pos_nai);vec_sigma_nai.push_back(sigma_nai); vec_theta.push_back(theta); vec_theta_tan.push_back(tan); vec_resA.push_back(resA); vec_respower.push_back(respower);
        vec_caliba.push_back(caliba); vec_calibb.push_back(calibb);
        }
        else{
          PARISangleliste >> name >> pos_labr >> sigma_labr >> pos_nai >> sigma_nai >> theta >> cos >> sin >> tan;
          names.push_back(name);vec_pos_labr.push_back(pos_labr);vec_sigma_labr.push_back(sigma_labr);vec_pos_nai.push_back(pos_nai);vec_sigma_nai.push_back(sigma_nai); vec_theta.push_back(theta); vec_theta_tan.push_back(tan);
          //cout << "\t" << name << "\t" << pos_labr << "\t" << sigma_labr << "\t" << pos_nai << "\t" << sigma_nai << "\t" << theta << "\t" << cos << "\t" << sin << "\t" << tan << endl;
        }
      }

    }
  }

  ifstream PIDliste(inputfilename,ios::in) ;
  if ( !PIDliste.good())
  {
    cout << SetBOLD << SetForeRED << endl;
    cout << "Problem in " << PIDfilename << endl ;
    cout << RESETTEXT << endl;
    return 0;
  }

  // Preparing the file ID.dat as an output
  ofstream iddat(outputfilename,ios::out);

  // Now I loop on lines
  while(!PIDliste.eof())
  {
    getline(PIDliste,buffer_main);
    if(buffer_main != "") nbrline++;

    // Little function to skip a line
    TString inputfile = buffer_main.c_str() ;
    if ( inputfile.Index("#") == -1)
    {
      // cout << inputfile << endl;
      it1 = inputfile.Index(":",1,inputfile.kExact);
      TString buff = inputfile(0,it1);
      detindex_main = buff.Atoi();
      TString bufftype = inputfile(it1+1,inputfile.Length());
      it1 = bufftype.Index(":",1,inputfile.kExact);
      TString Detname = bufftype(it1+1,inputfile.Length());
      //cout << "# Case of detector named " << Detname << " labeled " << detindex_main << endl ;

      // the type ADC, QDC, RF
      TString comp=bufftype(0);
      //cout << comp << endl;
      //cout << comp << endl;
      if(comp == "A" || comp == "T" || comp == "C")//---------------------------------- ADC
      {
        Bool_t isbgo_data = false;
        //dettype = 1;
        // cout << "New ADC # " << localADCnbr << endl;
        // cout << "Detname " << Detname << endl;
        // cout << "detindex_main " << detindex_main << endl;
        // cout << "localADCnbr " << localADCnbr << endl;

        // Now I compute a new general tag for the position
        it1=Detname.Index("_",1,Detname.kExact);
        TString ring=Detname(it1-3,it1-3);
        TString alveoletest = Detname(it1-2,it1-2);
        TString alveole = alveoletest(0);
        if(alveole=="A")alveole= Detname(it1-1,it1-1);
        else if(alveole !="A"){alveole= Detname(it1-2,it1-1);ring=Detname(it1-4,it1-4);}
        TString color = Detname(it1+1,Detname.Length());
        int ringnbr=ring.Atoi();
        int alveolenbr=alveole.Atoi();
        int colornbr=0;
        int type = 0;
        if(color == "red") colornbr = 1;
        else if(color == "black") colornbr = 2;
        else if(color == "green") colornbr = 3;
        else if (color == "blue") colornbr = 4;
        else if (color =="GE") colornbr = 0;
        else if (color =="BGO") {type = 1;colornbr=0;isbgo_data=true;}
        else if (color =="BGO1") {type = 1;colornbr=1;isbgo_data=true;}
        else if (color =="BGO2") {type = 1;colornbr=2;isbgo_data=true;}

        //cout << "isbgo " << isbgo_data << endl;

        // Defining the IDDat for future needs
        int channeltag = ringnbr*10000+alveolenbr*100+type*10+colornbr;
        //cout << detindex_main << " \t " << channeltag << endl;
        iddat << detindex_main << " \t " << channeltag << endl;

        // Creating my ADC and adding it to the ADC list
        ADC *temp_ADC = new ADC(Detname,detindex_main);// = new ADC();
        temp_ADC->SetDetectorName(Detname);
        if(detindex_main > maxLabel) maxLabel = detindex_main;
        temp_ADC->SetDetectorlabel(detindex_main);
        temp_ADC->SetDetectorCARDType("ADC");
        if(!isbgo_data)temp_ADC->SetDetectorType("Ge");
        else temp_ADC->SetDetectorType("BGO");
        temp_ADC->SetADCIndex(localADCnbr);
        temp_ADC->SetCalib(0.,0.,1.,0.);
        temp_ADC->SetIsBGO(isbgo_data);
        temp_ADC->SetRing(ringnbr);
        temp_ADC->SetAlveole(alveolenbr);
        temp_ADC->SetPID(channeltag);
        temp_ADC->SetADCPosition(channeltag);

        tab_Detectors.push_back(temp_ADC);
        localADCnbr++;
      }
      if(comp == "Q")//---------------------------------- QDC with charge number gates
      {
        Bool_t isbgo_data = false;

        //cout << "QDC # " << qdctag[0] << endl;
        //dettype = 2;
        TString charge =bufftype(3);
        chargenumber =charge.Atoi();
        if (chargenumber > maxchargenumber) maxchargenumber = chargenumber;
        

        // Now I compute a new general tag for the position
        it1=inputfile.Index("_",1,inputfile.kExact);
        TString ring= inputfile(it1-3,it1-3);
        TString alveoletest = inputfile(it1-2,it1-2);
        TString alveole = alveoletest(0);
        //cout << alveole.Data() << endl;
        if(alveole=="A")alveole= inputfile(it1-1,it1-1);
        else if(alveole !="A"){alveole= inputfile(it1-2,it1-1);ring=inputfile(it1-4,it1-4);}
        TString color = inputfile(it1+1,inputfile.Length());
        TString color2 =inputfile(inputfile.Length(),inputfile.Length());
        int ringnbr=ring.Atoi();
        int alveolenbr=alveole.Atoi();
        int colornbr = 0;
        int type=0.;
        if(color == "LaBr1") colornbr = 1;
        if(color == "LaBr2") colornbr = 2;
        else if (color =="BGO") {type = 1;colornbr=0;isbgo_data=true;}
        else if (color =="BGO1") {type = 1;colornbr=1;isbgo_data=true;}
        else if (color =="BGO2") {type = 1;colornbr=2;isbgo_data=true;}

        int channeltag = ringnbr*10000+alveolenbr*100+type*10+colornbr;
        iddat << detindex_main << " \t " << channeltag << endl;

        QDC *temp_QDC = new QDC(Detname,detindex_main);// = new ADC();
        temp_QDC->SetDetectorName(Detname);
        if(detindex_main > maxLabel) maxLabel = detindex_main;
        temp_QDC->SetDetectorlabel(detindex_main);
        temp_QDC->SetDetectorCARDType("QDC");
        if(!isbgo_data) temp_QDC->SetDetectorType("LaBr");
        else if(isbgo_data) temp_QDC->SetDetectorType("BGO");
        temp_QDC->SetQDCIndex(localQDCnbr);
        temp_QDC->SetChargeNbr(chargenumber);
        temp_QDC->SetCalib(0.,0.,1.,0.);
        temp_QDC->SetRing(ringnbr);
        temp_QDC->SetAlveole(alveolenbr);
        temp_QDC->SetPID(channeltag);
        temp_QDC->SetQDCPosition(channeltag);

        // As I have PARIS in scanning all the PARIS Names to set the PSD parameters and rotation angle
        if(vec_theta.size() != 0)
        {
          for(int paris=0; paris < static_cast<int>(names.size()); paris++)
          {
            if(names.at(paris) == Detname)
            {
              temp_QDC->SetRotationAngle(vec_theta.at(paris));
              temp_QDC->SetRotationAngleTan(vec_theta_tan.at(paris));
              temp_QDC->SetLaBrDiscriPosition(vec_pos_labr.at(paris));
              temp_QDC->SetLaBrDiscriSigma(vec_sigma_labr.at(paris));
              temp_QDC->SetNaIDiscriPosition(vec_pos_nai.at(paris));
              temp_QDC->SetNaIDiscriSigma(vec_sigma_nai.at(paris));
              if (resolution_binning && calibrating){
                temp_QDC->SetResA(vec_resA.at(paris));
                temp_QDC->SetRespower(vec_respower.at(paris));
                temp_QDC->SetCaliba(vec_caliba.at(paris));
                temp_QDC->SetCalibb(vec_calibb.at(paris));
              }
            }
          }
        }
        
        //temp_QDC->Print();

        //tab_QDC.push_back(temp_QDC);
        tab_Detectors.push_back(temp_QDC);
        //tab_Detectors.back()->PrintDetInfo();
        localQDCnbr++;
      }
      if(comp == "R")//---------------------------------- RF
      {
        //dettype = 3;
        RF *temp_RF = new RF(Detname,detindex_main);
        temp_RF->SetDetectorName(Detname);
        if(detindex_main > maxLabel) maxLabel = detindex_main;
        temp_RF->SetDetectorlabel(detindex_main);
        temp_RF->SetDetectorCARDType("RF");
        temp_RF->SetDetectorType("RF");
        //temp_RF->SetRFIndex(1);
        iddat << detindex_main << " \t " << detindex_main << endl;

        //tab_RF.push_back(temp_RF);
        tab_Detectors.push_back(temp_RF);
        localRFnbr++;
      }
      if(comp == "S")//---------------------------------- SAMPLER
      {
        int channeltag = detindex_main;

        iddat << detindex_main << " \t " << channeltag << endl;

        SAMPLER *temp_SAMPLER = new SAMPLER(Detname,detindex_main);// = new ADC();
        temp_SAMPLER->SetDetectorName(Detname);
        if(detindex_main > maxLabel) maxLabel = detindex_main;
        temp_SAMPLER->SetDetectorlabel(detindex_main);
        temp_SAMPLER->SetDetectorCARDType("SAMPLER");
        temp_SAMPLER->SetDetectorType("SAMPLER");
        temp_SAMPLER->SetSAMPLERIndex(localSAMPLERnbr);
        temp_SAMPLER->SetChargeNbr(chargenumber);
        temp_SAMPLER->SetCalib(0.,0.,1.,0.);
        temp_SAMPLER->SetRing(0);
        temp_SAMPLER->SetAlveole(0);
        temp_SAMPLER->SetPID(channeltag);
        temp_SAMPLER->SetSAMPLERPosition(channeltag);

        //tab_QDC.push_back(temp_QDC);
        tab_Detectors.push_back(temp_SAMPLER);

        localSAMPLERnbr++;
      }
      continue ;
    }
  }

  ifstream Typeliste(Typefilename,ios::in) ;
  if ( !Typeliste.good())
  {
    cout << SetBOLD << SetForeRED << endl;
    cout << "Problem in " << Typefilename << endl ;
    cout << RESETTEXT << endl;
    return 0;
  }
  // Now I loop on lines
  SetLabel2Detnbrs();
  int locallabel;
  TString localtype;
  while(!Typeliste.eof())
  {
    Typeliste >> locallabel >> localtype;
    tab_Detectors.at(Label2Detnbr[locallabel])->SetDetectorType(localtype);
  }

  cout << "Building map of Detector names and numbers" << endl;
  SetDetName2Detnbrs();
  cout << "DONE" << endl;

  bool ADC_OK     = kFALSE;
  bool QDC_OK     = kFALSE;
  bool RF_OK      = kFALSE;
  bool SAMPLER_OK = kFALSE;
  if(localADCnbr == NbrofADC)           ADC_OK = kTRUE;
  if(localQDCnbr == NbrofQDC)           QDC_OK = kTRUE;
  if(localRFnbr  == NbrofRF)            RF_OK  = kTRUE;
  if(localSAMPLERnbr  == NbrofSAMPLER)  SAMPLER_OK  = kTRUE;
  if(ADC_OK && QDC_OK && RF_OK && SAMPLER_OK)
  {
    cout << "We registered " << tab_Detectors.size() << " Detectors" << endl;
    SetQDCNbr(maxchargenumber);
    return 1;
  }
  else
  {
    cout << "Some detectors were not properly initialized " << endl;
    cout << NbrofADC << " ADC were expected and " << localADCnbr << " were initialized" << endl;
    cout << NbrofQDC << " QDC were expected and " << localQDCnbr << " were initialized" << endl;
    cout << NbrofRF << " RF were expected and " << localRFnbr << " were initialized" << endl;
    cout << NbrofSAMPLER << " SAMPLER were expected and " << localSAMPLERnbr << " were initialized" << endl;
    return 0;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CExperiment::SetQDCNbr(const int &nbr)
{
  switch(static_cast<int>(nbr))
  {
    case 2:
    cout << " The data will get 2 energy labels as QDC2 are present in the experiment" << endl;
    isQDC2 = kTRUE;
    break;

    case 3:
    cout << " The data will get 3 energy labels as QDC3 are present in the experiment" << endl;
    isQDC2 = kTRUE;
    isQDC3 = kTRUE;
    break;

    case 4:
    cout << " The data will get 4 energy labels as QDC4 are present in the experiment" << endl;
    isQDC2 = kTRUE;
    isQDC3 = kTRUE;
    isQDC4 = kTRUE;
    break;

    case 1:
    //cout << " The data will get 1 energy labels as QDC1 are present in the experiment" << endl;
    isQDC2 = kFALSE;
    isQDC3 = kFALSE;
    isQDC4 = kFALSE;
    break;

    default:
    cout << "Weird Charge number has been defined. By default only 1 energy label will be considered for this QDC" << endl;
    isQDC2 = kFALSE;
    isQDC3 = kFALSE;
    isQDC4 = kFALSE;

  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CExperiment::SetReferenceDetector(const Int_t &label)
{
  cout << "For the Time analysis the Reference detector is " << endl;
  Ref_detector = tab_Detectors.at(Label2Detnbr[label]);
  Ref_detector->SetIsReferenceDetector(kTRUE);
  tab_Detectors.at(Label2Detnbr[label])->SetIsReferenceDetector(kTRUE);
  Ref_detector->PrintDetInfo();
  //tab_Detectors.at(Label2Detnbr[label])->PrintDetInfo();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int CExperiment::LoadCalibration(TString CalibrationFileName)
{
  TString det_name(0),buffa,buffb,buffc,buffcste,buffamp,buffpower;
  Double_t a3 = 0.;
  Double_t a2 = 0.;
  Double_t a1 = 0.;
  Double_t a0 = 0.;
  Double_t cste = 0.;
  Double_t amp = 0.;
  Double_t power = 0.;
  Bool_t firstparis = true;
  cout << "Now loading detectors informations : " << endl;

  ifstream Calibliste(CalibrationFileName,ios::in) ;
  if ( !Calibliste.good())
  {
    cout << SetBOLD << SetForeRED << endl;
    cout << "Problem in Calibration File named: " << CalibrationFileName << endl ;
    cout << RESETTEXT << endl;
    return 0;
  }


  // Now I loop on lines
  Calibliste >> det_name >> buffa >> buffb >> buffc >> buffcste;// >> buffamp >> buffpower;
  // cout << "Reading the following from calibration file: " << endl;
  // cout << det_name  << "\t" <<  buffa  << "\t" <<  buffb  << "\t" <<  buffc  << "\t" <<  buffd << endl;;
  while(!Calibliste.eof())
  {
    Calibliste >> det_name >> a3 >> a2 >> a1 >> a0;  // >> cste >> amp >> power ;

    TString shortname = det_name(0,5);

    //cout << det_name << "\t" << shortname << "\t" << a1 << "\t" << a0 << endl;
    // cout << "DetName2nbr["<< det_name << "])=" << DetName2nbr[det_name] << endl;

    if(shortname != "PARIS")GetDetector(DetName2nbr[det_name])->SetCalib(0,a2,a1,a0);
    else
    {
      GetDetector(DetName2nbr[det_name])->SetCalib(0,a2,a1,a0);
      if(det_name.EndsWith("_Qs")) {GetDetector(DetName2nbr[det_name])->SetCalib(0,a2,a1,a0);}
      else if(det_name.EndsWith("_Ql")) {GetDetector(DetName2nbr[det_name])->SetSecondaryCalib(0,a2,a1,a0);}
    }
    cout << det_name << "\t" << GetDetector(DetName2nbr[det_name])->GetDetectorName() << "\t" << GetDetector(DetName2nbr[det_name])->GetCalibx3() << "\t" << GetDetector(DetName2nbr[det_name])->GetCalibx2() << "\t" << GetDetector(DetName2nbr[det_name])->GetCalibx1() << "\t" << GetDetector(DetName2nbr[det_name])->GetCalibx0() << endl;
    // int det_label = GetDetector(DetName2nbr[det_name])->GetDetectorlabel();
    // cout << "Det Label " << det_label << "\t"<< GetDetector(Label2Detnbr[det_label])->GetCaliba() << "\t" << GetDetector(Label2Detnbr[det_label])->GetCalibb() << "\t" << GetDetector(Label2Detnbr[det_label])->GetCalibc() << "\t" << GetDetector(Label2Detnbr[det_label])->GetCalibd() << endl;
  }

  return 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int CExperiment::LoadTimeCalibration(TString CalibrationFileName)
{
  Int_t det_label(0);
  Double_t dt(0.),rez(0.);
  cout << "Now loading detectors informations : " << endl;

  ifstream Calibliste(CalibrationFileName,ios::in) ;
  if ( !Calibliste.good())
  {
    cout << SetBOLD << SetForeRED << endl;
    cout << "Problem in " << CalibrationFileName << " for Time alignement"<< endl ;
    cout << RESETTEXT << endl;
    return 0;
  }

  // Now I loop on lines
  while(!Calibliste.eof())
  {
    Calibliste >> det_label >> dt >> rez;
    GetDetector(Label2Detnbr[det_label])->SetTimeShift(dt*1000.);
    //cout << "Det Label " << det_label << "\t"<< GetDetector(Label2Detnbr[det_label])->GetTimeShift() << endl;
  }

  return 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CExperiment::PrintallDetectorsInfo()
{
  for(int i = 0; i < (int)tab_Detectors.size(); i++)
  {
    cout << endl;
    cout << "----------------------------" << endl;
    tab_Detectors.at(i)->PrintDetInfo();
    cout << "----------------------------" << endl;
    cout << endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CExperiment::SetRUNlistname(const TString &name )
{
  cout << endl << "---------------------------------------------------------------------" << endl;
  cout << "Reading " << name << " file to prepare analysis" << endl;
  cout << "---------------------------------------------------------------------" << endl << endl;
  Runlist_filename = name;

  string buffer_main;

  // Reading of the run list
  ifstream liste(Runlist_filename,ios::in) ;
  if ( !liste.good())
  {
    cout << SetBOLD << SetForeRED << endl;
    cout << "Problem in " << Runlist_filename << endl ;
    cout << RESETTEXT << endl;
    return;
  }

  cout<<"#######################################################################" <<endl ;
  cout <<"#     Parameters for the Graphes Creation : "<<endl <<"#"<<endl ;

  getline(liste,buffer_main);
  Parameter_folder = buffer_main;
  cout<<"# All parameter files will be in : "<< Parameter_folder.Data() << endl;

  getline(liste,buffer_main);
  NRJCalibration_filename = buffer_main ;
  cout<<"# Will use calibration file : "<< NRJCalibration_filename.Data() << endl;

  getline(liste,buffer_main);
  TimeCalibration_filename = buffer_main ;
  cout<<"# Will use time gate for coincidence file : "<< TimeCalibration_filename.Data() << endl ;

  getline(liste,buffer_main);
  PARISangles_filename = buffer_main ;
  cout<<"# Will use PARIS Angles file : "<< PARISangles_filename.Data() << endl ;

  getline(liste,buffer_main);
  FileDirectory_IN = buffer_main ;
  cout << "# The root files are in : " << FileDirectory_IN << endl << "#" << endl ;

  getline(liste,buffer_main);
  FileDirectory_OUT = buffer_main ;
  cout << "# The root output files will be in : " << FileDirectory_OUT << endl << "#" << endl ;
  cout << "#######################################################################" <<endl ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CExperiment::SetLabel2Detnbrs()
{
  cout << "Max Label = "  << maxLabel << endl;
  // I fill the vector with bullshit numbers
  for (auto i = 0; i < maxLabel+1; i++) Label2Detnbr.push_back(-666);

  for(int i = 0; i < (int)tab_Detectors.size(); i++)
  {
    // I change the value associated to the label
    Label2Detnbr.at(tab_Detectors.at(i)->GetDetectorlabel()) = i;
    cout<<i<<endl;
  }
}
/*void CExperiment::SetLabel2Detnbrs()
{
    cout << "Max Label = " << maxLabel << endl;

    // Ensure Label2Detnbr is initialized with -1 (for invalid entries)
    Label2Detnbr.assign(maxLabel + 1, -1);  // Initialize with -1

    for (int i = 0; i < static_cast<int>(tab_Detectors.size()); i++)
    {
        int label = tab_Detectors.at(i)->GetDetectorlabel();
        if (label >= 0 && label <= maxLabel)  // Ensure label is within bounds
        {
            cout << "Processing detector index: " << i << ", label: " << label << endl;
            Label2Detnbr.at(tab_Detectors.at(i)->GetDetectorlabel()) = i;
            cout << "Assigned detector index " << i << " to Label2Detnbr[" << label << "]" << endl;
        }
        else
        {
            cout << "Skipping invalid label " << label << " for detector index " << i << endl;
        }
    }
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CExperiment::PrintLabel2Detnbrs()
{
  // for (const auto &[k, v] : Label2Detnbr)
  // {
  //     std::cout <<"Label " << k << "  is Detector " << v << "\n";
  // }
  cout<<"all good "<<endl;
  for (auto i=0; i<(int)Label2Detnbr.size(); i++) std::cout <<"Label " << i << "  is Detector " << Label2Detnbr.at(i) << "\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Int_t CExperiment::ChainBuilder(Bool_t isMultiChain, Bool_t per_Run)
{
  // I want one unique TChain for running the same stuff on all Data
  if(!isMultiChain)
  {
    SetNbrofTChains(1);
    string buffer_main;
    std::string inputfile;
    std::string outputname;
    TChain *localTChain = new TChain(GetDataTreeName().data());

    cout << endl << "---------------------------------------" << endl;
    cout << "Reading the Runs that should belong to the Chain" << endl;

    // Reading of the run list
    if(Runlist_filename == "UNKNW")
    {
      cout << SetBOLD << SetForeRED << endl << "ERROR : File name to get the list of run has nos been defined " << endl;
      cout << RESETTEXT << endl;
      return 0;
    }
    ifstream liste(Runlist_filename,ios::in);

    // Skipping the first 4 lines that have already been read.
    getline(liste,buffer_main);
    getline(liste,buffer_main);
    getline(liste,buffer_main);
    getline(liste,buffer_main);
    getline(liste,buffer_main);
    getline(liste,buffer_main);

    // Now looping on all the runs that should be included in the TChain
    while(!liste.eof())
    {
      // Reading the name of the run that should be added
      getline(liste,buffer_main);cout << buffer_main << endl;
      // Checking it is not an empty line
      if(buffer_main != "")
      {
        inputfile = buffer_main.c_str() ;
        // Little function to skip a file if it has been commented with a *
        if ( inputfile.find("#") != std::string::npos)
        {
          //cout << "# Skiping the " << inputfile << endl ;
          continue ;
        }
        inputfile = FileDirectory_IN;
        inputfile += buffer_main;
        cout << "# Inclusion of " << inputfile << " to the TChain >>" << endl;
        localTChain -> Add(inputfile.data());
        Datafilenames.push_back(inputfile.data());
        outputname = buffer_main;
      }
    }
    chained_oak.push_back(localTChain);
  }
  // I want several TChains to send to various threads
  else
  {
    // First I need to know how many TChains I'll have
    SetNbrofTChains(ChainCounter(per_Run));
    cout << endl << "---------------------------------------" << endl;
    cout << "We will build " << GetNbrofTChains() << "  TChains "  << endl;

    string buffer_main;
    TString inputfile;
    TString chainfile;
    TString outputname;
    int runct(0), runnbr(0);
    int memory_runbr(0);
    if(DataTreeName == "UNKNW")
    {
      cout << SetBOLD << SetForeRED << endl << "ERROR : Tree name to get the TChain has not been defined " << endl;
      cout << RESETTEXT << endl;
      return 0;
    }
    TChain **local_chained_oak;
    if(NbrofTChains != 0)
    {
      local_chained_oak = new TChain*[NbrofTChains];
      for(auto co=0; co < NbrofTChains; co++)
      {
        local_chained_oak[co] =  new TChain(DataTreeName.data());
      }
    }
    else
    {
      cout << SetBOLD << SetForeRED << endl << "ERROR : The number of TChain has not been properly defined " << endl;
      cout << RESETTEXT << endl;
      return 0;
    }

    // Creation of the TChain to process all the data
    cout << "Distributing run files in different chains" << endl;
    // Reading of the run list
    if(Runlist_filename == "UNKNW")
    {
      cout << SetBOLD << SetForeRED << endl << "ERROR : File name to get the list of run has not been defined " << endl;
      cout << RESETTEXT << endl;
      return 0;
    }

    ifstream liste(Runlist_filename,ios::in);
    getline(liste,buffer_main);
    getline(liste,buffer_main);
    getline(liste,buffer_main);
    getline(liste,buffer_main);
    getline(liste,buffer_main);
    getline(liste,buffer_main);

    while(!liste.eof())
    {
      // Vreation of the TChain
      getline(liste,buffer_main);
      if(buffer_main != "")
      {
        // Little function to skip a file
        inputfile = buffer_main.c_str() ;
        if ( inputfile.Index("#") != -1)
        {
          //cout << "# Skiping the " << inputfile << endl ;
          continue ;
        }
        if(per_Run)
        {
          int it1=inputfile.Index("n",1,inputfile.kExact);
          int it2=inputfile.Index("_",1,inputfile.kExact);
          TString run = inputfile(0,it1+1);
          if(run == "Run" || run =="run")
          {
            TString srunnbr = inputfile(it1+1,it2);
            runnbr = srunnbr.Atoi();
          }
          if(runnbr == memory_runbr)
          {
            chainfile = GetFileDirectory_IN();
            chainfile += inputfile;
            local_chained_oak[runct] -> Add(chainfile);
            local_chained_oak[runct] -> SetName(chainfile);
          }

          if(runnbr != memory_runbr)
          {
            // I place the TChain in the adequate chain
            runct++;
            chainfile = GetFileDirectory_IN();
            chainfile += inputfile;
            local_chained_oak[runct] -> Add(chainfile);
            local_chained_oak[runct] -> SetName(chainfile);
            memory_runbr = runnbr;
          }
        }
        else
        {

          chainfile = GetFileDirectory_IN();
          chainfile += inputfile;
          local_chained_oak[runct] -> Add(chainfile);
          local_chained_oak[runct] -> SetName(chainfile);
          runct++;
        }
      }
    }

    // Now I should have all the different TChains and Associated outpufile names
    // I put them into the experiment adequate variable
    //cout << chained_oak.size() << endl;
    for(auto co=0; co < NbrofTChains; co++)
    {
      chained_oak.push_back(local_chained_oak[co]);
    }
  }
  Bool_t Check_ChainNbr = kFALSE;
  if((int)chained_oak.size() == GetNbrofTChains()) Check_ChainNbr = kTRUE;
  if(Check_ChainNbr)
  {
    cout << "We registered " << chained_oak.size() << " TChains properly" << endl;
    cout << "---------------------------------------" << endl;
    return 1;
  }
  else
  {
    cout << SetForeCYN << endl << "                WARNING                    " << endl;
    cout << "Some TChains were not properly initialized " << endl;
    cout << GetNbrofTChains() << " TChain were expected and " << (int)chained_oak.size() << " were initialized" << endl;
    cout << RESETTEXT << endl;
    return 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Int_t CExperiment::ChainCounter(Bool_t per_Run)
{
  string buffer_main;
  TString inputfile;
  TString outputname;
  int runct(1), runnbr(0);
  int memory_runbr(0);

  // Reading of the run list
  if(Runlist_filename == "UNKNW")
  {
    cout << "ERROR : File name to get the list of run has nos been defined " << endl;
    return 0;
  }
  ifstream liste(Runlist_filename,ios::in);

  // Skipping the first 4 lines that have already been read.
  getline(liste,buffer_main);
  getline(liste,buffer_main);
  getline(liste,buffer_main);
  getline(liste,buffer_main);
  getline(liste,buffer_main);
  getline(liste,buffer_main);

  // Now looping on all the runs that should be included in the TChain
  while(!liste.eof())
  {
    // Reading the name of the run that should be added
    getline(liste,buffer_main);
    if(buffer_main != "")
    {
      // Little function to skip a file
      inputfile = buffer_main.c_str() ;
      if ( inputfile.Index("#") != -1)
      {
        //cout << "# Skiping the " << inputfile << endl ;
        continue ;
      }
      int it1=inputfile.Index("n",1,inputfile.kExact);
      int it2=inputfile.Index("_",1,inputfile.kExact);
      TString run = inputfile(0,it1+1);

      if(per_Run)
      {
        if(run == "Run" || run =="run")
        {
          TString srunnbr = inputfile(it1+1,it2);
          runnbr = srunnbr.Atoi();
        }

        if(runnbr != memory_runbr)
        {
          runct++;
          memory_runbr = runnbr;
        }
      }
      else
      {
        runct++;
      }
    }
  }
  return runct-1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//ClassImp(CExperiment);
