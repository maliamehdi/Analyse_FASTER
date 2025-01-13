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

// Call for the principal parts of the functions
#include "CAnalyse.h"

// Definition of the I/O function
using namespace std;
using namespace ROOT;

// Definition of the standard types
// using label_type = UInt_t;
// using tm_Caltype = Double_t;
// using tm_Uncaltype = ULong64_t;
// using nrj_Ucaltype = UInt_t;
// using nrj_Caltype = Double_t;
// using pu_type = Bool_t;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//          Peak CORRELATED FUNCTIONS
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#define PEAKS_C_FIT_AREAS 1 // To get the peak areas and not just the height
// Definition of a Gaussian with a linear backgrd
Double_t fpeaks(Double_t *x, Double_t *par, Int_t npeaks) {
   Double_t result = par[0] + par[1]*x[0];
   for (Int_t p=0;p<npeaks;p++) {
      Double_t norm  = par[3*p+2]; // "height" or "area"
      Double_t mean  = par[3*p+3];
      Double_t sigma = par[3*p+4];
#if defined(__PEAKS_C_FIT_AREAS__)
      norm /= sigma * (TMath::Sqrt(TMath::TwoPi())); // "area"
#endif /* defined(__PEAKS_C_FIT_AREAS__) */
      result += norm*TMath::Gaus(x[0],mean,sigma);
   }
   return result;
}

Double_t fpeaks2(Double_t *x, Double_t *par, Int_t npeaks) {
   Double_t result = 0.1;
   for (Int_t p=0;p<npeaks;p++) {
      Double_t norm   = par[5*p+0];
      Double_t mean1  = par[5*p+1];
      Double_t sigma1 = par[5*p+2];
      Double_t mean2  = par[5*p+3];
      Double_t sigma2 = par[5*p+4];
      result += norm*TMath::Gaus(x[0],mean1,sigma1)*TMath::Gaus(x[1],mean2,sigma2);
   }
   return result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t DoubleTailedStepedGaussian(Double_t *xx,Double_t*pp)
{
  Double_t f_tot = 0.;

  Double_t Back_const = pp[1];
  Double_t Back_slope = pp[2];
  Double_t Back_Exp = pp[3];

  f_tot += (Back_const + (xx[0])*Back_slope)*exp((xx[0])*Back_Exp);

  Double_t Ampli     = pp[4];
  Double_t Mean      = pp[5];
  Double_t Sigma     = pp[6]*1./sqrt(8.*log(2.));
  Double_t Lambda    = pp[7];
  Double_t Rho       = pp[8];
  Double_t S         = pp[9];

  Double_t U         = (xx[0]-Mean)/Sigma;
  Double_t f_g       = Ampli*TMath::Exp(-U*U*0.5);
  Double_t f_lambda  = Ampli*TMath::Exp(-0.5*Lambda*(2.*U-Lambda));
  Double_t f_rho     = Ampli*TMath::Exp(-0.5*Rho*(2.*U-Rho));
  Double_t f_S       = Ampli*S*1./((1+TMath::Exp(U))*(1+TMath::Exp(U)));

  if(U<Lambda) f_tot += f_lambda;
  else if(U>Rho) f_tot += f_rho;
  else f_tot += f_g;

  f_tot += f_S;

  return f_tot;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//          ENERGY CORRELATED FUNCTIONS
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t Faster2bitsNRJConverter(nrj_Rawtype enrj, const CExperiment &experiment, Int_t detectorindex)
{
  Double_t alea = rand()*1./(RAND_MAX*1.);
  Double_t NRJ = (enrj/(Double_t)experiment.GetDetector(static_cast<int>(detectorindex))->GetMaxchNumber())*(Double_t)experiment.GetDetector((int)detectorindex)->GetNbrChannels()+alea;

  return NRJ;
}

UInt_t CompressFASTERValue(const nrj_Rawtype &enrj, const int &originalBits, const int &targetBits)
{
  
  if (targetBits >= originalBits) {
        std::cerr << "Le nombre de bits cible doit être inférieur au nombre de bits d'origine." << std::endl;
        return enrj;
    }
  
  // Calcul de la mise à l'échelle
    uint32_t maxValueOriginal = (1 << originalBits) - 1; // Valeur maximale avec originalBits bits
    uint32_t maxValueTarget = (1 << targetBits) - 1; // Valeur maximale avec targetBits bits

    // Initialisation du générateur de nombres aléatoires pour le dithering
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    // Ajouter du dithering
    double dither = dis(gen) - 0.5; // Valeur aléatoire entre -0.5 et 0.5

    // Mise à l'échelle de la valeur avec dithering
    double scaledValue = static_cast<double>(enrj) * maxValueTarget / maxValueOriginal + dither;

    // Clamping pour s'assurer que la valeur est dans les limites
    if (scaledValue < 0.0) scaledValue = 0.0;
    if (scaledValue > maxValueTarget) scaledValue = maxValueTarget;

    uint32_t compressedValue = static_cast<uint32_t>(scaledValue);

    return compressedValue;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int CalculateCalibrationCoefficients(const CExperiment &experiment, TString sourcename, Bool_t Calib_BGO, Bool_t Calib_Ge, Bool_t Calib_LaBr)
{
  //Definition du germe pour le tirage aleatoire
  srand48(time(NULL));

  std::vector<TH1F*> NRJSpectra;

  // Define all the TH1F for the analysis
  TString outputfilename = experiment.GetFileDirectory_OUT();
  outputfilename += "UncalibratedEnergyspectra_all.root";
  TFile *spectrafile = new TFile(outputfilename,"READ");
  if(!spectrafile || spectrafile->IsZombie() || spectrafile->GetNkeys()==0){
    cout << "File with all uncalibrated spectra do not exists; creating it ..." << endl;
    NRJSpectra = DrawAllEnergyUncalibratedSpectra(experiment);
    cout << "Number of spectra to analyse = " << NRJSpectra.size() << endl;
  }
  //NRJpectra = DrawAllEnergyUncalibratedSpectra(experiment);
  else
  {
    cout << FORECYN << "File already exists" << endl;
    cout << "Loading spectra from file..." << endl;
    for(auto sindex = 0; sindex < (int)experiment.GetDetectors().size();sindex++)
    {
      TString spectrumname = "nrjspectrum";
      spectrumname+=experiment.GetDetectors().at(sindex)->GetDetectorName();
      TH1F *temp_h1 = (TH1F*)spectrafile->Get(spectrumname);
      NRJSpectra.push_back(temp_h1);
    }
    cout << "Spectra Loaded ..." << RESETTEXT << endl;
    cout << "Number of spectra to analyse = " << NRJSpectra.size() << endl;
  }

  if(NRJSpectra.size() == 0)
  {
    cout << FORERED << SetBOLD << endl;
    cout << "Problem with Uncalibrated Energy spectra generation" << endl;
    cout << RESETTEXT << endl;
    return 0;
  }

  // I define the energies I need to calibrate
  std::vector<Double_t> Ref_nrj;
  std::vector<Double_t> Ref_nrj_err;
  std::vector<Double_t> Ref_nrj_intensity;

  tie(Ref_nrj,Ref_nrj_err,Ref_nrj_intensity) = SetCalibrationSource(sourcename,Calib_BGO,Calib_Ge,Calib_LaBr);

  // I declare two files.
  // One to store the TGraph with the linear fit..
  TString output_fit_name = experiment.GetFileDirectory_OUT();
  std::vector<TChain *> tab_chained_oak = experiment.GettheTChain();
  TChain *chained_oak = tab_chained_oak.at(0);
  TString temp = chained_oak->GetFile()->GetName();
  int it1 = temp.Index("Calibration_",12,1,temp.kExact);temp = temp(it1,temp.Length());
  it1 = temp.Index(".",1,temp.kExact);temp = temp(0,it1);
  output_fit_name += temp;
  output_fit_name +="_AllFits";
  if(Calib_Ge) output_fit_name +="_Ge.root";
  if(Calib_BGO) output_fit_name +="_BGO.root";
  if(Calib_LaBr) output_fit_name +="_LaBr.root";
  cout << FOREGRN << "The Calibration TGraph will be saved in " << output_fit_name << endl;
  TFile *output_fit = new TFile(output_fit_name,"RECREATE");
  TString calibfilename=experiment.GetFileDirectory_OUT();

  calibfilename += temp;//experiment.GetNRJCalibration_filename();
  if(Calib_Ge) calibfilename +="_Calib_Ge";
  if(Calib_BGO) calibfilename +="_Calib_BGO";
  if(Calib_LaBr) calibfilename +="_Calib_LaBr";
  calibfilename +=".data";
  cout << "The Calibration parameters will be saved in " << calibfilename << RESETTEXT<< endl;
  ofstream calibfile(calibfilename,ios::out);
  calibfile.precision(7);

  // Now I proceed to the Analysis of each spectrum
  // I will MT the process
  // Now I proceed to the Analysis of each spectrum
    int bouik = NRJSpectra.size();
    // I define the energies I need to calibrate
    for(int sindex = 0; sindex < (int)bouik; sindex++)
    {
      Bool_t DoIcalib = kFALSE;

      // I look at a detector type
      if(experiment.GetDetectors().at(sindex)->GetDetectorType() == "BGO"  && Calib_BGO)  DoIcalib = kTRUE;
      if(experiment.GetDetectors().at(sindex)->GetDetectorType() == "Ge"   && Calib_Ge)   DoIcalib = kTRUE;
      if(experiment.GetDetectors().at(sindex)->GetDetectorType() == "LaBr" && Calib_LaBr) DoIcalib = kTRUE;

      if(NRJSpectra.at(sindex)->GetEntries() < 15000) {DoIcalib = kFALSE;}

      if(DoIcalib)
      {
        cout << RESETTEXT << FOREGRN << endl << "Calibrating Detector : " << experiment.GetDetectors().at(sindex)->GetDetectorName() << RESETTEXT << endl;
        // First I define a new spectrum
        TSpectrum *s = new TSpectrum(Ref_nrj.size()*4); // I do not allow to find more than four times the number of reference peaks

        // I try to determinate the background with compton edges
        Int_t nbins = NRJSpectra.at(sindex)->GetXaxis()->GetNbins();
        Double_t xmin  = 0;
        Double_t xmax  = NRJSpectra.at(sindex)->GetXaxis()->GetXmax();
        Double_t * source = new Double_t[nbins];
        Double_t * dest = new Double_t[nbins];
        //Double_t sigma = 6;
        for (int i = 0; i < nbins; i++) source[i]=NRJSpectra.at(sindex)->GetBinContent(i + 1);

        // Now I search for the background
        //int niteration = 10;
        if(Calib_LaBr)s->Background(source,nbins,50,TSpectrum::kBackDecreasingWindow,
                      TSpectrum::kBackOrder4,kFALSE,
                      TSpectrum::kBackSmoothing15,kFALSE);
        if(Calib_Ge) s->Background(source,nbins,30,TSpectrum::kBackDecreasingWindow,
                      TSpectrum::kBackOrder2,kTRUE,
                      TSpectrum::kBackSmoothing15,kTRUE);

        // I prepare a TH1 to store the background
        TString name = NRJSpectra.at(sindex)->GetName();name +="_1";
        TH1F* myspectrum = (TH1F*)NRJSpectra.at(sindex)->Clone(name);
        TH1F *d1 = new TH1F("d1","",nbins,xmin,xmax);

        // I fill the spectrum witht the background
        for (int i = 0; i < nbins; i++) d1->SetBinContent(i + 1,source[i]);

        // Now I clean the spectrum from the background
        myspectrum->Add(myspectrum,d1,1,-1);
        //for (int i = 0; i < nbins; i++) source[i]=NRJSpectra.at(sindex)->GetBinContent(i + 1);


        Int_t nfound(0);
        if(experiment.GetDetector(sindex)->GetDetectorType() == "LaBr" && NRJSpectra.at(sindex)->GetEntries() > 1000)
        {
          Double_t intensity = 0.05;
          while(nfound < (int)Ref_nrj.size())
          {
            nfound = s->Search(NRJSpectra.at(sindex),8,"",intensity);
            intensity = intensity - 0.01;
          }
        }
        if(experiment.GetDetector(sindex)->GetDetectorType() == "BGO" && NRJSpectra.at(sindex)->GetEntries() > 1000 ) nfound = s->Search(NRJSpectra.at(sindex),40,"",0.01);
        if(experiment.GetDetector(sindex)->GetDetectorType() == "Ge"  && NRJSpectra.at(sindex)->GetEntries() > 1000 )
        {
          Double_t intensity = 0.05;
          while(nfound < (int)Ref_nrj.size() && intensity > 0.01)
          {
            nfound = s->Search(NRJSpectra.at(sindex),24,"",intensity);//s->SearchHighRes(source, dest, nbins, 20, intensity, kTRUE, 2, kTRUE, 15);//Search(NRJSpectra.at(sindex),4,"",intensity);//
            intensity = intensity - 0.01;
          }
        }

        if (nfound != 0 || nfound < (int)Ref_nrj.size()) {cout <<"Found " << nfound << " candidate peaks to fit" << endl;DoIcalib= kFALSE;}
        else
        {
          cout << BACKRED << FOREWHT << "No peaks were found in the spectrum" << endl;
          if(NRJSpectra.at(sindex)->GetEntries() == 0) cout << "The spectrum is empty for detector " << experiment.GetDetectors().at(sindex)->GetDetectorName() << endl;
          else cout << "A problem occured with detector " << experiment.GetDetectors().at(sindex)->GetDetectorName() <<BACKBLK <<  RESETTEXT << endl;
        }

        // I loop on the peaks
        // I load from TSpectrum the peak list
        Double_t *xpeaks = s->GetPositionX();
        Double_t *ypeaks = s->GetPositionY();
        Int_t npeaks = 0;
        Double_t sigma=30.;

        // I declare the Tables I need to fit:
        std::vector<Double_t> weight;
        std::vector<Double_t> pos;
        std::vector<Double_t> pos_err;
        std::vector<Double_t> reso;
        std::vector<Double_t> reso_err;
        std::vector<Double_t> integral;
        std::vector<Double_t> integral_err;
        //std::vector<Double_t> FWTM;
        // I calculate the ratio of the first two energies from the source
        Double_t initial_NRJ_ratio = Ref_nrj.at(1)/Ref_nrj.at(0);

        // I'm gonna fit all the found peaks
        if(Calib_Ge)
        {
          std::vector<Double_t> xpeak;
          std::vector<Double_t> ypeak;

          for(auto p = 0; p<nfound;p++)
          {
            if(xpeaks[p]>2500.)
            {
              //xpeak.push_back(NRJSpectra.at(sindex)->GetXaxis()->GetBinCenter(SimpleRoundTo(xpeaks[p])+1));
              xpeak.push_back(SimpleRoundTo(xpeaks[p]));
              ypeak.push_back(ypeaks[p]);
            }
          }
          // Double_t memory_intense_nrj;
          // if(xpeak.at(0) > 60) memory_intense_nrj =  xpeak.at(0);
          // else memory_intense_nrj =  xpeak.at(1);
          //
          sortVectors(xpeak, less<Double_t>(), xpeak,ypeak);
          printVector(xpeak);
          // Now they are sorted I need to identify the right peaks...
          // Now i search for this energy ratio in the peak that were found:
          int pos_1stpeak=0;int pos_2ndpeak = 0;
          Bool_t leaveloop = kFALSE;
          for(auto first = 0; first < (int)xpeak.size() && !leaveloop;first++)
          {
            for(auto second = first+1; second < (int)xpeak.size()&& !leaveloop;second++)
            {
              Double_t data_ratio = xpeak.at(second)/xpeak.at(first);
              //cout << xpeak.at(first) << "\t" << xpeak.at(second) << "\t" << data_ratio << " for " << initial_NRJ_ratio << " & memory_intense_nrj = " << memory_intense_nrj << endl;

              if(TMath::Abs(initial_NRJ_ratio-data_ratio) < 0.09)
              {
                pos_1stpeak = first;pos_2ndpeak=second;
                leaveloop = kTRUE;
                //cout << "Pair found in " << pos_1stpeak << " and " << pos_2ndpeak << endl;
              }
            }
          }

          // Now I've identified the first peaks, I get the realistic list of peaks
          // If I don't start from the right one I erase lower energy peaks
          if(pos_1stpeak !=0)
          {
            xpeak.erase(xpeak.begin(),xpeak.begin()+pos_1stpeak);
            ypeak.erase(ypeak.begin(),ypeak.begin()+pos_1stpeak);
          }
          printVector(xpeak);

          nfound=xpeak.size();

          for (auto p=0;p<nfound;p++)//nfound;p++)
          {
            // First I get bin center
            Double_t xp = xpeak.at(p);
            //cout << endl << endl << "One peak found @ " << xp << endl;
            std::vector<Double_t> temp_amplitude;
            std::vector<Double_t> temp_pos;
            std::vector<Double_t> temp_sigma;

            // Define my gaussian around my peak
            TF1 *f1 = new TF1("f1", "gaus", xp-(2*sigma), xp+(2*sigma));
            f1->SetParameter(1,xp);
            //f1->SetParLimits(1,xp-sigma,xp+sigma);

            // Then I get amplitude
            Int_t bin = NRJSpectra.at(sindex)->GetXaxis()->FindBin(xp);
            Double_t yp = NRJSpectra.at(sindex)->GetBinContent(bin);
            f1->SetParameter(0,yp);
            //f1->SetParLimits(0,yp-TMath::Sqrt(yp),yp+TMath::Sqrt(yp));

            // I fix the sigma
            f1->SetParameter(2,sigma);
            //f1->SetParLimits(2,0.,(Double_t)3*sigma);

            // Printing for debug
            //cout << xp-(3*sigma) << "\t" << xp+(3*sigma) << endl;
            //cout << xp << " p/m " << xp-sigma << "\t"<<xp+sigma<< endl;
            //cout << yp << " p/m " << yp-TMath::Sqrt(yp) << "\t"<<yp+TMath::Sqrt(yp)<< endl;
            //cout << sigma << " p/m " << 0 << "\t"<<3*sigma<< endl;

            // Then I fit my spectrum with a gaussian
            //cout << "Fitting with f1" << endl;
            NRJSpectra.at(sindex)->Fit("f1","RIQE");

            // I store the first fit info
            temp_amplitude.push_back(f1->GetParameter(0));
            temp_pos.push_back(f1->GetParameter(1));
            weight.push_back(10);
            temp_sigma.push_back(f1->GetParameter(2));

            // Printing the result out
            //cout << "temp Peak amplitude " << temp_amplitude << endl;
            //cout << "temp Peak mean value " << f1->GetParameter(1) << endl;
            //cout << "temp Peak sigma " << temp_sigma << endl;

            // I get a second fit with a linear bckground
            TF1 *f2 = new TF1("f2", "gaus(0)+pol1(3)", xp-(2*sigma), xp+(2*sigma));
            // cout << xp-(3*sigma) << "\t" << xp+(3*sigma) << endl;
            // cout << xp << " p/m " << xp-sigma << "\t"<<xp+sigma<< endl;
            // cout << yp << " p/m " << yp-TMath::Sqrt(yp) << "\t"<<yp+TMath::Sqrt(yp)<< endl;
            // cout << sigma << " p/m " << 0 << "\t"<<3*sigma<< endl << endl << endl;
            f2->SetParameter(0,temp_amplitude.at(0));
            f2->SetParameter(1,temp_pos.at(0));
            f2->SetParameter(2,temp_sigma.at(0));
            // f2->SetParameter(0,yp);
            // f2->SetParameter(1,xp);
            // f2->SetParameter(2,sigma);
            // f2->SetParLimits(0,yp-TMath::Sqrt(yp),yp+TMath::Sqrt(yp));f2->SetParameter(0,temp_amplitude);
            // f2->SetParLimits(1,xp-sigma,xp+sigma);f2->SetParameter(1,temp_pos);
            // f2->SetParLimits(2,0,3*sigma);f2->SetParameter(2,temp_sigma);
            //cout << "Fitting with f2" << endl;
            NRJSpectra.at(sindex)->Fit("f2","RIQE");
            temp_amplitude.push_back(f2->GetParameter(0));
            temp_pos.push_back(f2->GetParameter(1));
            weight.push_back(20);
            temp_sigma.push_back(f2->GetParameter(2));

            //cout << f2->GetParameter(1) << endl;

            // And now I Fit with a skewed gaussian
            // TF1 * f3 = new TF1("sgf","2.*gaus(x,[0],[1],[2])*ROOT::Math::normal_cdf([3]*x,1,0)",xp-(2*sigma), xp+(2*sigma));
            // f3->SetParameter(0,temp_amplitude.at(1));
            // f3->SetParameter(1,temp_pos.at(1));
            // f3->SetParameter(2,temp_sigma.at(1));
            // NRJSpectra.at(sindex)->Fit("sgf","RIQE");
            // temp_amplitude.push_back(f3->GetParameter(0));
            // temp_pos.push_back(f3->GetParameter(1));
            // temp_sigma.push_back(f3->GetParameter(2));

            //cout << "Fitting with f3" << endl;
            TF1 * fFitFunction = new TF1("MyFit",DoubleTailedStepedGaussian,xp-(2*sigma), xp+(2*sigma),10);
            fFitFunction->SetParName(0, "NumberOfPeaks");
            fFitFunction->SetParName(1, "BkgConst");
            fFitFunction->SetParName(2, "BkgSlope");
            fFitFunction->SetParName(3, "BkgExp");
            fFitFunction->SetParName(4+0, "Height");
            fFitFunction->SetParName(4+1, "Position");
            fFitFunction->SetParName(4+2, "FWHM");
            fFitFunction->SetParName(4+3, "LeftTail");
            fFitFunction->SetParName(4+4, "RightTail");
            fFitFunction->SetParName(4+5, "AmplitudeStep");
            fFitFunction->FixParameter(0, 1);
            fFitFunction->SetParameter(1,f2->GetParameter(3));
            fFitFunction->SetParameter(2,f2->GetParameter(4));
            fFitFunction->FixParameter(3, 0.);
            fFitFunction->SetParameter(4,f2->GetParameter(0));
            fFitFunction->SetParameter(5,f2->GetParameter(1));
            fFitFunction->SetParameter(6,f2->GetParameter(2));
            fFitFunction->SetParameter(7,-2.);
            fFitFunction->SetParLimits(7,-5.,-0.1);
            fFitFunction->SetParameter(8,2.);
            fFitFunction->SetParLimits(8,0.1,5);
            fFitFunction->SetParameter(9,0.01);
            fFitFunction->SetParLimits(9,-1.,1.);

            // cout << "Before fitting " << endl;
            // for(auto p= 0; p<10; p++)
            // {
            //   cout << fFitFunction->GetParName(p) << "\t" << fFitFunction->GetParameter(p) << endl;
            // }

            NRJSpectra.at(sindex)->Fit("MyFit","R0Q");
            // cout << "After fitting " << endl;
            // for(auto p= 0; p<10; p++)
            // {
            //   cout << fFitFunction->GetParName(p) << "\t" << fFitFunction->GetParameter(p) << endl;
            // }
            temp_amplitude.push_back(fFitFunction->GetParameter(4));
            temp_pos.push_back(fFitFunction->GetParameter(5));
            weight.push_back(1);
            temp_sigma.push_back(fFitFunction->GetParameter(6));

            //cout << fFitFunction->GetParameter(5) << endl;

            // Now I calculate the final parameters
            Double_t final_amplitude = VectorMean(temp_amplitude);
            Double_t final_pos = VectorWeightedMean(temp_pos,weight);
            Double_t final_pos_err = VectorStdDev(temp_pos);
            Double_t final_sigma = VectorMean(temp_sigma);
            Double_t final_sigma_err = VectorStdDev(temp_sigma);

            // Printing the result out
            // cout << "Peak amplitude " << final_amplitude << endl;
            //cout << "Peak mean value " << final_pos << "p/m " << final_pos_err << endl;
            // cout << "Peak sigma " << final_sigma << "p/m " << final_sigma_err<< endl<< endl << endl;

            // I calculate the integral of the peak
            Double_t final_Int = final_amplitude* final_sigma*TMath::Sqrt(2*TMath::Pi());
            Double_t a = final_amplitude; Double_t c = final_sigma; Double_t Delta_a = TMath::Sqrt(final_amplitude);Double_t Delta_c = final_sigma_err;
            Double_t cons = 2*TMath::Pi();
            Double_t final_Int_syst_err = TMath::Sqrt(TMath::Power(c,2)*TMath::Power(Delta_a,2)*cons+TMath::Power(a,2)*TMath::Power(Delta_c,2)*cons);
            Double_t final_Int_err = TMath::Sqrt(TMath::Power(final_Int_syst_err,2)+final_Int);

            // I fill the sortVectors
            pos.push_back(final_pos);
            pos_err.push_back(final_pos_err);
            reso.push_back(final_sigma*2.35482);
            reso_err.push_back(final_sigma_err*2.35482);
            integral.push_back(final_Int);
            integral_err.push_back(final_Int_err);

            // I stored all the peaks... Now I have to assign them
            // cout << " and fitted @ " << final_pos << " with a difference of " << TMath::Abs(xp-final_pos) << endl;
            npeaks++;

          }
          d1->Delete();
        }

        if(Calib_LaBr)
        {
          std::vector<Double_t> xpeak;
          std::vector<Double_t> ypeak;

          for(auto p = 0; p<nfound;p++)
          {
            xpeak.push_back(xpeaks[p]);
            ypeak.push_back(ypeaks[p]);
          }
          Double_t memory_intense_nrj;
          if(xpeak.at(0) > 60) memory_intense_nrj =  xpeak.at(0);
          else memory_intense_nrj =  xpeak.at(1);
          //printVector(xpeak);
          sortVectors(xpeak, less<Double_t>(), xpeak,ypeak);
          //printVector(xpeak);

          // Now they are sorted I need to identify the right peaks...
          // Now i search for this energy ratio in the peak that were found:
          int pos_1stpeak=0;int pos_2ndpeak=0;
          Bool_t leaveloop = kFALSE;
          for(auto first = 0; first < (int)xpeak.size() && !leaveloop;first++)
          {
            for(auto second = first+1; second < (int)xpeak.size()&& !leaveloop;second++)
            {
              Double_t data_ratio = xpeak.at(second)/xpeak.at(first);
              //cout << xpeak.at(first) << "\t" << xpeak.at(second) << "\t" << data_ratio << " for " << initial_NRJ_ratio << " & memory_intense_nrj = " << memory_intense_nrj << endl;

              if(TMath::Abs(initial_NRJ_ratio-data_ratio) < 0.15 && xpeak.at(first) == memory_intense_nrj)
              {
                pos_1stpeak = first;pos_2ndpeak=second;
                leaveloop = kTRUE;
                //cout << "Pair found in " << pos_1stpeak << " and " << pos_2ndpeak << endl;
              }
            }
          }

          // Now I've identified the first peaks, I get the realistic list of peaks
          // If I don't start from the right one I erase lower energy peaks
          if(pos_1stpeak !=0)
          {
            xpeak.erase(xpeak.begin(),xpeak.begin()+pos_1stpeak);
            ypeak.erase(ypeak.begin(),ypeak.begin()+pos_1stpeak);
          }
          //printVector(xpeak);

          nfound=xpeak.size();
          for (auto p=0;p<nfound;p++)//nfound;p++)
          {
            if(xpeak.at(p)/xpeak.at(0) < 9.91)
            {
              // First I get bin center
              Double_t xp = xpeak.at(p);
              //cout << "One peak found @ " << xp << endl;

              // Define my gaussian around my peak
              TF1 *f1 = new TF1("f1", "gaus", xp-(2*sigma), xp+(2*sigma));
              f1->SetParameter(1,xp);
              //f1->SetParLimits(1,xp-sigma,xp+sigma);

              // Then I get amplitude
              Int_t bin = NRJSpectra.at(sindex)->GetXaxis()->FindBin(xp);
              Double_t yp = NRJSpectra.at(sindex)->GetBinContent(bin);
              f1->SetParameter(0,yp);
              //f1->SetParLimits(0,yp-TMath::Sqrt(yp),yp+TMath::Sqrt(yp));

              // I fix the sigma
              f1->SetParameter(2,sigma);
              //f1->SetParLimits(2,0.,(Double_t)3*sigma);

              // Printing for debug
              //cout << xp-(3*sigma) << "\t" << xp+(3*sigma) << endl;
              //cout << xp << " p/m " << xp-sigma << "\t"<<xp+sigma<< endl;
              //cout << yp << " p/m " << yp-TMath::Sqrt(yp) << "\t"<<yp+TMath::Sqrt(yp)<< endl;
              //cout << sigma << " p/m " << 0 << "\t"<<3*sigma<< endl;

              // Then I fit my spectrum with a gaussian
              NRJSpectra.at(sindex)->Fit("f1","RIQE");

              // I store the first fit info
              Double_t temp_amplitude = f1->GetParameter(0);
              Double_t temp_pos = f1->GetParameter(1);
              Double_t temp_sigma = f1->GetParameter(2);

              // Printing the result out
              // cout << "temp Peak amplitude " << temp_amplitude << endl;
              // cout << "temp Peak mean value " << temp_pos << endl;
              // cout << "temp Peak sigma " << temp_sigma << endl;

              // I get a second fit with a linear bckground
              TF1 *f2 = new TF1("f2", "gaus(0)+pol1(3)", xp-(3*sigma), xp+(3*sigma));
              // cout << xp-(3*sigma) << "\t" << xp+(3*sigma) << endl;
              // cout << xp << " p/m " << xp-sigma << "\t"<<xp+sigma<< endl;
              // cout << yp << " p/m " << yp-TMath::Sqrt(yp) << "\t"<<yp+TMath::Sqrt(yp)<< endl;
              // cout << sigma << " p/m " << 0 << "\t"<<3*sigma<< endl << endl << endl;
              f2->SetParameter(0,temp_amplitude);
              f2->SetParameter(1,temp_pos);
              f2->SetParameter(2,temp_sigma);
              // f2->SetParameter(0,yp);
              // f2->SetParameter(1,xp);
              // f2->SetParameter(2,sigma);
              // f2->SetParLimits(0,yp-TMath::Sqrt(yp),yp+TMath::Sqrt(yp));f2->SetParameter(0,temp_amplitude);
              // f2->SetParLimits(1,xp-sigma,xp+sigma);f2->SetParameter(1,temp_pos);
              // f2->SetParLimits(2,0,3*sigma);f2->SetParameter(2,temp_sigma);
              NRJSpectra.at(sindex)->Fit("f2","RIQE");

              // Now I calculate the final parameters
              Double_t final_amplitude = (f2->GetParameter(0)+temp_amplitude)/2.;
              Double_t final_pos = (f2->GetParameter(1)+temp_pos)/2.;
              Double_t final_pos_err = TMath::Abs(final_pos-temp_pos);
              Double_t final_sigma = (f2->GetParameter(2)+temp_sigma)/2.;
              Double_t final_sigma_err = TMath::Abs(final_sigma-temp_sigma);

              // Printing the result out
              // cout << "Peak amplitude " << final_amplitude << endl;
              // cout << "Peak mean value " << final_pos << "p/m " << final_pos_err << endl;
              // cout << "Peak sigma " << final_sigma << "p/m " << final_sigma_err<< endl<< endl << endl;

              // I calculate the integral of the peak
              Double_t final_Int = final_amplitude* final_sigma*TMath::Sqrt(2*TMath::Pi());
              Double_t a = final_amplitude; Double_t c = final_sigma; Double_t Delta_a = TMath::Sqrt(final_amplitude);Double_t Delta_c = final_sigma_err;
              Double_t cons = 2*TMath::Pi();
              Double_t final_Int_syst_err = TMath::Sqrt(TMath::Power(c,2)*TMath::Power(Delta_a,2)*cons+TMath::Power(a,2)*TMath::Power(Delta_c,2)*cons);
              Double_t final_Int_err = TMath::Sqrt(TMath::Power(final_Int_syst_err,2)+final_Int);

              // I fill the sortVectors
              pos.push_back(final_pos);
              pos_err.push_back(final_pos_err);
              reso.push_back(final_sigma*2.35482);
              reso_err.push_back(final_sigma_err*2.35482);
              integral.push_back(final_Int);
              integral_err.push_back(final_Int_err);

              // I stored all the peaks... Now I have to assign them
              //cout << " and fitted @ " << final_pos << endl;
              npeaks++;
            }
            else // Normally I adress the question of the 1408 keV plus intrinsic activity in case resolution is shitty
            {
              // First I get bin center
              Double_t xp = xpeak.at(p);
              //cout << "One peak found @ " << xp << endl;

              // Define my gaussian around my peak
              TF1 *f1 = new TF1("f1", "gaus", xp-(2*sigma), xp+(2*sigma));
              f1->SetParameter(1,xp);
              //f1->SetParLimits(1,xp-sigma,xp+sigma);

              // Then I get amplitude
              Int_t bin = NRJSpectra.at(sindex)->GetXaxis()->FindBin(xp);
              Double_t yp = NRJSpectra.at(sindex)->GetBinContent(bin);
              f1->SetParameter(0,yp);
              //f1->SetParLimits(0,yp-TMath::Sqrt(yp),yp+TMath::Sqrt(yp));

              // I fix the sigma
              f1->SetParameter(2,sigma);
              //f1->SetParLimits(2,0.,(Double_t)3*sigma);

              // Printing for debug
              // cout << xp-(3*sigma) << "\t" << xp+(3*sigma) << endl;
              // cout << xp << " p/m " << xp-sigma << "\t"<<xp+sigma<< endl;
              // cout << yp << " p/m " << yp-TMath::Sqrt(yp) << "\t"<<yp+TMath::Sqrt(yp)<< endl;
              // cout << sigma << " p/m " << 0 << "\t"<<3*sigma<< endl;

              // Then I fit my spectrum with a gaussian
              NRJSpectra.at(sindex)->Fit("f1","RIQE");

              // I store the first fit info
              Double_t temp_amplitude = f1->GetParameter(0);
              Double_t temp_pos = f1->GetParameter(1);
              Double_t temp_sigma = f1->GetParameter(2);

              // Printing the result out
              // cout << "temp Peak amplitude " << temp_amplitude << endl;
              // cout << "temp Peak mean value " << temp_pos << endl;
              // cout << "temp Peak sigma " << temp_sigma << endl;

              //I get a second fit with a linear bckground
              TF1 *f2 = new TF1("f2", "gaus(0)+pol1(3)", xp-(2*sigma), xp+(2*sigma));
              // cout << xp-(3*sigma) << "\t" << xp+(3*sigma) << endl;
              // cout << xp << " p/m " << xp-sigma << "\t"<<xp+sigma<< endl;
              // cout << yp << " p/m " << yp-TMath::Sqrt(yp) << "\t"<<yp+TMath::Sqrt(yp)<< endl;
              // cout << sigma << " p/m " << 0 << "\t"<<3*sigma<< endl << endl << endl;
              f2->SetParameter(0,temp_amplitude);
              f2->SetParameter(1,temp_pos);
              f2->SetParameter(2,temp_sigma);
              // f2->SetParameter(0,yp);
              // f2->SetParameter(1,xp);
              // f2->SetParameter(2,sigma);
              // f2->SetParLimits(0,yp-TMath::Sqrt(yp),yp+TMath::Sqrt(yp));f2->SetParameter(0,temp_amplitude);
              // f2->SetParLimits(1,xp-sigma,xp+sigma);f2->SetParameter(1,temp_pos);
              // f2->SetParLimits(2,0,3*sigma);f2->SetParameter(2,temp_sigma);
              NRJSpectra.at(sindex)->Fit("f2","RIQE");

              // Now I calculate the final parameters
              Double_t final_amplitude = (f2->GetParameter(0)+temp_amplitude)/2.;
              Double_t final_pos = (f2->GetParameter(1)+temp_pos)/2.;
              Double_t final_pos_err = TMath::Abs(final_pos-temp_pos);
              Double_t final_sigma = (f2->GetParameter(2)+temp_sigma)/2.;
              Double_t final_sigma_err = TMath::Abs(final_sigma-temp_sigma);

              // Printing the result out
              // cout << "Peak amplitude " << final_amplitude << endl;
              // cout << "Peak mean value " << final_pos << "p/m " << final_pos_err << endl;
              // cout << "Peak sigma " << final_sigma << "p/m " << final_sigma_err<< endl<< endl << endl;

              // I calculate the integral of the peak
              Double_t final_Int = final_amplitude* final_sigma*TMath::Sqrt(2*TMath::Pi());
              Double_t a = final_amplitude; Double_t c = final_sigma; Double_t Delta_a = TMath::Sqrt(final_amplitude);Double_t Delta_c = final_sigma_err;
              Double_t cons = 2*TMath::Pi();
              Double_t final_Int_syst_err = TMath::Sqrt(TMath::Power(c,2)*TMath::Power(Delta_a,2)*cons+TMath::Power(a,2)*TMath::Power(Delta_c,2)*cons);
              Double_t final_Int_err = TMath::Sqrt(TMath::Power(final_Int_syst_err,2)+final_Int);

              // I fill the sortVectors
              pos.push_back(final_pos);
              pos_err.push_back(final_pos_err);
              reso.push_back(final_sigma*2.35482);
              reso_err.push_back(final_sigma_err*2.35482);
              integral.push_back(final_Int);
              integral_err.push_back(final_Int_err);

              // I stored all the peaks... Now I have to assign them
              //cout << " and fitted @ " << final_pos << endl;
              npeaks++;
            }
          }
          d1->Delete();
        }

        // Now I search the right peaks
        //printVector(pos);
        std::vector<int> Toremove;
        std::vector<Double_t> NRJratios;
        for(auto index = 1; index < (int)Ref_nrj.size(); index++)
        {
          cout << index << "\t" << Ref_nrj.at(index) << "\t" <<Ref_nrj.at(0) << "\t" << Ref_nrj.at(index)/Ref_nrj.at(0) << endl;
          NRJratios.push_back(Ref_nrj.at(index)/Ref_nrj.at(0));
        }
        int locali = 0;
        for(auto index = 1; index  <  (int)pos.size(); index++)
        {
          Double_t data_ratio = pos.at(index)/pos.at(0);
          //Double_t integral_ratio = integral.at(index)/pos.at(0);
          //if(locali < (int)NRJratios.size())cout << "NRJratios.at("<<locali<<") = " << NRJratios.at(locali) << "\t with " << pos.at(index) << " & " << pos.at(0) << " giving " << data_ratio << " & int ratio = " << integral_ratio << endl;
          //
          // //Because of non-linearity I have to consider a second order term to correct the threshold
          // Double_t Ratio_corrected;
          // if(index > 2 && data_ratio > 9.5)
          // {
          //   Double_t data_ratio0 = pos.at(index-2)/pos.at(0);
          //   Double_t data_ratio_i_1 = pos.at(index-1)/pos.at(0);
          //   Double_t corr_slope = (data_ratio0-data_ratio_i_1)/(pos.at(0)-pos.at(index-1));
          //   Double_t corr_origin = data_ratio0 - corr_slope*pos.at(0);
          //   cout << corr_slope << "\t" << corr_origin << endl;
          //   Ratio_corrected = corr_slope*pos.at(index);
          //   cout << Ratio_corrected << endl;
          //   //data_ratio+=Ratio_corrected;
          // }
          // if(locali < (int)NRJratios.size())cout << "NRJratios.at("<<locali<<") = " << NRJratios.at(locali) << "\t with " << pos.at(index) << " & " << pos.at(0) << " giving after correction: " << data_ratio << endl << endl;
          Double_t precision = 0.26;
          if(Calib_LaBr) precision = 0.35;
          if(locali < (int)NRJratios.size())
          {
            if(TMath::Abs(NRJratios.at(locali)-data_ratio) < precision) locali++;
            else if(Calib_LaBr && NRJratios.at(locali) > 6. && TMath::Abs(1-data_ratio/NRJratios.at(locali)) < 0.12)locali++; // To consider non linearity that builds up at hiher energy
            else Toremove.push_back(index);
          }
          else
          {
            Toremove.push_back(index);
            locali++;
          }
        }

        // Now I remove the useless peaks
        for (int rmv =Toremove.size()-1; rmv > -1; rmv--)
        {
          pos.erase(pos.begin()+Toremove.at(rmv));
          pos_err.erase(pos_err.begin()+Toremove.at(rmv));
          reso.erase(reso.begin()+Toremove.at(rmv));
          reso_err.erase(reso_err.begin()+Toremove.at(rmv));
          integral.erase(integral.begin()+Toremove.at(rmv));
          integral_err.erase(integral_err.begin()+Toremove.at(rmv));
        }

        // I define the TGraph to plot calibration, resolution & efficiency
        TString outcalibmeasure=experiment.GetFileDirectory_OUT();outcalibmeasure+=experiment.GetDetectors().at(sindex)->GetDetectorName();
        outcalibmeasure +="_allpeaks.txt";
        ofstream calibfileperdetect(outcalibmeasure,ios::out);
        calibfileperdetect.precision(7);
        // I print to the file the peak position and errors
        for(int p=0; p < (int)pos.size();p++)
        {
          calibfileperdetect << Ref_nrj.at(p) << "\t"<< pos.at(p) << "\t"<< Ref_nrj_err.at(p) << "\t"<< pos_err.at(p) << endl;
        }
        calibfileperdetect.close();
        printVector(pos);
        printVector(Ref_nrj);
        printVector(pos_err);
        printVector(Ref_nrj_err);
        auto calibration = new TGraphErrors(pos.size(),&(pos[0]),&(Ref_nrj[0]),&(pos_err[0]),&(Ref_nrj_err[0]));
        TString title = "Calibration of "; title += experiment.GetDetectors().at(sindex)->GetDetectorName();
        TString title2 = "Residues for "; title2 += experiment.GetDetectors().at(sindex)->GetDetectorName();
        calibration->SetTitle(title);
        calibration->SetName(title);
        // auto resolution = new TGraphErrors((int)pos.size(),Ref_nrj,reso,Ref_nrj_err,reso_err);
        // auto efficiency = new TGraphErrors((int)pos.size(),Ref_nrj,integral,Ref_nrj_err,integral_err);
        TGraph *residue;

        // Now I fit the calibration
        TF1 *f1_1 = new TF1("f1_1","pol1");
        TF1 *f2_1 = new TF1("f2_1","pol2");
        TF1 *f3_1 = new TF1("f3", "pol3");
        TString fitname = "f_"; fitname += experiment.GetDetector(sindex)->GetDetectorName();
        TF1 *f3 = new TF1(fitname, "pol3");
        TF1 *f4 = new TF1(fitname,"pol1");
        cout << FOREGRN <<  "Peaks identified now fitting calibration curve " << RESETTEXT << endl;
        std::vector<Double_t> residues;
        Double_t a, b, c, d;
        if(Calib_Ge)
        {
          calibration->Fit(f4,"MQW");
          std::cout.precision(7);
          a = 0;
          b = 0;
          c = f4->GetParameter(1);
          d = f4->GetParameter(0);
          calibfile << experiment.GetDetector(sindex)->GetDetectorlabel() << "\t" << a << "\t" << b << "\t" << c << "\t" << d << endl;

          // Now I calculate the Residues
          cout << FOREGRN << "Now Calculating Residues" << RESETTEXT << endl;
          for(auto nrj = 0 ; nrj < (int)pos.size(); nrj++)
          {
            residues.push_back(Ref_nrj.at(nrj)-(c*pos.at(nrj)+d));
          }
        }
        residue = new TGraph(residues.size(),&(Ref_nrj[0]),&(residues[0]));
        residue -> SetTitle(title2);
        residue -> SetName(title2);

        if(Calib_LaBr)
        {
          calibration->Fit(f1_1,"MQ");
          for(auto param = 0; param < 2; param++)f2_1->FixParameter(param,f1_1->GetParameter(param));
          calibration->Fit(f2_1,"MQ");
          for(auto param = 0; param < 3; param++)f2_1->SetParameter(param,f2_1->GetParameter(param));
          calibration->Fit(f2_1,"MQ");
          for(auto param = 0; param < 3; param++)f3_1->FixParameter(param,f2_1->GetParameter(param));
          calibration->Fit(f3_1,"MQ");
          for(auto param = 0; param < 4; param++)f3->SetParameter(param,f3_1->GetParameter(param));
          calibration->Fit(f3,"MQ");
          // cout << "The linear fit result: " << endl;

          // cout << "Energy (kev) = " << f3->GetParameter(1) << "* x + " << f3->GetParameter(0) << endl;
          // cout << "With a probability " << f3->GetProb() << endl;
          calibfile << experiment.GetDetector(sindex)->GetDetectorlabel() << "\t" << f3->GetParameter(3) << "\t" << f3->GetParameter(2) << "\t" << f3->GetParameter(1) << "\t" << f3->GetParameter(0) << endl;
        }
        output_fit->cd();
        calibration->Write();
        residue->Write();
        if(Calib_Ge)f4->Write();
        if(Calib_LaBr)f3->Write();
        cout << endl;
      }
    } // End of the loop on all spectra

  output_fit->Close();
  cout << "Calibration curve and residues have been written to " << output_fit->GetName() << endl;
  cout << "Calibration parameters have been written to " << calibfilename << endl;
  calibfile.close();

  return 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<TH1F*> DrawAllEnergyUncalibratedSpectra(const CExperiment &experiment)
{
  //Definition du germe pour le tirage aleatoire
  srand48(time(NULL));

  //Declaration of all variables for names
  TString title;
  TString spectrumname;

  ROOT::EnableThreadSafety();

  // Declaration of time spectra
  Int_t nbrofspectra = experiment.GetDetectors().size();
  std::vector<TH1F*> NRJspectra;
  Int_t highestdetlabel(0),highestnbrchannels(0),highestEmax(0);
  Int_t nbrchannels(0);
  Int_t Emin(0);
  Int_t Emax(0);

  // Defining all the energy spectra
  for(int sindex = 0; sindex < nbrofspectra; sindex++)
  {
    TH1F* localNRJspectrum;
    // Defining the title and name of the spectrum
    title = "Energy Spectrum of detector ";
    spectrumname = "nrjspectrum";
    title += experiment.GetDetector(sindex)->GetDetectorName();
    spectrumname +=experiment.GetDetector(sindex)->GetDetectorName();
    if(experiment.GetDetector(sindex)->GetDetectorlabel()>highestdetlabel) highestdetlabel = experiment.GetDetectors().at(sindex)->GetDetectorlabel();

    // Getting the right range
    if(experiment.GetDetector(sindex)->GetDetectorType()!="RF")
    {
      nbrchannels = experiment.GetDetector(sindex)->GetNbrChannels();
      //cout << nbrchannels << "\t" << Emax << endl;
      if(nbrchannels > highestnbrchannels) highestnbrchannels=nbrchannels;
      Emax = experiment.GetDetector(sindex)->GetMaxchNumber();//GetMaxchNumber();
      if(Emax > highestEmax) highestEmax=Emax;
      localNRJspectrum = new TH1F(spectrumname,title,nbrchannels,Emin,Emax);

      // Storing the NRJsectrum
      NRJspectra.push_back(localNRJspectrum);
    }

    // Cleaning the names for the next iteration
    spectrumname.Clear();
    title.Clear();
  }

  // I declare la Time matrix to check alignement later
  TH2F* NRJmatrix = new TH2F("NRJalignementmatrix","NRJ spectra of all detectors",highestdetlabel,1,highestdetlabel,highestnbrchannels,Emin,highestEmax);

  // Creation of the file to save all the data
  TString outputfilename = experiment.GetFileDirectory_OUT();
  outputfilename += "UncalibratedEnergyspectra_all.root";
  TFile *outputfile = new TFile(outputfilename,"RECREATE");

  // Loading the TTree for reading the Data
  std::vector<TChain *> tab_chained_oak = experiment.GettheTChain();
  TChain *chained_oak = tab_chained_oak.at(0);

  label_Rawtype index;
  nrj_Rawtype enrj;
  std::cout << FOREBLU << "Loading Tree ..." << std::endl;
  /*auto cachesize = 10000000; // 10 MBytes
  chained_oak -> SetCacheSize(cachesize);
  chained_oak -> SetCacheLearnEntries(100000);*/
  chained_oak -> SetBranchStatus("*",0);
  chained_oak -> SetBranchAddress("label",&index);   // Detector number
  chained_oak -> SetBranchAddress("nrj",&enrj);      // NRJ of the hit
  std::cout << "Loaded .." << RESETTEXT << std::endl;

  // Another Timer
  TStopwatch localTimer;

  //Starting of the chronometer
  localTimer.Reset();
  localTimer.Start();

  //---------------------------------------------------------------------------//
  //                                                                           //
  //                            Reading the TTree                              //
  //                                                                           //
  //---------------------------------------------------------------------------//
  // // Determination of the number of event that have to be treated
  // const auto max_workers = std::thread::hardware_concurrency();
  auto n_workers = experiment.GetThreadsNbr();
  // if(n_workers > max_workers){
  //   cout << "Current machine only supports " << max_workers << " parallel threads. I will be using this number" << endl;
  //   n_workers = max_workers;
  // }
  const auto chainentries = chained_oak -> GetEntries();
  // const auto range = chainentries/n_workers;
  // int percent = (int)(0.05*chainentries);

  std::cout << "Beginning of Energy spectra building research" << std::endl;
  std::cout << "On " << chainentries << " entries" << std::endl;

  // Creating a ProgressBar which is nice in MT mode
  using namespace indicators;
  ProgressBar bar{
    option::BarWidth{60},
    option::Start{"["},
    option::Fill{"="},
    option::Lead{">"},
    option::Remainder{"-"},
    option::End{"]"},
    option::PostfixText{"Reading"},
    option::ShowElapsedTime{true},
    option::ShowRemainingTime{true},
    option::ForegroundColor{Color::grey},
    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
  };

  // I read all the file of the TChain
  int fileindex  = 0;
  int filenumber =  experiment.GetDataFileNames().size();
  for(auto &file : experiment.GetDataFileNames())
  {
    fileindex++;
    //PrintVector(experiment.GetDataFileNames());
    bar.set_progress((int)((double)fileindex/(double)filenumber*100.));
    // Defining a TTreeProcessor object to handle the TTree in the File in a MT mode
    ROOT::TTreeProcessorMT TP(file, experiment.GetDataTreeName(), n_workers);

    // Scanning TTree
    // Launch the parallel processing of the tree
    auto loop_and_fill = [&] (TTreeReader &myReader){
      TTreeReaderValue<label_Rawtype> labelRV(myReader,"label");
      TTreeReaderValue<nrj_Rawtype> enrjRV(myReader,"nrj");


      while(myReader.Next())
      {
        auto label = *labelRV;
        auto NRJ   = *enrjRV;

        int spectrumindex = experiment.GetLabel2Detnbrs(static_cast<int>(label));
        //Double_t NRJ = enrj;//Faster2bitsNRJConverter(enrj, experiment,spectrumindex); // To convert to a sort of integer (it should remove the "non-linearity")

        // Calculating the right energy
        if(NRJ > 100) NRJspectra.at(spectrumindex)->Fill(NRJ);

        // And the energy Matrix
        if (NRJ > 100) NRJmatrix->Fill(label,NRJ);
      }
    };
    TP.Process(loop_and_fill);
  }


  std::cout << endl << RESETTEXT << "Saving in " << outputfilename << std::endl;
  outputfile->cd();

  for(int i = 0; i < (int) NRJspectra.size(); i++) NRJspectra.at(i)->Write();

  NRJmatrix->Write();
  outputfile->Close();

  cout << "NRJ spectra are saved " << endl;

  // Printing of chronometer measurement
  localTimer.Stop();
  Double_t rtime2 = localTimer.RealTime();
  Double_t ctime2 = localTimer.CpuTime();
  std::cout << std::endl;
  std::cout << "End of Energy Spectra Plotting" << std::endl;
  std::cout << "# RealTime=" << rtime2 << " seconds, CpuTime="<< ctime2 << " seconds" <<std::endl;
  std::cout << std::endl;

  return NRJspectra;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int DrawAllEnergyCalibratedSpectra(const CExperiment &experiment, Double_t Emin, Double_t Emax, Int_t binning)
{
  //Definition du germe pour le tirage aleatoire
  srand48(time(NULL));

  //Declaration of all variables for names
  TString title;
  TString spectrumname;

  ROOT::EnableThreadSafety();

  // Declaration of time spectra
  Int_t nbrofspectra = experiment.GetDetectors().size();
  std::vector<TH1F*> NRJspectra;
  Int_t highestdetlabel(0),highestnbrchannels(0),highestEmax(0);
  Int_t nbrchannels(0);
  

  // Defining all the energy spectra
  for(int sindex = 0; sindex < nbrofspectra; sindex++)
  {
    TH1F* localNRJspectrum;
    // Defining the title and name of the spectrum
    title = "Energy Spectrum of detector ";
    spectrumname = "nrjspectrum";
    title += experiment.GetDetector(sindex)->GetDetectorName();
    spectrumname +=experiment.GetDetector(sindex)->GetDetectorName();
    if(experiment.GetDetector(sindex)->GetDetectorlabel()>highestdetlabel) highestdetlabel = experiment.GetDetectors().at(sindex)->GetDetectorlabel();

    // Getting the right range
    if(experiment.GetDetector(sindex)->GetDetectorType()!="RF")
    {
      nbrchannels = experiment.GetDetector(sindex)->GetNbrChannels();
      //cout << nbrchannels << "\t" << Emax << endl;
      if(nbrchannels > highestnbrchannels) highestnbrchannels=nbrchannels;
      //Emax = experiment.GetDetector(sindex)->GetMaxchNumber();//GetMaxchNumber();
      if(Emax > highestEmax) highestEmax=Emax;
      localNRJspectrum = new TH1F(spectrumname,title,nbrchannels,Emin,Emax);

      // Storing the NRJsectrum
      NRJspectra.push_back(localNRJspectrum);
    }

    // Cleaning the names for the next iteration
    spectrumname.Clear();
    title.Clear();
  }

  // I declare la Time matrix to check alignement later
  TH2F* NRJmatrix = new TH2F("NRJalignementmatrix","NRJ spectra of all detectors",highestdetlabel,1,highestdetlabel,highestnbrchannels,Emin,highestEmax);

  // Creation of the file to save all the data
  TString outputfilename = experiment.GetFileDirectory_OUT();
  outputfilename += "CalibratedEnergyspectra_all.root";
  TFile *outputfile = new TFile(outputfilename,"RECREATE");

  // Loading the TTree for reading the Data
  std::vector<TChain *> tab_chained_oak = experiment.GettheTChain();
  TChain *chained_oak = tab_chained_oak.at(0);

  label_Rawtype index;
  nrj_type enrj;
  std::cout << FOREBLU << "Loading Tree ..." << std::endl;
  /*auto cachesize = 10000000; // 10 MBytes
  chained_oak -> SetCacheSize(cachesize);
  chained_oak -> SetCacheLearnEntries(100000);*/
  chained_oak -> SetBranchStatus("*",0);
  chained_oak -> SetBranchAddress("label",&index);   // Detector number
  chained_oak -> SetBranchAddress("nrj",&enrj);      // NRJ of the hit
  std::cout << "Loaded .." << RESETTEXT << std::endl;

  // Another Timer
  TStopwatch localTimer;

  //Starting of the chronometer
  localTimer.Reset();
  localTimer.Start();

  //---------------------------------------------------------------------------//
  //                                                                           //
  //                            Reading the TTree                              //
  //                                                                           //
  //---------------------------------------------------------------------------//
  // // Determination of the number of event that have to be treated
  // const auto max_workers = std::thread::hardware_concurrency();
  auto n_workers = experiment.GetThreadsNbr();
  // if(n_workers > max_workers){
  //   cout << "Current machine only supports " << max_workers << " parallel threads. I will be using this number" << endl;
  //   n_workers = max_workers;
  // }
  const auto chainentries = chained_oak -> GetEntries();
  // const auto range = chainentries/n_workers;
  // int percent = (int)(0.05*chainentries);

  std::cout << "Beginning of Energy spectra building research" << std::endl;
  std::cout << "On " << chainentries << " entries" << std::endl;

  // Creating a ProgressBar which is nice in MT mode
  using namespace indicators;
  ProgressBar bar{
    option::BarWidth{60},
    option::Start{"["},
    option::Fill{"="},
    option::Lead{">"},
    option::Remainder{"-"},
    option::End{"]"},
    option::PostfixText{"Reading"},
    option::ShowElapsedTime{true},
    option::ShowRemainingTime{true},
    option::ForegroundColor{Color::grey},
    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
  };

  // I read all the file of the TChain
  int fileindex  = 0;
  int filenumber =  experiment.GetDataFileNames().size();
  for(auto &file : experiment.GetDataFileNames())
  {
    fileindex++;
    //PrintVector(experiment.GetDataFileNames());
    bar.set_progress((int)((double)fileindex/(double)filenumber*100.));
    // Defining a TTreeProcessor object to handle the TTree in the File in a MT mode
    ROOT::TTreeProcessorMT TP(file, experiment.GetDataTreeName(), n_workers);

    // Scanning TTree
    // Launch the parallel processing of the tree
    auto loop_and_fill = [&] (TTreeReader &myReader){
      TTreeReaderValue<label_type> labelRV(myReader,"label");
      TTreeReaderValue<nrj_type> enrjRV(myReader,"nrj");


      while(myReader.Next())
      {
        label_type label = *labelRV;
        nrj_type NRJ   = *enrjRV;

        int spectrumindex = experiment.GetLabel2Detnbrs(static_cast<int>(label));
        //Double_t NRJ = enrj;//Faster2bitsNRJConverter(enrj, experiment,spectrumindex); // To convert to a sort of integer (it should remove the "non-linearity")

        // Calculating the right energy
        if(NRJ > 10) NRJspectra.at(spectrumindex)->Fill(NRJ);

        // And the energy Matrix
        if (NRJ > 10) NRJmatrix->Fill(label,NRJ);
      }
    };
    TP.Process(loop_and_fill);
  }


  std::cout << endl << RESETTEXT << "Saving in " << outputfilename << std::endl;
  outputfile->cd();

  for(int i = 0; i < (int) NRJspectra.size(); i++) NRJspectra.at(i)->Write();

  NRJmatrix->Write();
  outputfile->Close();

  cout << "NRJ spectra are saved " << endl;

  // Printing of chronometer measurement
  localTimer.Stop();
  Double_t rtime2 = localTimer.RealTime();
  Double_t ctime2 = localTimer.CpuTime();
  std::cout << std::endl;
  std::cout << "End of Energy Spectra Plotting" << std::endl;
  std::cout << "# RealTime=" << rtime2 << " seconds, CpuTime="<< ctime2 << " seconds" <<std::endl;
  std::cout << std::endl;

  return 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int DrawOneDetectorTypeEnergyCalibratedSpectra(const CExperiment &experiment, const TString DetectorType, Double_t Emin, Double_t Emax, Int_t binning)
{
  //Definition du germe pour le tirage aleatoire
  srand48(time(NULL));

  //Declaration of all variables for names
  TString title;
  TString spectrumname;

  ROOT::EnableThreadSafety();

  // Declaration of time spectra
  Int_t nbrofspectra = experiment.GetDetectors().size();
  std::vector<TH1F*> NRJspectra;
  Int_t highestdetlabel(0),highestnbrchannels(0),highestEmax(0);
  Int_t nbrchannels(0);
  
TH1F* localNRJspectrum;
// Defining the title and name of the spectrum
title = "Energy Spectrum of "; title+=DetectorType; title+= " detector ";
spectrumname = DetectorType; spectrumname += "cumulatednrjspectrum";
localNRJspectrum = new TH1F(spectrumname,title,binning,Emin,Emax);
  // I declare la Time matrix to check alignement later
  //TH2F* NRJmatrix = new TH2F("NRJalignementmatrix","NRJ spectra of all detectors",highestdetlabel,1,highestdetlabel,highestnbrchannels,Emin,highestEmax);

  // Creation of the file to save all the data
  TString outputfilename = experiment.GetFileDirectory_OUT();
  outputfilename += "CumulatedCalibratedEnergyspectra_all.root";
  TFile *outputfile = new TFile(outputfilename,"UPDATE");

  // Loading the TTree for reading the Data
  std::vector<TChain *> tab_chained_oak = experiment.GettheTChain();
  TChain *chained_oak = tab_chained_oak.at(0);

  label_Rawtype index;
  nrj_type enrj;
  std::cout << FOREBLU << "Loading Tree ..." << std::endl;
  /*auto cachesize = 10000000; // 10 MBytes
  chained_oak -> SetCacheSize(cachesize);
  chained_oak -> SetCacheLearnEntries(100000);*/
  chained_oak -> SetBranchStatus("*",0);
  chained_oak -> SetBranchAddress("label",&index);   // Detector number
  chained_oak -> SetBranchAddress("nrj",&enrj);      // NRJ of the hit
  std::cout << "Loaded .." << RESETTEXT << std::endl;

  // Another Timer
  TStopwatch localTimer;

  //Starting of the chronometer
  localTimer.Reset();
  localTimer.Start();

  //---------------------------------------------------------------------------//
  //                                                                           //
  //                            Reading the TTree                              //
  //                                                                           //
  //---------------------------------------------------------------------------//
  // // Determination of the number of event that have to be treated
  // const auto max_workers = std::thread::hardware_concurrency();
  auto n_workers = experiment.GetThreadsNbr();
  // if(n_workers > max_workers){
  //   cout << "Current machine only supports " << max_workers << " parallel threads. I will be using this number" << endl;
  //   n_workers = max_workers;
  // }
  const auto chainentries = chained_oak -> GetEntries();
  // const auto range = chainentries/n_workers;
  // int percent = (int)(0.05*chainentries);

  std::cout << "Beginning of Energy spectra building research" << std::endl;
  std::cout << "On " << chainentries << " entries" << std::endl;

  // Creating a ProgressBar which is nice in MT mode
  using namespace indicators;
  ProgressBar bar{
    option::BarWidth{60},
    option::Start{"["},
    option::Fill{"="},
    option::Lead{">"},
    option::Remainder{"-"},
    option::End{"]"},
    option::PostfixText{"Reading"},
    option::ShowElapsedTime{true},
    option::ShowRemainingTime{true},
    option::ForegroundColor{Color::grey},
    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
  };

  // I read all the file of the TChain
  int fileindex  = 0;
  int filenumber =  experiment.GetDataFileNames().size();
  for(auto &file : experiment.GetDataFileNames())
  {
    fileindex++;
    //PrintVector(experiment.GetDataFileNames());
    bar.set_progress((int)((double)fileindex/(double)filenumber*100.));
    // Defining a TTreeProcessor object to handle the TTree in the File in a MT mode
    ROOT::TTreeProcessorMT TP(file, experiment.GetDataTreeName(), n_workers);

    // Scanning TTree
    // Launch the parallel processing of the tree
    auto loop_and_fill = [&] (TTreeReader &myReader){
      TTreeReaderValue<label_type> labelRV(myReader,"label");
      TTreeReaderValue<nrj_type> enrjRV(myReader,"nrj");


      while(myReader.Next())
      {
        label_type label = *labelRV;
        nrj_type NRJ   = *enrjRV;

        int spectrumindex = experiment.GetLabel2Detnbrs(static_cast<int>(label));
        TString localdettype = experiment.GetDetector(spectrumindex)->GetDetectorType();
        //Double_t NRJ = enrj;//Faster2bitsNRJConverter(enrj, experiment,spectrumindex); // To convert to a sort of integer (it should remove the "non-linearity")
        //cout << "Label = " << label << " corresponding to a " << localdettype << endl;
        // Calculating the right energy
        if(NRJ > 10 && localdettype == DetectorType) localNRJspectrum->Fill(NRJ);
      }
    };
    TP.Process(loop_and_fill);
  }


  std::cout << endl << RESETTEXT << "Saving in " << outputfilename << std::endl;
  outputfile->cd();
  localNRJspectrum->Write();
  //for(int i = 0; i < (int) NRJspectra.size(); i++) NRJspectra.at(i)->Write();

  //NRJmatrix->Write();
  outputfile->Close();

  cout << "NRJ spectra are saved " << endl;

  // Printing of chronometer measurement
  localTimer.Stop();
  Double_t rtime2 = localTimer.RealTime();
  Double_t ctime2 = localTimer.CpuTime();
  std::cout << std::endl;
  std::cout << "End of Energy Spectra Plotting" << std::endl;
  std::cout << "# RealTime=" << rtime2 << " seconds, CpuTime="<< ctime2 << " seconds" <<std::endl;
  std::cout << std::endl;

  return 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>> SetCalibrationSource(TString sourcename, Bool_t Calib_BGO, Bool_t Calib_Ge, Bool_t Calib_LaBr)
{
  // I define the energies I need to calibrate
  std::vector<Double_t> Ref_nrj;
  std::vector<Double_t> Ref_nrj_err;
  std::vector<Double_t> Ref_nrj_intensity;

  // Now I complete these as a function of the source:
  std::map<TString, Bool_t> isSource;
  isSource.insert(pair<TString, Bool_t>("Eu", kFALSE));
  isSource.insert(pair<TString, Bool_t>("Co", kFALSE));
  isSource.insert(pair<TString, Bool_t>("Cs", kFALSE));
  isSource.insert(pair<TString, Bool_t>("Ba", kFALSE));
  isSource.insert(pair<TString, Bool_t>("Am", kFALSE));
  isSource.insert(pair<TString, Bool_t>("Th", kFALSE));
  isSource.insert(pair<TString, Bool_t>("LaBr", kFALSE));
  isSource.insert(pair<TString, Bool_t>("Perso", kFALSE));

  // Now I need to decrypt the source name
  // It is a mono-source
  Bool_t isMono = kFALSE;
  if(sourcename.Index("+") == -1)
  {
    isSource[sourcename] = kTRUE;
    isMono = kTRUE;
  }
  else
  {
    while(sourcename.Index("+") != -1)
    {
      // I find the first source name
      int it1 = sourcename.Index("+",1,sourcename.kExact);
      isSource[sourcename(0,it1)] = kTRUE;
      sourcename = sourcename(it1+1,sourcename.Length());
    }
    isSource[sourcename] = kTRUE;
  }

  if(isSource["Eu"])
  {
    cout << FORECYN << "Adding 152Eu energies for calibration" << RESETTEXT << endl;
    if(Calib_Ge)
    {
      // int nrj_nbr = 8;
      // Double_t nrj[8]={121.7817,244.6975,344.2785,778.904,867.378,964.079,1112.074,1408.006};
      // Double_t nrj_err[8]={0.0003,.0008,.0012,.0024,.003,.018,.003,.003};
      // Double_t intensite[8]={0.2858,0.07583,0.265,0.12942,0.04245,0.14605,0.13644,0.21005};

      int nrj_nbr = 7;
      Double_t nrj[7]={121.7817,244.6975,344.2785,778.904,964.079,1408.006};
      Double_t nrj_err[7]={0.0003,.0008,.0012,.0024,.018,.003};
      Double_t intensite[7]={0.2858,0.07583,0.265,0.12942,0.14605,0.21005};

      // I initialize the table of energies
      for(int i = 0; i < nrj_nbr; i++)
      {
        Ref_nrj.push_back(nrj[i]);
        Ref_nrj_err.push_back(nrj_err[i]);
        Ref_nrj_intensity.push_back(intensite[i]);
      }
    }
    if(Calib_LaBr)
    {
      // int nrj_nbr = 7;
      // Double_t nrj[7]={121.7817,244.6975,344.2785,778.904,867.378,964.079,1408.006};
      // Double_t nrj_err[7]={0.0003,.0008,.0012,.0024,.003,.018,.003};
      // Double_t intensite[7]={0.2858,0.07583,0.265,0.12942,0.04245,0.14605,0.21005};
      int nrj_nbr = 5;
      Double_t nrj[5]={121.7817,244.6975,344.2785,778.904,1408.006};
      Double_t nrj_err[5]={0.0003,.0008,.0012,.0024,.003};
      Double_t intensite[5]={0.2858,0.07583,0.265,0.12942,0.21005};

      // I initialize the table of energies
      for(int i = 0; i < nrj_nbr; i++)
      {
        Ref_nrj.push_back(nrj[i]);
        Ref_nrj_err.push_back(nrj_err[i]);
        Ref_nrj_intensity.push_back(intensite[i]);
      }
    }
  }

  if(isSource["Co"])
  {
    cout << FORECYN << "Adding 60Co energies for calibration" << RESETTEXT << endl;
    int nrj_nbr = 2;
    Double_t nrj[2]={1173.237,1332.501};
    Double_t nrj_err[2]={0.004,.005};
    Double_t intensite[2]={0.999736,0.999856};

    // I initialize the table of energies
    for(int i = 0; i < nrj_nbr; i++)
    {
      Ref_nrj.push_back(nrj[i]);
      Ref_nrj_err.push_back(nrj_err[i]);
      Ref_nrj_intensity.push_back(intensite[i]);
    }
  }

  if(isSource["Cs"])
  {
    cout << FORECYN << "Adding 137Cs energies for calibration" << RESETTEXT << endl;
    int nrj_nbr = 1;
    Double_t nrj[1]={661.657};
    Double_t nrj_err[1]={0.003};
    Double_t intensite[1]={0.851};

    // I initialize the table of energies
    for(int i = 0; i < nrj_nbr; i++)
    {
      Ref_nrj.push_back(nrj[i]);
      Ref_nrj_err.push_back(nrj_err[i]);
      Ref_nrj_intensity.push_back(intensite[i]);
    }
  }

  if(isSource["Ba"])
  {
    cout << FORECYN << "Adding 133Ba energies for calibration" << RESETTEXT << endl;
    int nrj_nbr = 3;
    Double_t nrj[3]={80.9971,302.853,356.017};
    Double_t nrj_err[3]={0.0014,0.001,0.002};
    Double_t intensite[3]={0.3406,0.1833,0.6205};

    // I initialize the table of energies
    for(int i = 0; i < nrj_nbr; i++)
    {
      Ref_nrj.push_back(nrj[i]);
      Ref_nrj_err.push_back(nrj_err[i]);
      Ref_nrj_intensity.push_back(intensite[i]);
    }
  }

  if(isSource["Am"])
  {
    cout << FORECYN << "Adding 241Am energies for calibration" << RESETTEXT << endl;
    int nrj_nbr = 1;
    Double_t nrj[1]={59.5412};
    Double_t nrj_err[1]={0.0002};
    Double_t intensite[1]={0.359};

    // I initialize the table of energies
    for(int i = 0; i < nrj_nbr; i++)
    {
      Ref_nrj.push_back(nrj[i]);
      Ref_nrj_err.push_back(nrj_err[i]);
      Ref_nrj_intensity.push_back(intensite[i]);
    }
  }

  if(isSource["AmBe"])
  {
    cout << "Adding 241Am energies for calibration" << endl;
    int nrj_nbr = 4;
    Double_t nrj[4]={59.5412,3416.91,3927.91,4438.91};
    Double_t nrj_err[4]={0.0002,1,1,1};
    Double_t intensite[4]={0.359,0,0,0.56};

    // I initialize the table of energies
    for(int i = 0; i < nrj_nbr; i++)
    {
      Ref_nrj.push_back(nrj[i]);
      Ref_nrj_err.push_back(nrj_err[i]);
      Ref_nrj_intensity.push_back(intensite[i]);
    }
  }

  if(isSource["Th"])
  {
    cout << FORECYN << "Adding natTh energies for calibration" << RESETTEXT << endl;
    int nrj_nbr = 10;
    Double_t nrj[10]={238.632,338.320,463.004,510.77,583.191,727.330,794.947,911.204,968.971,2614.533};
    Double_t nrj_err[10]={0.002,0.003,0.006,0.1,0.002,0.009,0.005,0.004,0.017,0.013};
    Double_t intensite[10]={0.4330,0.1127,0.0044,0.226,0.317,0.0658,0.258,0.158,0.3725};

    // I initialize the table of energies
    for(int i = 0; i < nrj_nbr; i++)
    {
      Ref_nrj.push_back(nrj[i]);
      Ref_nrj_err.push_back(nrj_err[i]);
      Ref_nrj_intensity.push_back(intensite[i]);
    }
  }

  if(isSource["Perso"])
  {
    TString filename = "BOUIK";
    cout << FORECYN << "Adding specific energies from file" << filename << RESETTEXT << endl;
    int nrj_nbr = 1;
    Double_t nrj[1]={666};
    Double_t nrj_err[1]={0.666};
    Double_t intensite[1]={0.666};

    // I initialize the table of energies
    for(int i = 0; i < nrj_nbr; i++)
    {
      Ref_nrj.push_back(nrj[i]);
      Ref_nrj_err.push_back(nrj_err[i]);
      Ref_nrj_intensity.push_back(intensite[i]);
    }
  }

  // If I calibrate LaBr3 I need to add some peaks associated to intrinsic activity
  if(isSource["LaBr"])
  {
    cout << FORECYN << "Adding LaBr3 intrinsic energies for calibration" << RESETTEXT << endl;
    int nrj_nbr = 2;
    Double_t nrj[10]={1435.795,1460.830};
    Double_t nrj_err[10]={0.01,0.003};
    Double_t intensite[10]={0.66,0.11};

    // I initialize the table of energies
    for(int i = 0; i < nrj_nbr; i++)
    {
      Ref_nrj.push_back(nrj[i]);
      Ref_nrj_err.push_back(nrj_err[i]);
      Ref_nrj_intensity.push_back(intensite[i]);
    }
  }

  // I order the energies
  if(!isMono)
  {
    //auto p = sort_permutation(Ref_nrj,[](T const& a, T const& b){ /*some comparison*/ });
    // auto p = sort_permutation(Ref_nrj,[](T const& a, T const& b));
    //
    // Ref_nrj = apply_permutation(Ref_nrj, p);
    // Ref_nrj_err = apply_permutation(Ref_nrj_err, p);
    // Ref_nrj_intensity = apply_permutation(Ref_nrj_intensity, p);
    sortVectors(Ref_nrj, less<Double_t>(), Ref_nrj, Ref_nrj_err, Ref_nrj_intensity);
  }

  return std::make_tuple(Ref_nrj, Ref_nrj_err, Ref_nrj_intensity);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::tuple<int, TGraphErrors*, TF1*, TF1*, TGraph*, Double_t, Double_t , Double_t ,Double_t> OneSpectrumPeakSearchandAnalysis(const CExperiment &experiment, Bool_t Calib_BGO, Bool_t Calib_Ge, Bool_t Calib_LaBr, TH1F *spectrum, int sindex, std::vector<Double_t> Ref_nrj, std::vector<Double_t> Ref_nrj_err, std::vector<Double_t> Ref_nrj_intensity)
{
  cout << RESETTEXT << FOREGRN << endl << "Calibrating Detector : " << experiment.GetDetectors().at(sindex)->GetDetectorName() << RESETTEXT << endl;
  // First I define a new spectrum
  TSpectrum *s = new TSpectrum(Ref_nrj.size()*4); // I do not allow to find more than four times the number of reference peaks

  // I try to determinate the background with compton edges
  Int_t nbins = spectrum->GetXaxis()->GetNbins();
  Double_t xmin  = 0;
  Double_t xmax  = spectrum->GetXaxis()->GetXmax();
  Double_t * source = new Double_t[nbins];
  //Double_t * dest = new Double_t[nbins];
  //Double_t sigma = 6;
  for (int i = 0; i < nbins; i++) source[i]=spectrum->GetBinContent(i + 1);

  // Now I search for the background
  //int niteration = 10;
  if(Calib_LaBr)s->Background(source,nbins,50,TSpectrum::kBackDecreasingWindow,
                              TSpectrum::kBackOrder4,kFALSE,
                              TSpectrum::kBackSmoothing15,kFALSE);
  if(Calib_Ge) s->Background(source,nbins,30,TSpectrum::kBackDecreasingWindow,
                             TSpectrum::kBackOrder2,kTRUE,
                             TSpectrum::kBackSmoothing15,kTRUE);

  Coucou(1);
  // I prepare a TH1 to store the background
  TString d1name; d1name.Form("d1_%i",sindex);
  TH1F *d1 = new TH1F(d1name,"",nbins,xmin,xmax);

  Coucou(2);
  // I fill the spectrum witht the background
  for (int i = 0; i < nbins; i++) d1->SetBinContent(i + 1,source[i]);
Coucou(3);
  // Now I clean the spectrum from the background
  spectrum->Add(spectrum,d1,1,-1);
  delete d1;
Coucou(4);
  // I search for the peaks
  Int_t nfound(0);
  if(experiment.GetDetector(sindex)->GetDetectorType() == "LaBr" && spectrum->GetEntries() > 1000)
  {
    Double_t intensity = 0.05;
    while(nfound < (int)Ref_nrj.size())
    {
      nfound = s->Search(spectrum,8,"",intensity);
      intensity = intensity - 0.01;
    }
  }
  if(experiment.GetDetector(sindex)->GetDetectorType() == "BGO" && spectrum->GetEntries() > 1000 ) nfound = s->Search(spectrum,40,"",0.01);
  if(experiment.GetDetector(sindex)->GetDetectorType() == "Ge"  && spectrum->GetEntries() > 1000 )
  {
    Double_t intensity = 0.05;
    while(nfound < (int)Ref_nrj.size() && intensity > 0.01)
    {
      nfound = s->Search(spectrum,4,"",intensity);//
      intensity = intensity - 0.01;
    }
  }

  if (nfound != 0 || nfound < (int)Ref_nrj.size()) {cout <<"Found " << nfound << " candidate peaks to fit" << endl;}
  else
  {
    cout << BACKRED << FOREWHT << "No peaks were found in the spectrum" << endl;
    if(spectrum->GetEntries() == 0) cout << "The spectrum is empty for detector " << experiment.GetDetectors().at(sindex)->GetDetectorName() << endl;
    else cout << "A problem occured with detector " << experiment.GetDetectors().at(sindex)->GetDetectorName() <<BACKBLK <<  RESETTEXT << endl;
  }

  // I loop on the peaks
  // I load from TSpectrum the peak list
  Double_t *xpeaks = s->GetPositionX();
  Double_t *ypeaks = s->GetPositionY();
  Int_t npeaks = 0;
  Double_t sigma=6.;

  // I declare the Tables I need to fit:
  std::vector<Double_t> weight;
  std::vector<Double_t> pos;
  std::vector<Double_t> pos_err;
  std::vector<Double_t> reso;
  std::vector<Double_t> reso_err;
  std::vector<Double_t> integral;
  std::vector<Double_t> integral_err;
  //std::vector<Double_t> FWTM;
  // I calculate the ratio of the first two energies from the source
  Double_t initial_NRJ_ratio = Ref_nrj.at(1)/Ref_nrj.at(0);
  TString f1name;
  // I'm gonna fit all the found peaks
  if(Calib_Ge)
  {
    std::vector<Double_t> xpeak;
    std::vector<Double_t> ypeak;

    for(auto p = 0; p<nfound;p++)
    {
      xpeak.push_back(xpeaks[p]);
      ypeak.push_back(ypeaks[p]);
    }
    // Double_t memory_intense_nrj;
    // if(xpeak.at(0) > 60) memory_intense_nrj =  xpeak.at(0);
    // else memory_intense_nrj =  xpeak.at(1);
    //
    sortVectors(xpeak, less<Double_t>(), xpeak,ypeak);
    PrintVector(xpeak);
    // Now they are sorted I need to identify the right peaks...
    // Now i search for this energy ratio in the peak that were found:
    int pos_1stpeak=0;int pos_2ndpeak = 0;
    Bool_t leaveloop = kFALSE;
    for(auto first = 0; first < (int)xpeak.size() && !leaveloop;first++)
    {
      for(auto second = first+1; second < (int)xpeak.size()&& !leaveloop;second++)
      {
        Double_t data_ratio = xpeak.at(second)/xpeak.at(first);
        //cout << xpeak.at(first) << "\t" << xpeak.at(second) << "\t" << data_ratio << " for " << initial_NRJ_ratio << " & memory_intense_nrj = " << memory_intense_nrj << endl;

        if(TMath::Abs(initial_NRJ_ratio-data_ratio) < 0.15 && xpeak.at(first) > 5000.)
        {
          pos_1stpeak = first;pos_2ndpeak=second;
          leaveloop = kTRUE;
          //cout << "Pair found in " << pos_1stpeak << " and " << pos_2ndpeak << endl;
        }
      }
    }

    // Now I've identified the first peaks, I get the realistic list of peaks
    // If I don't start from the right one I erase lower energy peaks
    if(pos_1stpeak !=0)
    {
      xpeak.erase(xpeak.begin(),xpeak.begin()+pos_1stpeak);
      ypeak.erase(ypeak.begin(),ypeak.begin()+pos_1stpeak);
    }
    printVector(xpeak);

    nfound=xpeak.size();

    for (auto p=0;p<nfound;p++)//nfound;p++)
    {
      // First I get bin center
      Double_t xp = xpeak.at(p);
      //cout << endl << endl << "One peak found @ " << xp << endl;
      std::vector<Double_t> temp_amplitude;
      std::vector<Double_t> temp_pos;
      std::vector<Double_t> temp_sigma;

      // Define my gaussian around my peak
      f1name.Form("f1_%i",sindex);
      TF1 *f1 = new TF1(f1name, "gaus", xp-(2*sigma), xp+(2*sigma));
      f1->SetParameter(1,xp);
      //f1->SetParLimits(1,xp-sigma,xp+sigma);

      // Then I get amplitude
      Int_t bin = spectrum->GetXaxis()->FindBin(xp);
      Double_t yp = spectrum->GetBinContent(bin);
      f1->SetParameter(0,yp);
      //f1->SetParLimits(0,yp-TMath::Sqrt(yp),yp+TMath::Sqrt(yp));

      // I fix the sigma
      f1->SetParameter(2,sigma);
      //f1->SetParLimits(2,0.,(Double_t)3*sigma);

      // Printing for debug
      //cout << xp-(3*sigma) << "\t" << xp+(3*sigma) << endl;
      //cout << xp << " p/m " << xp-sigma << "\t"<<xp+sigma<< endl;
      //cout << yp << " p/m " << yp-TMath::Sqrt(yp) << "\t"<<yp+TMath::Sqrt(yp)<< endl;
      //cout << sigma << " p/m " << 0 << "\t"<<3*sigma<< endl;

      // Then I fit my spectrum with a gaussian
      //cout << "Fitting with f1" << endl;
      spectrum->Fit(f1,"RIQE");

      // I store the first fit info
      temp_amplitude.push_back(f1->GetParameter(0));
      temp_pos.push_back(f1->GetParameter(1));
      weight.push_back(1);
      temp_sigma.push_back(f1->GetParameter(2));

      // Printing the result out
      // cout << "temp Peak amplitude " << temp_amplitude << endl;
      //cout << "temp Peak mean value " << f1->GetParameter(1) << endl;
      // cout << "temp Peak sigma " << temp_sigma << endl;

      // I get a second fit with a linear bckground
      f1name.Form("f2_%i",sindex);
      TF1 *f2 = new TF1(f1name, "gaus(0)+pol1(3)", xp-(2*sigma), xp+(2*sigma));
      // cout << xp-(3*sigma) << "\t" << xp+(3*sigma) << endl;
      // cout << xp << " p/m " << xp-sigma << "\t"<<xp+sigma<< endl;
      // cout << yp << " p/m " << yp-TMath::Sqrt(yp) << "\t"<<yp+TMath::Sqrt(yp)<< endl;
      // cout << sigma << " p/m " << 0 << "\t"<<3*sigma<< endl << endl << endl;
      f2->SetParameter(0,temp_amplitude.at(0));
      f2->SetParameter(1,temp_pos.at(0));
      f2->SetParameter(2,temp_sigma.at(0));
      // f2->SetParameter(0,yp);
      // f2->SetParameter(1,xp);
      // f2->SetParameter(2,sigma);
      // f2->SetParLimits(0,yp-TMath::Sqrt(yp),yp+TMath::Sqrt(yp));f2->SetParameter(0,temp_amplitude);
      // f2->SetParLimits(1,xp-sigma,xp+sigma);f2->SetParameter(1,temp_pos);
      // f2->SetParLimits(2,0,3*sigma);f2->SetParameter(2,temp_sigma);
      //cout << "Fitting with f2" << endl;
      spectrum->Fit(f2,"RIQE");
      temp_amplitude.push_back(f2->GetParameter(0));
      temp_pos.push_back(f2->GetParameter(1));
      weight.push_back(2);
      temp_sigma.push_back(f2->GetParameter(2));

      //cout << f2->GetParameter(1) << endl;

      // And now I Fit with a skewed gaussian
      // TF1 * f3 = new TF1("sgf","2.*gaus(x,[0],[1],[2])*ROOT::Math::normal_cdf([3]*x,1,0)",xp-(2*sigma), xp+(2*sigma));
      // f3->SetParameter(0,temp_amplitude.at(1));
      // f3->SetParameter(1,temp_pos.at(1));
      // f3->SetParameter(2,temp_sigma.at(1));
      // NRJSpectra.at(sindex)->Fit("sgf","RIQE");
      // temp_amplitude.push_back(f3->GetParameter(0));
      // temp_pos.push_back(f3->GetParameter(1));
      // temp_sigma.push_back(f3->GetParameter(2));

      //cout << "Fitting with f3" << endl;
      f1name.Form("MyFit_%i",sindex);
      TF1 * fFitFunction = new TF1(f1name,DoubleTailedStepedGaussian,xp-(2*sigma), xp+(2*sigma),10);
      fFitFunction->SetParName(0, "NumberOfPeaks");
      fFitFunction->SetParName(1, "BkgConst");
      fFitFunction->SetParName(2, "BkgSlope");
      fFitFunction->SetParName(3, "BkgExp");
      fFitFunction->SetParName(4+0, "Height");
      fFitFunction->SetParName(4+1, "Position");
      fFitFunction->SetParName(4+2, "FWHM");
      fFitFunction->SetParName(4+3, "LeftTail");
      fFitFunction->SetParName(4+4, "RightTail");
      fFitFunction->SetParName(4+5, "AmplitudeStep");
      fFitFunction->FixParameter(0, 1);
      fFitFunction->SetParameter(1,f2->GetParameter(3));
      fFitFunction->SetParameter(2,f2->GetParameter(4));
      fFitFunction->FixParameter(3, 0.);
      fFitFunction->SetParameter(4,f2->GetParameter(0));
      fFitFunction->SetParameter(5,f2->GetParameter(1));
      fFitFunction->SetParameter(6,f2->GetParameter(2));
      fFitFunction->SetParameter(7,-2.);
      fFitFunction->SetParLimits(7,-5.,-0.1);
      fFitFunction->SetParameter(8,2.);
      fFitFunction->SetParLimits(8,0.1,5);
      fFitFunction->SetParameter(9,0.01);
      fFitFunction->SetParLimits(9,-1.,1.);

      // cout << "Before fitting " << endl;
      // for(auto p= 0; p<10; p++)
      // {
      //   cout << fFitFunction->GetParName(p) << "\t" << fFitFunction->GetParameter(p) << endl;
      // }

      spectrum->Fit(fFitFunction,"R0Q");
      // cout << "After fitting " << endl;
      // for(auto p= 0; p<10; p++)
      // {
      //   cout << fFitFunction->GetParName(p) << "\t" << fFitFunction->GetParameter(p) << endl;
      // }
      temp_amplitude.push_back(fFitFunction->GetParameter(4));
      temp_pos.push_back(fFitFunction->GetParameter(5));
      weight.push_back(10);
      temp_sigma.push_back(fFitFunction->GetParameter(6));

      //cout << f3->GetParameter(1) << endl;

      // Now I calculate the final parameters
      Double_t final_amplitude = VectorMean(temp_amplitude);
      Double_t final_pos = VectorWeightedMean(temp_pos,weight);
      Double_t final_pos_err = VectorStdDev(temp_pos);
      Double_t final_sigma = VectorMean(temp_sigma);
      Double_t final_sigma_err = VectorStdDev(temp_sigma);

      // Printing the result out
      // cout << "Peak amplitude " << final_amplitude << endl;
      //cout << "Peak mean value " << final_pos << "p/m " << final_pos_err << endl;
      // cout << "Peak sigma " << final_sigma << "p/m " << final_sigma_err<< endl<< endl << endl;

      // I calculate the integral of the peak
      Double_t final_Int = final_amplitude* final_sigma*TMath::Sqrt(2*TMath::Pi());
      Double_t a = final_amplitude; Double_t c = final_sigma; Double_t Delta_a = TMath::Sqrt(final_amplitude);Double_t Delta_c = final_sigma_err;
      Double_t cons = 2*TMath::Pi();
      Double_t final_Int_syst_err = TMath::Sqrt(TMath::Power(c,2)*TMath::Power(Delta_a,2)*cons+TMath::Power(a,2)*TMath::Power(Delta_c,2)*cons);
      Double_t final_Int_err = TMath::Sqrt(TMath::Power(final_Int_syst_err,2)+final_Int);

      // I fill the sortVectors
      pos.push_back(final_pos);
      pos_err.push_back(final_pos_err);
      reso.push_back(final_sigma*2.35482);
      reso_err.push_back(final_sigma_err*2.35482);
      integral.push_back(final_Int);
      integral_err.push_back(final_Int_err);

      // I stored all the peaks... Now I have to assign them
      // cout << " and fitted @ " << final_pos << " with a difference of " << TMath::Abs(xp-final_pos) << endl;
      npeaks++;

    }
    d1->Delete();
  }

  if(Calib_LaBr)
  {
    std::vector<Double_t> xpeak;
    std::vector<Double_t> ypeak;

    for(auto p = 0; p<nfound;p++)
    {
      xpeak.push_back(xpeaks[p]);
      ypeak.push_back(ypeaks[p]);
    }
    Double_t memory_intense_nrj;
    if(xpeak.at(0) > 60) memory_intense_nrj =  xpeak.at(0);
    else memory_intense_nrj =  xpeak.at(1);
    //printVector(xpeak);
    sortVectors(xpeak, less<Double_t>(), xpeak,ypeak);
    //printVector(xpeak);

    // Now they are sorted I need to identify the right peaks...
    // Now i search for this energy ratio in the peak that were found:
    int pos_1stpeak=0;int pos_2ndpeak=0;
    Bool_t leaveloop = kFALSE;
    for(auto first = 0; first < (int)xpeak.size() && !leaveloop;first++)
    {
      for(auto second = first+1; second < (int)xpeak.size()&& !leaveloop;second++)
      {
        Double_t data_ratio = xpeak.at(second)/xpeak.at(first);
        //cout << xpeak.at(first) << "\t" << xpeak.at(second) << "\t" << data_ratio << " for " << initial_NRJ_ratio << " & memory_intense_nrj = " << memory_intense_nrj << endl;

        if(TMath::Abs(initial_NRJ_ratio-data_ratio) < 0.15 && xpeak.at(first) == memory_intense_nrj)
        {
          pos_1stpeak = first;pos_2ndpeak=second;
          leaveloop = kTRUE;
          //cout << "Pair found in " << pos_1stpeak << " and " << pos_2ndpeak << endl;
        }
      }
    }

    // Now I've identified the first peaks, I get the realistic list of peaks
    // If I don't start from the right one I erase lower energy peaks
    if(pos_1stpeak !=0)
    {
      xpeak.erase(xpeak.begin(),xpeak.begin()+pos_1stpeak);
      ypeak.erase(ypeak.begin(),ypeak.begin()+pos_1stpeak);
    }
    //printVector(xpeak);

    nfound=xpeak.size();
    for (auto p=0;p<nfound;p++)//nfound;p++)
    {
      if(xpeak.at(p)/xpeak.at(0) < 9.91)
      {
        // First I get bin center
        Double_t xp = xpeak.at(p);
        //cout << "One peak found @ " << xp << endl;

        // Define my gaussian around my peak
        TF1 *f1 = new TF1("f1", "gaus", xp-(2*sigma), xp+(2*sigma));
        f1->SetParameter(1,xp);
        //f1->SetParLimits(1,xp-sigma,xp+sigma);

        // Then I get amplitude
        Int_t bin = spectrum->GetXaxis()->FindBin(xp);
        Double_t yp = spectrum->GetBinContent(bin);
        f1->SetParameter(0,yp);
        //f1->SetParLimits(0,yp-TMath::Sqrt(yp),yp+TMath::Sqrt(yp));

        // I fix the sigma
        f1->SetParameter(2,sigma);
        //f1->SetParLimits(2,0.,(Double_t)3*sigma);

        // Printing for debug
        //cout << xp-(3*sigma) << "\t" << xp+(3*sigma) << endl;
        //cout << xp << " p/m " << xp-sigma << "\t"<<xp+sigma<< endl;
        //cout << yp << " p/m " << yp-TMath::Sqrt(yp) << "\t"<<yp+TMath::Sqrt(yp)<< endl;
        //cout << sigma << " p/m " << 0 << "\t"<<3*sigma<< endl;

        // Then I fit my spectrum with a gaussian
        spectrum->Fit("f1","RIQE");

        // I store the first fit info
        Double_t temp_amplitude = f1->GetParameter(0);
        Double_t temp_pos = f1->GetParameter(1);
        Double_t temp_sigma = f1->GetParameter(2);

        // Printing the result out
        // cout << "temp Peak amplitude " << temp_amplitude << endl;
        // cout << "temp Peak mean value " << temp_pos << endl;
        // cout << "temp Peak sigma " << temp_sigma << endl;

        // I get a second fit with a linear bckground
        TF1 *f2 = new TF1("f2", "gaus(0)+pol1(3)", xp-(3*sigma), xp+(3*sigma));
        // cout << xp-(3*sigma) << "\t" << xp+(3*sigma) << endl;
        // cout << xp << " p/m " << xp-sigma << "\t"<<xp+sigma<< endl;
        // cout << yp << " p/m " << yp-TMath::Sqrt(yp) << "\t"<<yp+TMath::Sqrt(yp)<< endl;
        // cout << sigma << " p/m " << 0 << "\t"<<3*sigma<< endl << endl << endl;
        f2->SetParameter(0,temp_amplitude);
        f2->SetParameter(1,temp_pos);
        f2->SetParameter(2,temp_sigma);
        // f2->SetParameter(0,yp);
        // f2->SetParameter(1,xp);
        // f2->SetParameter(2,sigma);
        // f2->SetParLimits(0,yp-TMath::Sqrt(yp),yp+TMath::Sqrt(yp));f2->SetParameter(0,temp_amplitude);
        // f2->SetParLimits(1,xp-sigma,xp+sigma);f2->SetParameter(1,temp_pos);
        // f2->SetParLimits(2,0,3*sigma);f2->SetParameter(2,temp_sigma);
        spectrum->Fit("f2","RIQE");

        // Now I calculate the final parameters
        Double_t final_amplitude = (f2->GetParameter(0)+temp_amplitude)/2.;
        Double_t final_pos = (f2->GetParameter(1)+temp_pos)/2.;
        Double_t final_pos_err = TMath::Abs(final_pos-temp_pos);
        Double_t final_sigma = (f2->GetParameter(2)+temp_sigma)/2.;
        Double_t final_sigma_err = TMath::Abs(final_sigma-temp_sigma);

        // Printing the result out
        // cout << "Peak amplitude " << final_amplitude << endl;
        // cout << "Peak mean value " << final_pos << "p/m " << final_pos_err << endl;
        // cout << "Peak sigma " << final_sigma << "p/m " << final_sigma_err<< endl<< endl << endl;

        // I calculate the integral of the peak
        Double_t final_Int = final_amplitude* final_sigma*TMath::Sqrt(2*TMath::Pi());
        Double_t a = final_amplitude; Double_t c = final_sigma; Double_t Delta_a = TMath::Sqrt(final_amplitude);Double_t Delta_c = final_sigma_err;
        Double_t cons = 2*TMath::Pi();
        Double_t final_Int_syst_err = TMath::Sqrt(TMath::Power(c,2)*TMath::Power(Delta_a,2)*cons+TMath::Power(a,2)*TMath::Power(Delta_c,2)*cons);
        Double_t final_Int_err = TMath::Sqrt(TMath::Power(final_Int_syst_err,2)+final_Int);

        // I fill the sortVectors
        pos.push_back(final_pos);
        pos_err.push_back(final_pos_err);
        reso.push_back(final_sigma*2.35482);
        reso_err.push_back(final_sigma_err*2.35482);
        integral.push_back(final_Int);
        integral_err.push_back(final_Int_err);

        // I stored all the peaks... Now I have to assign them
        //cout << " and fitted @ " << final_pos << endl;
        npeaks++;
      }
      else // Normally I adress the question of the 1408 keV plus intrinsic activity in case resolution is shitty
      {
        // First I get bin center
        Double_t xp = xpeak.at(p);
        //cout << "One peak found @ " << xp << endl;

        // Define my gaussian around my peak
        TF1 *f1 = new TF1("f1", "gaus", xp-(2*sigma), xp+(2*sigma));
        f1->SetParameter(1,xp);
        //f1->SetParLimits(1,xp-sigma,xp+sigma);

        // Then I get amplitude
        Int_t bin = spectrum->GetXaxis()->FindBin(xp);
        Double_t yp = spectrum->GetBinContent(bin);
        f1->SetParameter(0,yp);
        //f1->SetParLimits(0,yp-TMath::Sqrt(yp),yp+TMath::Sqrt(yp));

        // I fix the sigma
        f1->SetParameter(2,sigma);
        //f1->SetParLimits(2,0.,(Double_t)3*sigma);

        // Printing for debug
        // cout << xp-(3*sigma) << "\t" << xp+(3*sigma) << endl;
        // cout << xp << " p/m " << xp-sigma << "\t"<<xp+sigma<< endl;
        // cout << yp << " p/m " << yp-TMath::Sqrt(yp) << "\t"<<yp+TMath::Sqrt(yp)<< endl;
        // cout << sigma << " p/m " << 0 << "\t"<<3*sigma<< endl;

        // Then I fit my spectrum with a gaussian
        spectrum->Fit("f1","RIQE");

        // I store the first fit info
        Double_t temp_amplitude = f1->GetParameter(0);
        Double_t temp_pos = f1->GetParameter(1);
        Double_t temp_sigma = f1->GetParameter(2);

        // Printing the result out
        // cout << "temp Peak amplitude " << temp_amplitude << endl;
        // cout << "temp Peak mean value " << temp_pos << endl;
        // cout << "temp Peak sigma " << temp_sigma << endl;

        //I get a second fit with a linear bckground
        TF1 *f2 = new TF1("f2", "gaus(0)+pol1(3)", xp-(2*sigma), xp+(2*sigma));
        // cout << xp-(3*sigma) << "\t" << xp+(3*sigma) << endl;
        // cout << xp << " p/m " << xp-sigma << "\t"<<xp+sigma<< endl;
        // cout << yp << " p/m " << yp-TMath::Sqrt(yp) << "\t"<<yp+TMath::Sqrt(yp)<< endl;
        // cout << sigma << " p/m " << 0 << "\t"<<3*sigma<< endl << endl << endl;
        f2->SetParameter(0,temp_amplitude);
        f2->SetParameter(1,temp_pos);
        f2->SetParameter(2,temp_sigma);
        // f2->SetParameter(0,yp);
        // f2->SetParameter(1,xp);
        // f2->SetParameter(2,sigma);
        // f2->SetParLimits(0,yp-TMath::Sqrt(yp),yp+TMath::Sqrt(yp));f2->SetParameter(0,temp_amplitude);
        // f2->SetParLimits(1,xp-sigma,xp+sigma);f2->SetParameter(1,temp_pos);
        // f2->SetParLimits(2,0,3*sigma);f2->SetParameter(2,temp_sigma);
        spectrum->Fit("f2","RIQE");

        // Now I calculate the final parameters
        Double_t final_amplitude = (f2->GetParameter(0)+temp_amplitude)/2.;
        Double_t final_pos = (f2->GetParameter(1)+temp_pos)/2.;
        Double_t final_pos_err = TMath::Abs(final_pos-temp_pos);
        Double_t final_sigma = (f2->GetParameter(2)+temp_sigma)/2.;
        Double_t final_sigma_err = TMath::Abs(final_sigma-temp_sigma);

        // Printing the result out
        // cout << "Peak amplitude " << final_amplitude << endl;
        // cout << "Peak mean value " << final_pos << "p/m " << final_pos_err << endl;
        // cout << "Peak sigma " << final_sigma << "p/m " << final_sigma_err<< endl<< endl << endl;

        // I calculate the integral of the peak
        Double_t final_Int = final_amplitude* final_sigma*TMath::Sqrt(2*TMath::Pi());
        Double_t a = final_amplitude; Double_t c = final_sigma; Double_t Delta_a = TMath::Sqrt(final_amplitude);Double_t Delta_c = final_sigma_err;
        Double_t cons = 2*TMath::Pi();
        Double_t final_Int_syst_err = TMath::Sqrt(TMath::Power(c,2)*TMath::Power(Delta_a,2)*cons+TMath::Power(a,2)*TMath::Power(Delta_c,2)*cons);
        Double_t final_Int_err = TMath::Sqrt(TMath::Power(final_Int_syst_err,2)+final_Int);

        // I fill the sortVectors
        pos.push_back(final_pos);
        pos_err.push_back(final_pos_err);
        reso.push_back(final_sigma*2.35482);
        reso_err.push_back(final_sigma_err*2.35482);
        integral.push_back(final_Int);
        integral_err.push_back(final_Int_err);

        // I stored all the peaks... Now I have to assign them
        //cout << " and fitted @ " << final_pos << endl;
        npeaks++;
      }
    }
    d1->Delete();
  }

  // Now I search the right peaks
  //printVector(pos);
  std::vector<int> Toremove;
  std::vector<Double_t> NRJratios;
  for(auto index = 1; index < (int)Ref_nrj.size(); index++)
  {
    //cout << index << "\t" << Ref_nrj.at(index) << "\t" <<Ref_nrj.at(0) << "\t" << Ref_nrj.at(index)/Ref_nrj.at(0) << endl;
    NRJratios.push_back(Ref_nrj.at(index)/Ref_nrj.at(0));
  }
  int locali = 0;
  for(auto index = 1; index  <  (int)pos.size(); index++)
  {
    Double_t data_ratio = pos.at(index)/pos.at(0);
    //Double_t integral_ratio = integral.at(index)/pos.at(0);
    //if(locali < (int)NRJratios.size())cout << "NRJratios.at("<<locali<<") = " << NRJratios.at(locali) << "\t with " << pos.at(index) << " & " << pos.at(0) << " giving " << data_ratio << " & int ratio = " << integral_ratio << endl;
    //
    // //Because of non-linearity I have to consider a second order term to correct the threshold
    // Double_t Ratio_corrected;
    // if(index > 2 && data_ratio > 9.5)
    // {
    //   Double_t data_ratio0 = pos.at(index-2)/pos.at(0);
    //   Double_t data_ratio_i_1 = pos.at(index-1)/pos.at(0);
    //   Double_t corr_slope = (data_ratio0-data_ratio_i_1)/(pos.at(0)-pos.at(index-1));
    //   Double_t corr_origin = data_ratio0 - corr_slope*pos.at(0);
    //   cout << corr_slope << "\t" << corr_origin << endl;
    //   Ratio_corrected = corr_slope*pos.at(index);
    //   cout << Ratio_corrected << endl;
    //   //data_ratio+=Ratio_corrected;
    // }
    // if(locali < (int)NRJratios.size())cout << "NRJratios.at("<<locali<<") = " << NRJratios.at(locali) << "\t with " << pos.at(index) << " & " << pos.at(0) << " giving after correction: " << data_ratio << endl << endl;
    Double_t precision = 0.26;
    if(Calib_LaBr) precision = 0.35;
    if(locali < (int)NRJratios.size())
    {
      if(TMath::Abs(NRJratios.at(locali)-data_ratio) < precision) locali++;
      else if(Calib_LaBr && NRJratios.at(locali) > 6. && TMath::Abs(1-data_ratio/NRJratios.at(locali)) < 0.12)locali++; // To consider non linearity that builds up at hiher energy
      else Toremove.push_back(index);
    }
    else
    {
      Toremove.push_back(index);
      locali++;
    }
  }

  // Now I remove the useless peaks
  for (int rmv =Toremove.size()-1; rmv > -1; rmv--)
  {
    pos.erase(pos.begin()+Toremove.at(rmv));
    pos_err.erase(pos_err.begin()+Toremove.at(rmv));
    reso.erase(reso.begin()+Toremove.at(rmv));
    reso_err.erase(reso_err.begin()+Toremove.at(rmv));
    integral.erase(integral.begin()+Toremove.at(rmv));
    integral_err.erase(integral_err.begin()+Toremove.at(rmv));
  }

  // I define the TGraph to plot calibration, resolution & efficiency
  printVector(pos);
  // printVector(Ref_nrj);
  // printVector(pos_err);
  // printVector(Ref_nrj_err);
  auto calibration = new TGraphErrors(pos.size(),&(pos[0]),&(Ref_nrj[0]),&(pos_err[0]),&(Ref_nrj_err[0]));
  TString title = "Calibration of "; title += experiment.GetDetectors().at(sindex)->GetDetectorName();
  TString title2 = "Residues for "; title2 += experiment.GetDetectors().at(sindex)->GetDetectorName();
  calibration->SetTitle(title);
  calibration->SetName(title);
  // auto resolution = new TGraphErrors((int)pos.size(),Ref_nrj,reso,Ref_nrj_err,reso_err);
  // auto efficiency = new TGraphErrors((int)pos.size(),Ref_nrj,integral,Ref_nrj_err,integral_err);
  TGraph *residue;

  // Now I fit the calibration
  f1name.Form("f1_1_%i",sindex);
  TF1 *f1_1 = new TF1(f1name,"pol1");
  f1name.Form("f2_1_%i",sindex);
  TF1 *f2_1 = new TF1(f1name,"pol2");
  f1name.Form("f3_1_%i",sindex);
  TF1 *f3_1 = new TF1(f1name, "pol3");
  TString fitname = "f_"; fitname += experiment.GetDetector(sindex)->GetDetectorName();
  TF1 *f3 = new TF1(fitname, "pol3");
  TF1 *f4 = new TF1(fitname,"pol1");
  cout << FOREGRN <<  "Peaks identified now fitting calibration curve " << RESETTEXT << endl;
  std::vector<Double_t> residues;
  Double_t a = 0;
  Double_t b = 0;
  Double_t c = 0;
  Double_t d = 0;
  if(Calib_Ge)
  {
    calibration->Fit(f4,"WQ");
    std::cout.precision(7);
    a = 0;
    b = 0;
    c = f4->GetParameter(1);
    d = f4->GetParameter(0);
    //calibfile << experiment.GetDetector(sindex)->GetDetectorlabel() << "\t" << a << "\t" << b << "\t" << c << "\t" << d << endl;

    // Now I calculate the Residues
    cout << FOREGRN << "Now Calculating Residues" << RESETTEXT << endl;
    for(auto nrj = 0 ; nrj < (int)pos.size(); nrj++)
    {
      residues.push_back(Ref_nrj.at(nrj)-(c*pos.at(nrj)+d));
    }
  }
  residue = new TGraph(residues.size(),&(Ref_nrj[0]),&(residues[0]));
  residue -> SetTitle(title2);
  residue -> SetName(title2);

  if(Calib_LaBr)
  {
    calibration->Fit(f1_1,"Q");
    for(auto param = 0; param < 2; param++)f2_1->FixParameter(param,f1_1->GetParameter(param));
    calibration->Fit(f2_1,"Q");
    for(auto param = 0; param < 3; param++)f2_1->SetParameter(param,f2_1->GetParameter(param));
    calibration->Fit(f2_1,"Q");
    for(auto param = 0; param < 3; param++)f3_1->FixParameter(param,f2_1->GetParameter(param));
    calibration->Fit(f3_1,"Q");
    for(auto param = 0; param < 4; param++)f3->SetParameter(param,f3_1->GetParameter(param));
    calibration->Fit(f3,"Q");
    // cout << "The linear fit result: " << endl;

    // cout << "Energy (kev) = " << f3->GetParameter(1) << "* x + " << f3->GetParameter(0) << endl;
    // cout << "With a probability " << f3->GetProb() << endl;
    a = f3->GetParameter(0);
    b = f3->GetParameter(1);
    c = f3->GetParameter(2);
    d = f3->GetParameter(3);
    //calibfile << experiment.GetDetector(sindex)->GetDetectorlabel() << "\t" << f3->GetParameter(3) << "\t" << f3->GetParameter(2) << "\t" << f3->GetParameter(1) << "\t" << f3->GetParameter(0) << endl;
  }



  return std::make_tuple(1,calibration,f4,f3,residue,a,b,c,d);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int DrawAllParisUncalibratedSpectra(const CExperiment &experiment)
{
  //Definition du germe pour le tirage aleatoire
  srand48(time(NULL));

  //Declaration of all variables for names
  TString title,title2,title3,title4;
  TString spectrumname,spectrumname2,spectrumname3,spectrumname4;

  ROOT::EnableThreadSafety();

  // Declaration of energy spectra
  // First I count the number of PARIS
  int nbrparis=0;
  for(int d=0; d < (int)experiment.GetDetectors().size();d++)
  {
    if(experiment.GetDetectors().at(d)->GetDetectorType() == "PARIS") nbrparis++;
  }

  Int_t nbrofspectra = (int)experiment.GetDetectors().size();
  std::vector<TH1F*> PSDSpectra;
  std::vector<TH2F*> PSDMatrix;
  std::vector<TH1F*> QDC1Spectra;
  std::vector<TH1F*> QDC2Spectra;
  Int_t highestdetlabel = 0;
  Int_t highestnbrchannels = 0;
  Int_t highestEmax = 0;
  Int_t nbrchannels = 0;
  Int_t Emin = 0;
  Int_t Emax = 0;

  // Defining all the energy spectra
  for(auto sindex = 0; sindex < nbrofspectra; sindex++)
  {
    //TH1F *localPSDSpectra;
    TH1F *localQDC1spectrum = nullptr;TH1F *localQDC2spectrum = nullptr; TH1F *localPSDSpectra = nullptr;
    TH2F *localPSDMatrix = nullptr;
    // Defining the title and name of the spectrum
    title = "LaBr Spectrum of detector ";title2 ="NaI Spectrum of detector ";title3 = " PSD Spectrum of detector ";title4 = " PSD Matrix of detector ";
    spectrumname = "nrjspectrum";spectrumname2 = "nrjspectrum2";spectrumname3 = "psdspectrum";spectrumname4 = "psdmatrix";
    title += experiment.GetDetector(sindex)->GetDetectorName();title2 += experiment.GetDetector(sindex)->GetDetectorName();title3 += experiment.GetDetector(sindex)->GetDetectorName();title4 += experiment.GetDetector(sindex)->GetDetectorName();
    spectrumname +=experiment.GetDetector(sindex)->GetDetectorName();spectrumname2 +=experiment.GetDetector(sindex)->GetDetectorName();spectrumname3 +=experiment.GetDetector(sindex)->GetDetectorName();spectrumname4 +=experiment.GetDetector(sindex)->GetDetectorName();

    if(experiment.GetDetector(sindex)->GetDetectorlabel()>highestdetlabel) highestdetlabel = experiment.GetDetectors().at(sindex)->GetDetectorlabel();

    // Getting the right range
    if(experiment.GetDetector(sindex)->GetDetectorType()!="RF")
    {
      nbrchannels = experiment.GetDetector(sindex)->GetNbrChannels();
      //cout << nbrchannels << "\t" << Emax << endl;
      if(nbrchannels > highestnbrchannels) highestnbrchannels=nbrchannels;
      Emax = experiment.GetDetector(sindex)->GetMaxchNumber();//GetMaxchNumber();
      if(Emax > highestEmax) highestEmax=Emax;
      //Emax = 200000;
      //nbrchannels = 200000;
      //if(experiment.GetDetector(sindex)->GetDetectorType()=="PARIS")
      {
        localQDC1spectrum = new TH1F(spectrumname,title,nbrchannels,Emin,Emax);
        localQDC2spectrum = new TH1F(spectrumname2,title2,nbrchannels,Emin,Emax);
        localPSDSpectra   = new TH1F(spectrumname3,title3,1600,0,TMath::Pi()/2.);
        localPSDMatrix    = new TH2F(spectrumname4,title4,nbrchannels,Emin,Emax,nbrchannels,Emin,Emax);
      }

      
      // Storing the NRJsectrum
      QDC1Spectra.push_back(localQDC1spectrum);
      QDC2Spectra.push_back(localQDC2spectrum);
      PSDSpectra.push_back(localPSDSpectra);
      PSDMatrix.push_back(localPSDMatrix);
      
    }

    // Cleaning the names for the next iteration
    spectrumname.Clear();spectrumname2.Clear();spectrumname3.Clear();spectrumname4.Clear();
    title.Clear();title2.Clear();title3.Clear();title4.Clear();
  }


  std::cout << FOREBLU << "We have " << nbrofspectra << " detectors among which " << nbrparis  << " PARIS phoswitches" << std::endl;
  std::cout << "We have created " << QDC1Spectra.size() << " QDC1 spectra among which " << nbrparis  << " PARIS phoswitches QDC1 spectra"  << std::endl;
  std::cout << "We have created " << QDC2Spectra.size() << " QDC2 spectra among which " << nbrparis  << " PARIS phoswitches QDC2 spectra"  << std::endl;
  std::cout << "We have created " << PSDSpectra.size() << " PSD spectra among which " << nbrparis  << " PARIS phoswitches PSD spectra"   << std::endl;
  std::cout << "We have created " << PSDMatrix.size() << " PSD Matrices among which " << nbrparis  << " PARIS phoswitches PSD Matrices"  << std::endl;
  std::cout << "We are Loaded .." << RESETTEXT << std::endl;

  // Creation of the file to save all the data
  TString outputfilename = experiment.GetFileDirectory_OUT();
  outputfilename += "UncalibratedPARISspectra_all.root";
  TFile *outputfile = new TFile(outputfilename,"RECREATE");

  // Loading the TTree for reading the Data
  std::vector<TChain *> tab_chained_oak = experiment.GettheTChain();
  TChain *chained_oak = tab_chained_oak.at(0);

  label_Rawtype index;
  nrj_Rawtype enrj;

   std::cout << FOREBLU << "Loading Tree ..." << std::endl;
  // /*auto cachesize = 10000000; // 10 MBytes
  // chained_oak -> SetCacheSize(cachesize);
  // chained_oak -> SetCacheLearnEntries(100000);*/
  // chained_oak -> SetBranchStatus("*",0);
  // chained_oak -> SetBranchAddress("label",&index);   // Detector number
  // chained_oak -> SetBranchAddress("nrj",&enrj);      // NRJ of the hit
  // std::cout << "Loaded .." << RESETTEXT << std::endl;

  // Another Timer
  TStopwatch localTimer;

  //Starting of the chronometer
  localTimer.Reset();
  localTimer.Start();

  //---------------------------------------------------------------------------//
  //                                                                           //
  //                            Reading the TTree                              //
  //                                                                           //
  //---------------------------------------------------------------------------//
  // // Determination of the number of event that have to be treated
  // const auto max_workers = std::thread::hardware_concurrency();
  auto n_workers = experiment.GetThreadsNbr();
  // if(n_workers > max_workers){
  //   cout << "Current machine only supports " << max_workers << " parallel threads. I will be using this number" << endl;
  //   n_workers = max_workers;
  // }
  const auto chainentries = chained_oak -> GetEntries();
  // const auto range = chainentries/n_workers;
  // int percent = (int)(0.05*chainentries);

  std::cout << "Beginning of Energy spectra building research" << std::endl;
  std::cout << "On " << chainentries << " entries" << std::endl;

  // Creating a ProgressBar which is nice in MT mode
  using namespace indicators;
  ProgressBar bar{
    option::BarWidth{60},
    option::Start{"["},
    option::Fill{"="},
    option::Lead{">"},
    option::Remainder{"-"},
    option::End{"]"},
    option::PostfixText{"Reading"},
    option::ShowElapsedTime{true},
    option::ShowRemainingTime{true},
    option::ForegroundColor{Color::grey},
    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
  };
  cout<<"All goodd"<<endl;

  // I read all the file of the TChain
  int fileindex  = 0;
  int filenumber =  experiment.GetDataFileNames().size();
  mutex forfilling;
  for(auto &file : experiment.GetDataFileNames())
  {
    fileindex++;
    //PrintVector(experiment.GetDataFileNames());
    bar.set_progress((int)((double)fileindex/(double)filenumber*100.));
    // Defining a TTreeProcessor object to handle the TTree in the File in a MT mode
    ROOT::TTreeProcessorMT TP(file, experiment.GetDataTreeName(), n_workers);

    // Scanning TTree
    // Launch the parallel processing of the tree
    int threadnbr = 0;
    auto loop_and_fill = [&] (TTreeReader &myReader){
      TTreeReaderValue<label_Rawtype> labelRV(myReader,"label");
      TTreeReaderValue<nrj_Rawtype> QDC1RV(myReader,"nrj");
      TTreeReaderValue<nrj_Rawtype> QDC2RV(myReader,"nrj2");

      ULong64_t hitnumber = 0;
      while(myReader.Next())
      {
        auto label  = *labelRV;
        auto NRJ    = *QDC1RV;  // Short Gate
        auto NRJ2   = *QDC2RV;  // Long Gate

        // I define a new hit a fill in the information
        hitnumber++;threadnbr++;
        if(label > 0) // To make sure I only consider PARIS detectors
        {
          CHit *hit = new CHit(hitnumber+1000*threadnbr);
          hit->SetHit(label, 0, NRJ, NRJ2, false);

          int spectrumindex = experiment.GetLabel2Detnbrs(static_cast<int>(label));
          
          Double_t PSD = hit->PerformPARISPSD();// PSD = atan(short/long)
          Bool_t isLaBr = experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->IsPureLaBr3(PSD,hit->GetHitE1(),hit->GetHitE2()); // I need the Long charge to get proper LaBr3 selection
          Bool_t isBeyondLaBr3andNaI = experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->IsBeyondLaBr3andNaI(PSD,hit->GetHitE1(),hit->GetHitE2()); // I need the Long charge to get proper LaBr3 selection

          //cout << " Detector " << (string)experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->GetDetectorName() << " discriLaBr_pos = " << experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->GetLaBrDiscriPosition() << endl;
          // Filling up the spectra
          //forfilling.lock();
          double NRJ_bf_ROT = hit->GetHitE1();
          if(experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->GetDetectorType() == "PARIS")
           {
             tie(NRJ2,NRJ) = Rotation(
                                      static_cast<double>(experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->GetRotationAngle()),
                                      static_cast<double>(experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->GetRotationAngleTan()),
                                      NRJ2,
                                      NRJ);
             //tie(NRJ2,NRJ) = Rotation(0.,0.,NRJ2,NRJ);
           }
          if(spectrumindex < PSDSpectra.size())
          {
            // cout << endl;
            // hit->PrintHit();
            // cout << "Hit # " << hitnumber+1000*threadnbr << "; Filling Spectra " << spectrumindex << endl;
            PSDSpectra.at(spectrumindex)  -> Fill(PSD);
            if(isLaBr) QDC1Spectra.at(spectrumindex) -> Fill(NRJ_bf_ROT);
            if(!isLaBr && ! isBeyondLaBr3andNaI) QDC2Spectra.at(spectrumindex) -> Fill(NRJ2);
          }
          
          if(spectrumindex < PSDMatrix.size())
          {
            //cout << "Hit # " << hitnumber+1000*threadnbr << "; Filling Matrix " << spectrumindex << endl;
            //if(isLaBr)PSDMatrix.at(spectrumindex)->Fill((Double_t)NRJ2,(Double_t)NRJ);
            if(! isBeyondLaBr3andNaI) PSDMatrix.at(spectrumindex)->Fill((Double_t)NRJ2,(Double_t)NRJ);
            //PSDMatrix.at(spectrumindex)->Fill((Double_t)NRJ2,(Double_t)NRJ);
          }
          //forfilling.unlock();

          delete hit;
        }
      }
    };
    TP.Process(loop_and_fill);
  }

  // All spectra has been built
  // Now analyzing the PSDSpectra to print to text files the parameters of PSD
  std::cout << endl << RESETTEXT << "Saving in " << outputfilename << std::endl;
  outputfile->cd();

  // Now I have the peak position for LaBr3 selection
  TString PSDoutputfilename = "newPSDParameter_PARIS.txt";
  ofstream PSDoutput(PSDoutputfilename, ios::out);
  if (!PSDoutput.is_open()) {
      std::cerr << "Error: Could not create or open newPSDParameter_PARIS.txt!" << std::endl;
  } 
  else {
      std::cout << "PSDPArameter_PARIS.txt has been created" << std::endl;
      PSDoutput << "Det Name \t LaBrPos \t LaBrSigma \t NaIPos \t NaISigma \t theta" << std::endl;
  }
  for(int i = 0; i < nbrofspectra; i++)
  {
    if(experiment.GetDetectors().at(i)->GetDetectorType()=="PARIS" && PSDSpectra.at(i)->GetEntries() !=0)
    {
      cout << "Saving Detector " << experiment.GetDetectors().at(i)->GetDetectorName() << endl;
      std::vector<Double_t> pos, sigma;
      Double_t theta=0;
      tie(pos,sigma) = PSDSpectrumAnalyzer(PSDSpectra.at(i));
      //theta = PSDMatrixAnalyzer(PSDMatrix.at(i));

      // I write it to file
      PSDoutput << experiment.GetDetectors().at(i)->GetDetectorName() << "\t" << pos.at(0) << "\t" << sigma.at(0) << "\t" << pos.at(1) << "\t" << sigma.at(1)<< "\t" << theta << endl;

      QDC1Spectra.at(i)->Write();
      QDC2Spectra.at(i)->Write();
      PSDSpectra.at(i)->Write();
      PSDMatrix.at(i)->Write();
    }
    else if(PSDSpectra.at(i)->GetEntries() !=0){
      cout << SetBOLD << SetForeRED << endl;
      cout << " PSD Spectrum for " <<experiment.GetDetectors().at(i)->GetDetectorName() << " is empty" << endl;
      cout << " Check data file(s) " << endl,
      cout << RESETTEXT << endl;
    }
  }
  outputfile->Close();
  PSDoutput.close();

  cout << "NRJ spectra are saved " << endl;

  // Printing of chronometer measurement
  localTimer.Stop();
  Double_t rtime2 = localTimer.RealTime();
  Double_t ctime2 = localTimer.CpuTime();
  std::cout << std::endl;
  std::cout << "End of Energy Spectra Plotting" << std::endl;
  std::cout << "# RealTime=" << rtime2 << " seconds, CpuTime="<< ctime2 << " seconds" <<std::endl;
  std::cout << std::endl;

  return 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int DrawAllParisUncalibratedSpectra_with_rotation(const CExperiment &experiment)
{
  //Definition du germe pour le tirage aleatoire
  srand48(time(NULL));

  //Declaration of all variables for names
  TString title,title2,title3,title4,title5,title6,title7;
  TString spectrumname,spectrumname2,spectrumname3,spectrumname4,spectrumname5,spectrumname6,spectrumname7;

  ROOT::EnableThreadSafety();

  // Declaration of energy spectra
  // First I count the number of PARIS
  int nbrparis=0;
  for(int d=0; d < (int)experiment.GetDetectors().size();d++)
  {
    if(experiment.GetDetectors().at(d)->GetDetectorType() == "PARIS") nbrparis++;
  }

  Int_t nbrofspectra = (int)experiment.GetDetectors().size();
  std::vector<TH1F*> PSDSpectra;
  std::vector<TH2F*> PSDMatrix;
  std::vector<TH2F*> PSDMatrixLaBr;
  std::vector<TH2F*> PSDMatrixNaI;
  std::vector<TH2F*> PSDMatrixrejected;
  std::vector<TH1F*> QDC1Spectra;
  std::vector<TH1F*> QDC2Spectra;
  Int_t highestdetlabel = 0;
  Int_t highestnbrchannels = 0;
  Int_t highestEmax = 0;
  Int_t nbrchannels = 0;
  Int_t Emin = 0;
  Int_t Emax = 0;

  // Defining all the energy spectra
  for(auto sindex = 0; sindex < nbrofspectra; sindex++)
  {
    //TH1F *localPSDSpectra;
    TH1F *localQDC1spectrum = nullptr;TH1F *localQDC2spectrum = nullptr; TH1F *localPSDSpectra = nullptr;
    TH2F *localPSDMatrix = nullptr;
    TH2F *localPSDMatrixLaBr= nullptr;
    TH2F *localPSDMatrixNaIcrosstalk = nullptr;
    //TH2F *localPSDMatrixBeyondLaBrNaI = nullptr;
    // Defining the title and name of the spectrum
    title = "LaBr Spectrum of detector ";
    title2 ="NaI Spectrum of detector ";
    title3 = " PSD Spectrum of detector ";
    title4 = " PSD Matrix of detector ";
    title5 = " PSD Matrix of CeBr "; 
    title6 = " PSD Matrix of NaI & crosstalk "; 
    //title7 = " PSD Matrix of the rest ";
    spectrumname = "nrjspectrum";
    spectrumname2 = "nrjspectrum2";
    spectrumname3 = "psdspectrum";
    spectrumname4 = "psdmatrix";
    spectrumname5 = "psdmatrixCeBr";
    spectrumname6 = "psdmatrixNaI";
    //spectrumname7 = "psdmatrixrejected";
    title += experiment.GetDetector(sindex)->GetDetectorName();
    title2 += experiment.GetDetector(sindex)->GetDetectorName();
    title3 += experiment.GetDetector(sindex)->GetDetectorName();
    title4 += experiment.GetDetector(sindex)->GetDetectorName();
    title5 += experiment.GetDetector(sindex)->GetDetectorName();
    title6 += experiment.GetDetector(sindex)->GetDetectorName();
    //title7 += experiment.GetDetector(sindex)->GetDetectorName();
    spectrumname +=experiment.GetDetector(sindex)->GetDetectorName();
    spectrumname2 +=experiment.GetDetector(sindex)->GetDetectorName();
    spectrumname3 +=experiment.GetDetector(sindex)->GetDetectorName();
    spectrumname4 +=experiment.GetDetector(sindex)->GetDetectorName();
    spectrumname5 +=experiment.GetDetector(sindex)->GetDetectorName();
    spectrumname6 +=experiment.GetDetector(sindex)->GetDetectorName();
    //spectrumname7 +=experiment.GetDetector(sindex)->GetDetectorName();

    if(experiment.GetDetector(sindex)->GetDetectorlabel()>highestdetlabel) highestdetlabel = experiment.GetDetectors().at(sindex)->GetDetectorlabel();

    // Getting the right range
    if(experiment.GetDetector(sindex)->GetDetectorType()!="RF")
    {
      nbrchannels = experiment.GetDetector(sindex)->GetNbrChannels();
      //cout << nbrchannels << "\t" << Emax << endl;
      if(nbrchannels > highestnbrchannels) highestnbrchannels=nbrchannels;
      Emax = experiment.GetDetector(sindex)->GetMaxchNumber();//GetMaxchNumber();
      if(Emax > highestEmax) highestEmax=Emax;
      //Emax = 200000;
      //nbrchannels = 200000;
      //if(experiment.GetDetector(sindex)->GetDetectorType()=="PARIS")
      {
        localQDC1spectrum = new TH1F(spectrumname,title,nbrchannels,Emin,Emax);
        localQDC2spectrum = new TH1F(spectrumname2,title2,nbrchannels,Emin,Emax);
        localPSDSpectra   = new TH1F(spectrumname3,title3,1600,0,TMath::Pi()/2.);
        localPSDMatrix    = new TH2F(spectrumname4,title4,nbrchannels,Emin,Emax,nbrchannels,Emin,Emax);
        localPSDMatrixLaBr    = new TH2F(spectrumname5,title5,nbrchannels,Emin,Emax,nbrchannels,Emin,Emax);
        localPSDMatrixNaIcrosstalk    = new TH2F(spectrumname6,title6,nbrchannels,Emin,Emax,nbrchannels,Emin,Emax);
        //localPSDMatrixBeyondLaBrNaI    = new TH2F(spectrumname7,title7,nbrchannels,Emin,Emax,nbrchannels,Emin,Emax);

        localQDC1spectrum->SetDirectory(0);
        localQDC2spectrum->SetDirectory(0);
        localPSDSpectra->SetDirectory(0);
        localPSDMatrix->SetDirectory(0);
        localPSDMatrixLaBr->SetDirectory(0);
        localPSDMatrixNaIcrosstalk->SetDirectory(0);
      }

      
      // Storing the NRJsectrum
      QDC1Spectra.push_back(localQDC1spectrum);
      QDC2Spectra.push_back(localQDC2spectrum);
      PSDSpectra.push_back(localPSDSpectra);
      PSDMatrix.push_back(localPSDMatrix);
      PSDMatrixLaBr.push_back(localPSDMatrixLaBr);
      PSDMatrixNaI.push_back(localPSDMatrixNaIcrosstalk);
      //PSDMatrixrejected.push_back(localPSDMatrixBeyondLaBrNaI);
    }

    // Cleaning the names for the next iteration
    spectrumname.Clear();spectrumname2.Clear();spectrumname3.Clear();spectrumname4.Clear();spectrumname5.Clear();spectrumname6.Clear();//spectrumname7.Clear();
    title.Clear();title2.Clear();title3.Clear();title4.Clear();title5.Clear();title6.Clear();//title7.Clear();
  }


  std::cout << FOREBLU << "We have " << nbrofspectra << " detectors among which " << nbrparis  << " PARIS phoswitches" << std::endl;
  std::cout << "We have created " << QDC1Spectra.size() << " QDC1 spectra among which " << nbrparis  << " PARIS phoswitches QDC1 spectra"  << std::endl;
  std::cout << "We have created " << QDC2Spectra.size() << " QDC2 spectra among which " << nbrparis  << " PARIS phoswitches QDC2 spectra"  << std::endl;
  std::cout << "We have created " << PSDSpectra.size() << " PSD spectra among which " << nbrparis  << " PARIS phoswitches PSD spectra"   << std::endl;
  std::cout << "We have created " << PSDMatrix.size() << " PSD Matrices among which " << nbrparis  << " PARIS phoswitches PSD Matrices"  << std::endl;
  std::cout << "We have created " << PSDMatrixLaBr.size() << " PSD CeBr Matrices among which " << nbrparis  << " PARIS CeBr3 PSD Matrices"  << std::endl;
  std::cout << "We have created " << PSDMatrixNaI.size() << " PSD NaI Matrices among which " << nbrparis  << " PARIS CeBr3 PSD Matrices"  << std::endl;
  //std::cout << "We have created " << PSDMatrixrejected.size() << " PSD beyond CebR and below NaI Matrices among which " << nbrparis  << " PARIS CeBr3 PSD Matrices"  << std::endl;
  std::cout << "We are Loaded .." << RESETTEXT << std::endl;

  // Creation of the file to save all the data
  TString outputfilename = experiment.GetFileDirectory_OUT();
  outputfilename += "3New_UncalibratedPARISspectraROTATED_all.root";
  TFile *outputfile = new TFile(outputfilename,"RECREATE");

  // Loading the TTree for reading the Data
  std::vector<TChain *> tab_chained_oak = experiment.GettheTChain();
  TChain *chained_oak = tab_chained_oak.at(0);

  label_Rawtype index;
  nrj_Rawtype enrj;

  // std::cout << FOREBLU << "Loading Tree ..." << std::endl;
  //auto cachesize = 500000000; // 500 MBytes
  //chained_oak -> SetCacheSize(cachesize);
  // chained_oak -> SetCacheLearnEntries(100000);*/
  // chained_oak -> SetBranchStatus("*",0);
  // chained_oak -> SetBranchAddress("label",&index);   // Detector number
  // chained_oak -> SetBranchAddress("nrj",&enrj);      // NRJ of the hit
  // std::cout << "Loaded .." << RESETTEXT << std::endl;

  // Another Timer
  TStopwatch localTimer;

  //Starting of the chronometer
  localTimer.Reset();
  localTimer.Start();

  //---------------------------------------------------------------------------//
  //                                                                           //
  //                            Reading the TTree                              //
  //                                                                           //
  //---------------------------------------------------------------------------//
  // // Determination of the number of event that have to be treated
  // const auto max_workers = std::thread::hardware_concurrency();
  auto n_workers = experiment.GetThreadsNbr();
  // if(n_workers > max_workers){
  //   cout << "Current machine only supports " << max_workers << " parallel threads. I will be using this number" << endl;
  //   n_workers = max_workers;
  // }
  const auto chainentries = chained_oak -> GetEntries();
  // const auto range = chainentries/n_workers;
  // int percent = (int)(0.05*chainentries);

  std::cout << "Beginning of Energy spectra building research" << std::endl;
  std::cout << "On " << chainentries << " entries" << std::endl;

  // Creating a ProgressBar which is nice in MT mode
  using namespace indicators;
  ProgressBar bar{
    option::BarWidth{60},
    option::Start{"["},
    option::Fill{"="},
    option::Lead{">"},
    option::Remainder{"-"},
    option::End{"]"},
    option::PostfixText{"Reading"},
    option::ShowElapsedTime{true},
    option::ShowRemainingTime{true},
    option::ForegroundColor{Color::grey},
    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
  };

  // I read all the file of the TChain
  int fileindex  = 0;
  int filenumber =  experiment.GetDataFileNames().size();
  mutex forfilling;
  for(auto &file : experiment.GetDataFileNames())
  {
    fileindex++;
    //PrintVector(experiment.GetDataFileNames());
    bar.set_progress((int)((double)fileindex/(double)filenumber*100.));
    // Defining a TTreeProcessor object to handle the TTree in the File in a MT mode
    ROOT::TTreeProcessorMT TP(file, experiment.GetDataTreeName(), n_workers);

    // Scanning TTree
    // Launch the parallel processing of the tree
    int threadnbr = 0;
    auto loop_and_fill = [&] (TTreeReader &myReader){
      TTreeReaderValue<label_Rawtype> labelRV(myReader,"label");
      TTreeReaderValue<nrj_Rawtype> QDC1RV(myReader,"nrj");
      TTreeReaderValue<nrj_Rawtype> QDC2RV(myReader,"nrj2");

      ULong64_t hitnumber = 0;
      while(myReader.Next())
      {
        auto label  = *labelRV;
        auto NRJ    = *QDC1RV;  // Short Gate
        auto NRJ2   = *QDC2RV;  // Long Gate

        // I define a new hit a fill in the information
        hitnumber++;threadnbr++;
        if(label > 19) // To make sure I only consider PARIS detectors
        {
          CHit *hit = new CHit(hitnumber+1000*threadnbr);
          hit->SetHit(label, 0, NRJ, NRJ2, false);

          int spectrumindex = experiment.GetLabel2Detnbrs(static_cast<int>(label));
          
          Double_t PSD = hit->PerformPARISPSD();// PSD = atan(short/long)
          //std::cout << "All good" << std::endl;
          Bool_t isLaBr = experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->IsPureLaBr3(PSD,hit->GetHitE1(),hit->GetHitE2()); // I need the Long charge to get proper LaBr3 selection
          Bool_t isBeyondLaBr3andNaI = experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->IsBeyondLaBr3andNaI(PSD,hit->GetHitE1(),hit->GetHitE2()); // I need the Long charge to get proper LaBr3 selection
          Bool_t isNaI = experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->IsPureNaI(PSD,hit->GetHitE1(),hit->GetHitE2()); // I need the Long charge to get proper NaI selection
          //cout<<isLaBr<<endl;
          //cout << " Detector " << (string)experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->GetDetectorName() << " discriLaBr_pos = " << experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->GetLaBrDiscriPosition() << endl;
          // Filling up the spectra
          //forfilling.lock();
          double NRJ_bf_ROT = hit->GetHitE1();
          double NRJ2_bf_ROT = hit->GetHitE2();
          if(experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->GetDetectorType() == "PARIS")
          {
            //cout << "Qs = " << NRJ2 << "; Ql = " << NRJ << endl;
            tie(NRJ2,NRJ) = Rotation(
                                      static_cast<double>(experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->GetRotationAngle()),
                                      static_cast<double>(experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->GetRotationAngleTan()),
                                      NRJ2,
                                      NRJ);
             //tie(NRJ2,NRJ) = Rotation(0.,0.,NRJ2,NRJ);
             //cout << "ROT_Qs = " << NRJ2 << "; ROT_Ql = " << NRJ << endl;
          }
          if(spectrumindex < PSDSpectra.size())
          {
            // cout << endl;
            // hit->PrintHit();
            // cout << "Hit # " << hitnumber+1000*threadnbr << "; Filling Spectra " << spectrumindex << endl;
            PSDSpectra.at(spectrumindex)  -> Fill(PSD);
            if(isLaBr && !isNaI) QDC1Spectra.at(spectrumindex) -> Fill(NRJ_bf_ROT);
            if(isNaI) QDC2Spectra.at(spectrumindex) -> Fill(NRJ2);
          }
          
          if(spectrumindex < PSDMatrix.size())
          {
            //cout << "Hit # " << hitnumber+1000*threadnbr << "; Filling Matrix " << spectrumindex << endl;
            if(isLaBr && !isNaI)PSDMatrixLaBr.at(spectrumindex)->Fill((Double_t)NRJ2,(Double_t)NRJ);
            //if(isNaI)PSDMatrixNaI.at(spectrumindex)->Fill((Double_t)NRJ2,(Double_t)NRJ);
            if(!isLaBr && !isBeyondLaBr3andNaI ) PSDMatrixNaI.at(spectrumindex)->Fill((Double_t)NRJ2,(Double_t)NRJ);
            if(!isBeyondLaBr3andNaI) PSDMatrix.at(spectrumindex)->Fill((Double_t)NRJ2,(Double_t)NRJ);
          }
          //forfilling.unlock();

          delete hit;
        }
      }
    };
    TP.Process(loop_and_fill);
  }

  // All spectra has been built
  // Now analyzing the PSDSpectra to print to text files the parameters of PSD
  std::cout << endl << RESETTEXT << "Saving in " << outputfilename << std::endl;
  outputfile->cd();

  // Now I have the peak position for LaBr3 selection
  TString PSDoutputfilename = "PSDParameter1_PARIS.txt";
  ofstream PSDoutput(PSDoutputfilename, ios::out);
  PSDoutput << "Det Name \t LaBrPos \t LaBrSigma \t NaIPos \t NaISigma \t theta"<< endl;
  for(int i = 0; i < nbrofspectra; i++)
  {
    if(experiment.GetDetectors().at(i)->GetDetectorType()=="PARIS" && PSDSpectra.at(i)->GetEntries() !=0)
    {
      cout << "Saving Detector " << experiment.GetDetectors().at(i)->GetDetectorName() << endl;
      std::vector<Double_t> pos, sigma;
      Double_t theta;
      tie(pos,sigma) = PSDSpectrumAnalyzer(PSDSpectra.at(i));
      //theta = PSDMatrixAnalyzer(PSDMatrix.at(i));

      // I write it to file
      PSDoutput << experiment.GetDetectors().at(i)->GetDetectorName() << "\t" << pos.at(0) << "\t" << sigma.at(0) << "\t" << pos.at(1) << "\t" << sigma.at(1)<< "\t" << theta << endl;

      QDC1Spectra.at(i)->Write();
      QDC2Spectra.at(i)->Write();
      PSDSpectra.at(i)->Write();
      PSDMatrix.at(i)->Write();
      PSDMatrixLaBr.at(i)->Write();
      PSDMatrixNaI.at(i)->Write();
      //PSDMatrixrejected.at(i)->Write();

    }
    else if(PSDSpectra.at(i)->GetEntries() !=0){
      cout << SetBOLD << SetForeRED << endl;
      cout << " PSD Spectrum for " <<experiment.GetDetectors().at(i)->GetDetectorName() << " is empty" << endl;
      cout << " Check data file(s) " << endl,
      cout << RESETTEXT << endl;
    }
  }
  outputfile->Close();
  PSDoutput.close();

  cout << "NRJ spectra are saved " << endl;
  // Clean up histograms
  for (auto hist : QDC1Spectra) delete hist;
  for (auto hist : QDC2Spectra) delete hist;
  for (auto hist : PSDSpectra) delete hist;
  for (auto hist : PSDMatrix) delete hist;
  for (auto hist : PSDMatrixLaBr) delete hist;
  for (auto hist : PSDMatrixNaI) delete hist;

  // Clear vectors
  QDC1Spectra.clear();
  QDC2Spectra.clear();
  PSDSpectra.clear();
  PSDMatrix.clear();
  PSDMatrixLaBr.clear();
  PSDMatrixNaI.clear();

  // Printing of chronometer measurement
  localTimer.Stop();
  Double_t rtime2 = localTimer.RealTime();
  Double_t ctime2 = localTimer.CpuTime();
  std::cout << std::endl;
  std::cout << "End of Energy Spectra Plotting" << std::endl;
  std::cout << "# RealTime=" << rtime2 << " seconds, CpuTime="<< ctime2 << " seconds" <<std::endl;
  std::cout << std::endl;

  return 1;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int DrawAllParisCalibratedSpectra(const CExperiment &experiment)
{
  //Definition du germe pour le tirage aleatoire
  srand48(time(NULL));

  //Declaration of all variables for names
  TString title,title2,title3,title4;
  TString spectrumname,spectrumname2,spectrumname3,spectrumname4;

  ROOT::EnableThreadSafety();

  // Declaration of energy spectra
  // First I count the number of PARIS
  int nbrparis=0;
  for(int d=0; d < (int)experiment.GetDetectors().size();d++)
  {
    if(experiment.GetDetectors().at(d)->GetDetectorType() == "PARIS") nbrparis++;
  }

  Int_t nbrofspectra = (int)experiment.GetDetectors().size();
  //std::vector<TH1F*> PSDSpectra;
  //std::vector<TH2F*> PSDMatrix;
  std::vector<TH1F*> QDC1Spectra;
  std::vector<TH1F*> QDC2Spectra;
  Int_t highestdetlabel = 0;
  Int_t highestnbrchannels = 0;
  Int_t highestEmax = 0;
  Int_t nbrchannels = 0;
  Int_t Emin = 0;
  Int_t Emax = 20000;

  // Defining all the energy spectra
  for(int sindex = 0; sindex < nbrofspectra; sindex++)
  {
    //TH1F *localPSDSpectra;
    TH1F *localQDC1spectrum,*localQDC2spectrum;//,*localPSDSpectra;
    //TH2F *localPSDMatrix;
    // Defining the title and name of the spectrum
    title = "LaBr Spectrum of detector ";title2 ="NaI Spectrum of detector ";title3 = " PSD Spectrum of detector ";title4 = " PSD Matrix of detector ";
    spectrumname = "nrjspectrum";spectrumname2 = "nrjspectrum2";spectrumname3 = "psdspectrum";spectrumname4 = "psdmatrix";
    title += experiment.GetDetector(sindex)->GetDetectorName();title2 += experiment.GetDetector(sindex)->GetDetectorName();title3 += experiment.GetDetector(sindex)->GetDetectorName();title4 += experiment.GetDetector(sindex)->GetDetectorName();
    spectrumname +=experiment.GetDetector(sindex)->GetDetectorName();spectrumname2 +=experiment.GetDetector(sindex)->GetDetectorName();spectrumname3 +=experiment.GetDetector(sindex)->GetDetectorName();spectrumname4 +=experiment.GetDetector(sindex)->GetDetectorName();

    if(experiment.GetDetector(sindex)->GetDetectorlabel()>highestdetlabel) highestdetlabel = experiment.GetDetectors().at(sindex)->GetDetectorlabel();

    // Getting the right range
    if(experiment.GetDetector(sindex)->GetDetectorType()!="RF")
    {
      nbrchannels = (Emax-Emin)/2;
      //cout << nbrchannels << "\t" << Emax << endl;
      if(nbrchannels > highestnbrchannels) highestnbrchannels=nbrchannels;
      Emax = 20000; // keV
      if(Emax > highestEmax) highestEmax=Emax;
      localQDC1spectrum = new TH1F(spectrumname,title,nbrchannels,Emin,Emax);
      localQDC2spectrum = new TH1F(spectrumname2,title2,nbrchannels,Emin,Emax);
      //localPSDSpectra   = new TH1F(spectrumname3,title3,1600,-1,1);
      //localPSDMatrix    = new TH2F(spectrumname4,title4,nbrchannels,Emin,Emax,nbrchannels,Emin,Emax);

      
      // Storing the NRJsectrum
      QDC1Spectra.push_back(localQDC1spectrum);
      QDC2Spectra.push_back(localQDC2spectrum);
      //PSDSpectra.push_back(localPSDSpectra);
      //PSDMatrix.push_back(localPSDMatrix);
      
    }

    // Cleaning the names for the next iteration
    spectrumname.Clear();spectrumname2.Clear();spectrumname3.Clear();spectrumname4.Clear();
    title.Clear();title2.Clear();title3.Clear();title4.Clear();
  }


  std::cout << FOREBLU << "We have " << nbrofspectra << " detectors among which " << nbrparis  << " PARIS phoswitches" << std::endl;
  std::cout << "We have created " << QDC1Spectra.size() << " QDC1 spectra among which " << nbrparis  << " PARIS phoswitches QDC1 spectra"  << std::endl;
  std::cout << "We have created " << QDC2Spectra.size() << " QDC2 spectra among which " << nbrparis  << " PARIS phoswitches QDC2 spectra"  << std::endl;
  //std::cout << "We have created " << PSDSpectra.size() << " PSD spectra among which " << nbrparis  << " PARIS phoswitches PSD spectra"   << std::endl;
  //std::cout << "We have created " << PSDMatrix.size() << " PSD Matrices among which " << nbrparis  << " PARIS phoswitches PSD Matrices"  << std::endl;
  std::cout << "We are Loaded .." << RESETTEXT << std::endl;

  // Creation of the file to save all the data
  TString outputfilename = experiment.GetFileDirectory_OUT();
  outputfilename += "CalibratedPARISspectra_all.root";
  TFile *outputfile = new TFile(outputfilename,"RECREATE");

  // Loading the TTree for reading the Data
  std::vector<TChain *> tab_chained_oak = experiment.GettheTChain();
  TChain *chained_oak = tab_chained_oak.at(0);

  label_Rawtype index;
  nrj_Rawtype enrj;

  // std::cout << FOREBLU << "Loading Tree ..." << std::endl;
  // /*auto cachesize = 10000000; // 10 MBytes
  // chained_oak -> SetCacheSize(cachesize);
  // chained_oak -> SetCacheLearnEntries(100000);*/
  // chained_oak -> SetBranchStatus("*",0);
  // chained_oak -> SetBranchAddress("label",&index);   // Detector number
  // chained_oak -> SetBranchAddress("nrj",&enrj);      // NRJ of the hit
  // std::cout << "Loaded .." << RESETTEXT << std::endl;

  // Another Timer
  TStopwatch localTimer;

  //Starting of the chronometer
  localTimer.Reset();
  localTimer.Start();

  //---------------------------------------------------------------------------//
  //                                                                           //
  //                            Reading the TTree                              //
  //                                                                           //
  //---------------------------------------------------------------------------//
  // // Determination of the number of event that have to be treated
  // const auto max_workers = std::thread::hardware_concurrency();
  auto n_workers = experiment.GetThreadsNbr();
  // if(n_workers > max_workers){
  //   cout << "Current machine only supports " << max_workers << " parallel threads. I will be using this number" << endl;
  //   n_workers = max_workers;
  // }
  const auto chainentries = chained_oak -> GetEntries();
  // const auto range = chainentries/n_workers;
  // int percent = (int)(0.05*chainentries);

  std::cout << "Beginning of Energy spectra building research" << std::endl;
  std::cout << "On " << chainentries << " entries" << std::endl;

  // Creating a ProgressBar which is nice in MT mode
  using namespace indicators;
  ProgressBar bar{
    option::BarWidth{60},
    option::Start{"["},
    option::Fill{"="},
    option::Lead{">"},
    option::Remainder{"-"},
    option::End{"]"},
    option::PostfixText{"Reading"},
    option::ShowElapsedTime{true},
    option::ShowRemainingTime{true},
    option::ForegroundColor{Color::grey},
    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
  };

  // I read all the file of the TChain
  int fileindex  = 0;
  int filenumber =  experiment.GetDataFileNames().size();
  mutex forfilling;
  for(auto &file : experiment.GetDataFileNames())
  {
    fileindex++;
    //PrintVector(experiment.GetDataFileNames());
    bar.set_progress((int)((double)fileindex/(double)filenumber*100.));
    // Defining a TTreeProcessor object to handle the TTree in the File in a MT mode
    ROOT::TTreeProcessorMT TP(file, experiment.GetDataTreeName(), n_workers);

    // Scanning TTree
    // Launch the parallel processing of the tree
    int threadnbr = 0;
    auto loop_and_fill = [&] (TTreeReader &myReader){
      TTreeReaderValue<label_Rawtype> labelRV(myReader,"label");
      TTreeReaderValue<nrj_type> QDC1RV(myReader,"nrj");
      TTreeReaderValue<nrj_type> QDC2RV(myReader,"nrj2");

      ULong64_t hitnumber = 0;
      while(myReader.Next())
      {
        auto label  = *labelRV;
        auto NRJ    = *QDC1RV;  // Short Gate
        auto NRJ2   = *QDC2RV;  // Long Gate

        // I define a new hit a fill in the information
        hitnumber++;threadnbr++;
        if(label > 0) // To make sure I only consider PARIS detectors
        {
          CHit *hit = new CHit(hitnumber+1000*threadnbr);
          hit->SetHit(label, 0, NRJ, NRJ2, false);

          int spectrumindex = experiment.GetLabel2Detnbrs(static_cast<int>(label));
          
          Bool_t isLaBr = kFALSE;
          if(NRJ2 == 0) isLaBr = kTRUE;

          // Filling up the spectra
          //forfilling.lock();
          if(spectrumindex < QDC1Spectra.size())
          {
            // cout << endl;
            // hit->PrintHit();
            // cout << "Hit # " << hitnumber+1000*threadnbr << "; Filling Spectra " << spectrumindex << endl;
            if(isLaBr) QDC1Spectra.at(spectrumindex) -> Fill(NRJ);
            QDC2Spectra.at(spectrumindex) -> Fill(NRJ2);
          }
          //forfilling.unlock();

          delete hit;
        }
      }
    };
    TP.Process(loop_and_fill);
  }

  // All spectra has been built
  // Now analyzing the PSDSpectra to print to text files the parameters of PSD
  std::cout << endl << RESETTEXT << "Saving in " << outputfilename << std::endl;
  outputfile->cd();

  // Now I have the peak position for LaBr3 selection
  for(int i = 0; i < nbrofspectra; i++)
  {
    if(experiment.GetDetectors().at(i)->GetDetectorType()=="PARIS" && QDC1Spectra.at(i)->GetEntries()!=0)
    {
      cout << "Saving Detector " << experiment.GetDetectors().at(i)->GetDetectorName() << endl;
      QDC1Spectra.at(i)->Write();
      QDC2Spectra.at(i)->Write();
    }
  }
  outputfile->Close();
  
  cout << "NRJ spectra are saved " << endl;

  // Printing of chronometer measurement
  localTimer.Stop();
  Double_t rtime2 = localTimer.RealTime();
  Double_t ctime2 = localTimer.CpuTime();
  std::cout << std::endl;
  std::cout << "End of Energy Spectra Plotting" << std::endl;
  std::cout << "# RealTime=" << rtime2 << " seconds, CpuTime="<< ctime2 << " seconds" <<std::endl;
  std::cout << std::endl;

  return 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


int DrawAllCalibrationSpectra(const CExperiment &experiment)
{
  //Definition du germe pour le tirage aleatoire
  srand48(time(NULL));

  //Declaration of all variables for names
  TString title;
  TString spectrumname;

  ROOT::EnableThreadSafety();

  // Declaration of time spectra
  Int_t nbrofspectra = experiment.GetDetectors().size();
  std::vector<TH1F*> NRJspectra;
  Int_t highestdetlabel(0),highestnbrchannels(0),highestEmax(0);
  Int_t nbrchannels(0);
  Int_t Emin(0);
  Int_t Emax(0);
  int binnumber_MOSAHR = static_cast<int>(TMath::Power(2,16));
  int binnumber_CARAS = static_cast<int>(TMath::Power(2,14));


  // Defining all the energy spectra
  for(int sindex = 0; sindex < nbrofspectra; sindex++)
  {
    TH1F* localNRJspectrum;
    // Defining the title and name of the spectrum
    title = "Energy Spectrum of detector ";
    spectrumname = "nrjspectrum";
    title += experiment.GetDetector(sindex)->GetDetectorName();
    spectrumname +=experiment.GetDetector(sindex)->GetDetectorName();
    if(experiment.GetDetector(sindex)->GetDetectorlabel()>highestdetlabel) highestdetlabel = experiment.GetDetectors().at(sindex)->GetDetectorlabel();

    // Getting the right range
    if(experiment.GetDetector(sindex)->GetDetectorType()!="RF")
    {
      nbrchannels = experiment.GetDetector(sindex)->GetNbrChannels();
      //cout << nbrchannels << "\t" << Emax << endl;
      if(nbrchannels > highestnbrchannels) highestnbrchannels=nbrchannels;
      Emax = experiment.GetDetector(sindex)->GetMaxchNumber();//GetMaxchNumber();
      if(Emax > highestEmax) highestEmax=Emax;
      localNRJspectrum = new TH1F(spectrumname,title,nbrchannels,Emin,Emax);

      // Storing the NRJsectrum
      NRJspectra.push_back(localNRJspectrum);
    }

    // Cleaning the names for the next iteration
    spectrumname.Clear();
    title.Clear();
  }

  // I declare la Time matrix to check alignement later
  TH2F* NRJmatrix = new TH2F("NRJalignementmatrix","NRJ spectra of all detectors",highestdetlabel,1,highestdetlabel,binnumber_MOSAHR,Emin,highestEmax);

  //Declaration of all variables for names
  TString title2,title3,title4;
  TString spectrumname2,spectrumname3,spectrumname4;

  ROOT::EnableThreadSafety();

  // Declaration of energy spectra
  // First I count the number of PARIS
  int nbrparis=0;
  for(int d=0; d < (int)experiment.GetDetectors().size();d++)
  {
    if(experiment.GetDetectors().at(d)->GetDetectorType() == "PARIS") nbrparis++;
  }

  nbrofspectra = (int)experiment.GetDetectors().size();
  std::vector<TH1F*> PSDSpectra;
  std::vector<TH2F*> PSDMatrix;
  std::vector<TH1F*> QDC1Spectra;
  std::vector<TH1F*> QDC2Spectra;
  //Int_t highestdetlabel(0),highestnbrchannels(0),highestEmax(0);
  ///Int_t nbrchannels(0);
  //Int_t Emin(0);
  //Int_t Emax(0);

  // Defining all the energy spectra
  for(int sindex = 0; sindex < nbrofspectra; sindex++)
  {
    //TH1F *localPSDSpectra;
    TH1F *localQDC1spectrum,*localQDC2spectrum,*localPSDSpectra;
    TH2F *localPSDMatrix;
    // Defining the title and name of the spectrum
    title = "LaBr Spectrum of detector ";title2 ="NaI Spectrum of detector ";title3 = " PSD Spectrum of detector ";title4 = " PSD Matrix of detector ";
    spectrumname = "qdcspectrum";spectrumname2 = "qdcspectrum2";spectrumname3 = "psdspectrum";spectrumname4 = "psdmatrix";
    title += experiment.GetDetector(sindex)->GetDetectorName();title2 += experiment.GetDetector(sindex)->GetDetectorName();title3 += experiment.GetDetector(sindex)->GetDetectorName();title4 += experiment.GetDetector(sindex)->GetDetectorName();
    spectrumname +=experiment.GetDetector(sindex)->GetDetectorName();spectrumname2 +=experiment.GetDetector(sindex)->GetDetectorName();spectrumname3 +=experiment.GetDetector(sindex)->GetDetectorName();spectrumname4 +=experiment.GetDetector(sindex)->GetDetectorName();

    if(experiment.GetDetector(sindex)->GetDetectorlabel()>highestdetlabel) highestdetlabel = experiment.GetDetectors().at(sindex)->GetDetectorlabel();

    // Getting the right range
    if(experiment.GetDetector(sindex)->GetDetectorType()!="RF")
    {
      nbrchannels = experiment.GetDetector(sindex)->GetNbrChannels();
      //cout << nbrchannels << "\t" << Emax << endl;
      if(nbrchannels > highestnbrchannels) highestnbrchannels=nbrchannels;
      Emax = experiment.GetDetector(sindex)->GetMaxchNumber();//GetMaxchNumber();
      if(Emax > highestEmax) highestEmax=Emax;
      localQDC1spectrum = new TH1F(spectrumname,title,nbrchannels,Emin,Emax);
      localQDC2spectrum = new TH1F(spectrumname2,title2,nbrchannels,Emin,Emax);
      localPSDSpectra   = new TH1F(spectrumname3,title3,1600,-1,1);
      localPSDMatrix    = new TH2F(spectrumname4,title4,nbrchannels,Emin,Emax,nbrchannels,Emin,Emax);

      // Storing the NRJsectrum
      QDC1Spectra.push_back(localQDC1spectrum);
      QDC2Spectra.push_back(localQDC2spectrum);
      PSDSpectra.push_back(localPSDSpectra);
      PSDMatrix.push_back(localPSDMatrix);
    }

    // Cleaning the names for the next iteration
    spectrumname.Clear();spectrumname2.Clear();spectrumname3.Clear();spectrumname4.Clear();
    title.Clear();title2.Clear();title3.Clear();title4.Clear();
  }


  std::cout << FOREBLU << "We have " << nbrofspectra << " detectors" << std::endl;
  std::cout << "We have created " << QDC1Spectra.size() << " QDC1 spectra" << std::endl;
  std::cout << "We have created " << QDC2Spectra.size() << " QDC2 spectra" << std::endl;
  std::cout << "We have created " << PSDSpectra.size() << " PSD spectra" << std::endl;
  std::cout << "We have created " << PSDMatrix.size() << " PSD Matrices" << std::endl;
  std::cout << "We are Loaded .." << RESETTEXT << std::endl;

  // Creation of the file to save all the data
  TString outputfilename = experiment.GetFileDirectory_OUT();
  outputfilename += "UncalibratedEnergyspectra_all.root";
  TFile *outputfile = new TFile(outputfilename,"RECREATE");


  // Loading the TTree for reading the Data
  std::vector<TChain *> tab_chained_oak = experiment.GettheTChain();
  TChain *chained_oak = tab_chained_oak.at(0);

  label_Rawtype index;
  nrj_Rawtype enrj;

  // std::cout << FOREBLU << "Loading Tree ..." << std::endl;
  // /*auto cachesize = 10000000; // 10 MBytes
  // chained_oak -> SetCacheSize(cachesize);
  // chained_oak -> SetCacheLearnEntries(100000);*/
  // chained_oak -> SetBranchStatus("*",0);
  // chained_oak -> SetBranchAddress("label",&index);   // Detector number
  // chained_oak -> SetBranchAddress("nrj",&enrj);      // NRJ of the hit
  // std::cout << "Loaded .." << RESETTEXT << std::endl;

  // Another Timer
  TStopwatch localTimer;

  //Starting of the chronometer
  localTimer.Reset();
  localTimer.Start();

  //---------------------------------------------------------------------------//
  //                                                                           //
  //                            Reading the TTree                              //
  //                                                                           //
  //---------------------------------------------------------------------------//
  // // Determination of the number of event that have to be treated
  // const auto max_workers = std::thread::hardware_concurrency();
  auto n_workers = experiment.GetThreadsNbr();
  // if(n_workers > max_workers){
  //   cout << "Current machine only supports " << max_workers << " parallel threads. I will be using this number" << endl;
  //   n_workers = max_workers;
  // }
  const auto chainentries = chained_oak -> GetEntries();
  // const auto range = chainentries/n_workers;
  // int percent = (int)(0.05*chainentries);

  std::cout << "Beginning of Energy spectra building research" << std::endl;
  std::cout << "On " << chainentries << " entries" << std::endl;

  // Creating a ProgressBar which is nice in MT mode
  using namespace indicators;
  ProgressBar bar{
    option::BarWidth{60},
    option::Start{"["},
    option::Fill{"="},
    option::Lead{">"},
    option::Remainder{"-"},
    option::End{"]"},
    option::PostfixText{"Reading"},
    option::ShowElapsedTime{true},
    option::ShowRemainingTime{true},
    option::ForegroundColor{Color::grey},
    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
  };

  // I read all the file of the TChain
  int fileindex  = 0;
  int filenumber =  experiment.GetDataFileNames().size();
  mutex forfilling;
  
  for(auto &file : experiment.GetDataFileNames())
  {
    fileindex++;
    //PrintVector(experiment.GetDataFileNames());
    bar.set_progress((int)((double)fileindex/(double)filenumber*100.));
    // Defining a TTreeProcessor object to handle the TTree in the File in a MT mode
    ROOT::TTreeProcessorMT TP(file, experiment.GetDataTreeName(), n_workers);

    // Scanning TTree
    // Launch the parallel processing of the tree
    int threadnbr = 0;
    auto loop_and_fill = [&] (TTreeReader &myReader){
      // std::random_device rd;
      // std::mt19937 gen(rd());
      // std::uniform_real_distribution<> dis(0.0,1.0);
      
      TTreeReaderValue<label_Rawtype> labelRV(myReader,"label");
      TTreeReaderValue<nrj_Rawtype> QDC1RV(myReader,"nrj");
      TTreeReaderValue<nrj_Rawtype> QDC2RV(myReader,"nrj2");

      ULong64_t hitnumber = 0;
      while(myReader.Next())
      {
        auto label  = *labelRV;
        auto NRJ    = *QDC1RV;  // Short Gate
        auto NRJ2   = *QDC2RV;  // Long Gate

        // I define a new hit a fill in the information
        hitnumber++;threadnbr++;
        if(label < 200)
        {
          int spectrumindex = experiment.GetLabel2Detnbrs(static_cast<int>(label));
          int orginalBits = experiment.GetDetector(static_cast<int>(spectrumindex))->GetMaxchNumber();
          //Double_t NRJ = enrj;//Faster2bitsNRJConverter(enrj, experiment,spectrumindex); // To convert to a sort of integer (it should remove the "non-linearity")
          //auto NRJ_compressed = CompressFASTERValue(NRJ,orginalBits,binnumber_MOSAHR);

          //auto alea = dis(gen);
          
          // Calculating the right energy
          if(NRJ > 10) NRJspectra.at(spectrumindex)->Fill(NRJ);

          // And the energy Matrix
          if (NRJ > 10) NRJmatrix->Fill(label,NRJ);
        }
        if(label > 0) // To make sure I only consider PARIS detectors
        {
          int spectrumindex = experiment.GetLabel2Detnbrs(static_cast<int>(label));
          int orginalBits = experiment.GetDetector(static_cast<int>(spectrumindex))->GetMaxchNumber();
          //auto NRJ_compressed = CompressFASTERValue(NRJ,orginalBits,binnumber_CARAS);
          CHit *hit = new CHit(hitnumber+1000*threadnbr);
          hit->SetHit(label, 0, NRJ, NRJ2, false);
          auto NRJmemory = NRJ;//NRJ_compressed;
          //int spectrumindex = experiment.GetLabel2Detnbrs(static_cast<int>(label));
          Double_t PSD = hit->PerformPARISPSD();// PSD = atan
          if(experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->GetDetectorType() == "PARIS")
          {
            tie(NRJ2,NRJ) = Rotation((double)experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->GetRotationAngle(),(double)experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->GetRotationAngleTan(),NRJ2,NRJ);
          }
          
          Bool_t isLaBr = experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->IsPureLaBr3(PSD,hit->GetHitE1(),hit->GetHitE2());
          Bool_t isBeyondLaBrandNaI = experiment.GetDetectors().at(experiment.GetLabel2Detnbrs(static_cast<int>(label)))->IsBeyondLaBr3andNaI(PSD,hit->GetHitE1(),hit->GetHitE2());
          
          // Filling up the spectra
          //forfilling.lock();
          if(spectrumindex < PSDSpectra.size())
          {
            // cout << endl;
            // hit->PrintHit();
            // cout << "Hit # " << hitnumber+1000*threadnbr << "; Filling Spectra " << spectrumindex << endl;
            PSDSpectra.at(spectrumindex)  -> Fill(PSD);
            if(isLaBr) QDC1Spectra.at(spectrumindex) -> Fill(NRJmemory);
            else if(!isLaBr && !isBeyondLaBrandNaI){QDC2Spectra.at(spectrumindex) -> Fill(NRJ2);}
          }
          if(spectrumindex < PSDMatrix.size())
          {
            //cout << "Hit # " << hitnumber+1000*threadnbr << "; Filling Matrix " << spectrumindex << endl;
            if(!isBeyondLaBrandNaI) PSDMatrix.at(spectrumindex)->Fill((Double_t)NRJ2,(Double_t)NRJ);
          }
          //forfilling.unlock();

          delete hit;
        }
      }
    };
    TP.Process(loop_and_fill);
  }

  // All spectra has been built
  // Now analyzing the PSDSpectra to print to text files the parameters of PSD
  std::cout << endl << RESETTEXT << "Saving in " << outputfilename << std::endl;
  outputfile->cd();

  // Now I have the peak position for LaBr3 selection
  //TString PSDoutputfilename = "PSDParameter_PARIS.txt";
  //ofstream PSDoutput(PSDoutputfilename, ios::out);
  //PSDoutput << "Det Name \t LaBrPos \t LaBrSigma \t NaIPos \t NaISigma \t theta"<< endl;
  for(int i = 0; i < nbrofspectra; i++)
  {
    if(experiment.GetDetectors().at(i)->GetDetectorType()=="PARIS")
    {
      //std::vector<Double_t> pos;std::vector<Double_t> sigma;Double_t theta(0.);
      //tie(pos,sigma) = PSDSpectrumAnalyzer(PSDSpectra.at(i));
      //theta = PSDMatrixAnalyzer(PSDMatrix.at(i));

      // I write it to file
      //PSDoutput << experiment.GetDetectors().at(i)->GetDetectorName() << "\t" << pos.at(0) << "\t" << sigma.at(0) << "\t" << pos.at(1) << "\t" << sigma.at(1)<< "\t" << theta << endl;

      QDC1Spectra.at(i)->Write();
      QDC2Spectra.at(i)->Write();
      PSDSpectra.at(i)->Write();
      PSDMatrix.at(i)->Write();
    }
  }

  for(int i = 0; i < (int) NRJspectra.size(); i++) NRJspectra.at(i)->Write();

  NRJmatrix->Write();


  outputfile->Close();

  cout << "NRJ spectra are saved " << endl;

  // Printing of chronometer measurement
  localTimer.Stop();
  Double_t rtime2 = localTimer.RealTime();
  Double_t ctime2 = localTimer.CpuTime();
  std::cout << std::endl;
  std::cout << "End of Energy Spectra Plotting" << std::endl;
  std::cout << "# RealTime=" << rtime2 << " seconds, CpuTime="<< ctime2 << " seconds" <<std::endl;
  std::cout << std::endl;

  return 1;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::tuple<std::vector<Double_t>,std::vector<Double_t>> PSDSpectrumAnalyzer(TH1F *spectrum)
{
  std::vector<Double_t> vec_pos;Double_t pos = 0.;
  std::vector<Double_t> vec_sigma;Double_t sigma = 0.02;
  TSpectrum *s = new TSpectrum(6);
  cout << "Searching for peaks in " << spectrum->GetName() << endl;
  int nfound = s->Search(spectrum,0.2,"",0.15);
  cout << nfound << " Peaks found, now fitting" << endl;
  if(nfound >= 2)
  {
    Double_t *xpeaks = s->GetPositionX();
    Double_t *ypeaks = s->GetPositionY();
    std::vector<Double_t> xpeak;
    std::vector<Double_t> ypeak;
    for(auto p = 0; p<nfound;p++)
    {
      xpeak.push_back(xpeaks[p]);
      ypeak.push_back(ypeaks[p]);
    }
    sortVectors(xpeak, less<Double_t>(), xpeak,ypeak);

    // Je ne fit que le deuxieme qui ne doit pas être loin de 1 // Sur PARIS il Correspond a LaBr3
    cout << "FITTING PEAK AT " << xpeak.at(1) << endl;
    Double_t xp = xpeak.at(1);
    TF1 *f1 = new TF1("peak", "gaus", xp-(2*sigma), xp+(2*sigma));
    f1->SetParameter(1,xp);
    //f1->SetParLimits(1,xp-sigma,xp+sigma);
    // Then I get amplitude
    Int_t bin = spectrum->GetXaxis()->FindBin(xp);
    Double_t yp = spectrum->GetBinContent(bin);
    f1->SetParameter(1,yp);
    // I fix the sigma
    f1->SetParameter(2,sigma);

    spectrum->Fit(f1,"RIQE");
    vec_pos.push_back(xp);
    vec_sigma.push_back(f1->GetParameter(2));

    cout << "FITTING PEAK AT " << xpeak.at(0) << endl;
    xp = xpeak.at(0);
    f1->SetParameter(0,xp);
    //f1->SetParLimits(1,xp-sigma,xp+sigma);
    // Then I get amplitude
    bin = spectrum->GetXaxis()->FindBin(xp);
    yp = spectrum->GetBinContent(bin);
    f1->SetParameter(1,yp);
    // I fix the sigma
    f1->SetParameter(2,0.02);

    spectrum->Fit(f1,"RIQE");
    //cout << "We fitted : " << xp << " with sigma " <<  f1->GetParameter(2) << endl;
    vec_pos.push_back(xp);
    vec_sigma.push_back(f1->GetParameter(2));
  }
  else
  {
    vec_pos.push_back(0);
    vec_sigma.push_back(0);
    vec_pos.push_back(0);
    vec_sigma.push_back(0);
  }
  return std::make_tuple(vec_pos,vec_sigma);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t PSDMatrixAnalyzer(TH2F *matrix)
{
  // std:vector<vector<Double_t>> matrice;
  //
  // // I Initialize the matrix:

  // for(int binx=0; binx < nbinx; binx++)
  // {
  //   std::vector<Double_t> localvector;
  //   for(int biny=0; biny < nbiny; biny++) localvector.push_back(matrix->GetBinContent(binx,biny));
  //   matrice.push_back(localvector);
  // }
  TSpectrum2 *s;
  Int_t nbinsx = matrix->GetXaxis()->GetNbins();
  Int_t nbinsy = matrix->GetYaxis()->GetNbins();
  Double_t xmin   = matrix->GetXaxis()->GetXmin();
  Double_t xmax   = matrix->GetXaxis()->GetXmax();
  Double_t ymin   = matrix->GetYaxis()->GetXmin();
  Double_t ymax   = matrix->GetYaxis()->GetXmax();
  Double_t dx = (xmax-xmin)/nbinsx;
  Double_t dy = (ymax-ymin)/nbinsy;

  cout << "Searching peaks ..." << endl;
  //now the real stuff: Finding the peaks
  Int_t nfound = s->Search(matrix,100,"");

  cout << "Found " << nfound << " peaks" << endl;
  //
  // TF2 *f2 = new TF2("f2",fpeaks2,xmin,xmax,ymin,ymax,5*nfound);
  // f2->SetNpx(100);
  // f2->SetNpy(100);
  // f2->SetParameters(par);

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int EnergyCalibrator(const CExperiment &experiment)
{
  //ROOT::DisableImplicitMT();
  // I will work in multithread to replace each file by an energy calibrated file
  // Initialisation of multithreads
  TThread::Initialize();

  using namespace indicators;
  int nOfThread = experiment.GetThreadsNbr();
  int chainnbr = experiment.GettheTChain().size();
  std::vector<TChain *> tab_chained_oak = experiment.GettheTChain();

  // I create one progress bar per thread
  if(chainnbr < nOfThread) nOfThread = chainnbr;
  cout << "Planned number of cores " << nOfThread << "\t number of files: " << chainnbr << endl;
  ProgressBar *bars = new ProgressBar[chainnbr];

  DynamicProgress<ProgressBar> Dbars;
  Dbars.set_option(option::HideBarWhenComplete{true});
  for (int i = 0; i < chainnbr; i++)
  {
    TChain *chained_oak = tab_chained_oak.at(i);
    string localfilename = chained_oak->GetName();
    cout << "Creating a progress bar for " << localfilename << endl;
    int lastposition = localfilename.find_last_of("\/");
    TString Posttext = "Calibrating "; Posttext += localfilename.substr(lastposition+1);
    bars[i].set_option(option::BarWidth{5});
    bars[i].set_option(option::Start{"["});
    bars[i].set_option(option::Fill{"="});
    bars[i].set_option(option::Lead{">"});
    bars[i].set_option(option::Remainder{"-"});
    bars[i].set_option(option::End{"]"});
    bars[i].set_option(option::PrefixText{Posttext.Data()});
    bars[i].set_option(option::ShowElapsedTime{true});
    bars[i].set_option(option::ShowRemainingTime{true});
    bars[i].set_option(option::ForegroundColor{Color::grey});
    bars[i].set_option(option::FontStyles{std::vector<FontStyle>{FontStyle::bold}});
  }


  int chainbr[1] = {0};
  mutex myMutex;
  mutex localMutex;


/*   Progress Bar should be declared in the "calibrate function" and added to Dbars within the function*/
/* if still not working then fuck off */


  //*********************Lambda function that will be launched by all threads**************************//
  auto Processfunction = [&]()
  {
    while(chainbr[0] < chainnbr)
    {
      chainbr[0]++;
      int localchainnbr = chainbr[0]-1;
      TChain *chained_oak = tab_chained_oak.at(localchainnbr);
      string localfilename = chained_oak->GetName();
      int lastposition = localfilename.find_last_of("\/");
      TString Posttext = "Calibrating "; Posttext += localfilename.substr(lastposition+1);

      // Getting the progressbar
      Dbars.push_back(bars[localchainnbr]);
      int check = NRJCalibrator(experiment,chained_oak,Dbars,localchainnbr,&localMutex);
      TString LastPosttext = "Completed : "; LastPosttext += localfilename.substr(lastposition+1);
      if(check == 1) {
        Dbars[localchainnbr].set_option(option::PostfixText{LastPosttext.Data()});
        Dbars[localchainnbr].mark_as_completed();
        break;
          //bars[localchainnbr].mark_as_completed();
        }
      //bar.set_progress((int)((double)chainbr[0]/(double)chainnbr*100.));
      //cout << "FUUUUUUUUUUCCCCKKK" << endl;
      if(check == 0)
      {
        cout << SetBOLD << SetForeRED << endl;
        cout << " Problem to apply the Calibration properly to the TChain #" << localchainnbr << endl;
        cout << " Check data file(s) " << endl,
        cout << RESETTEXT << endl;
        break;
        //return 0;
      }
    }

  };
  //---------------------------------------------------------------------------------------------------//
  vector<thread> threads;
  for(int i = 0; i < chainnbr; i++)
  {
    threads.emplace_back(Processfunction);
  }

  for(auto &&t : threads)
  {
    t.join(); //sleep(15);
  }

  cout << RESETTEXT << "All Done" << endl;

  ROOT::EnableImplicitMT();

  return 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int NRJCalibrator(const CExperiment &experiment,TChain *chained_oak, DynamicProgress<ProgressBar> &BAR, int localchainnbr, mutex *stop)
{

  //cout << "I DO HAVE QDC2" << experiment.GetisQDC2() << endl;
  //ProgressBar DBar = BAR[BARnumber];

  //BAR[BARnumber].tick();
  //myMutex->lock();
  //chained_oak->Print();
  //cout << endl << "Working on TChain based on file "<< chained_oak->GetName() << endl;
  TString inputfilename = chained_oak->GetName();
  //auto maple = chained_oak->GetTree();
  //myMutex->unlock();
  int it1 = inputfilename.Index(".root",5,1,inputfilename.kExact);
  TString outputfilename = inputfilename(0,it1); outputfilename += "_ECal.root";

  // Declaration of new tree variables
  nrj_type nrj = 0;
  nrj_type nrj2 = 0;
  nrj_type nrj3 = 0;
  nrj_type nrj4 = 0;
  label_type index =0;
  Bool_t pileup = false;
  tm_type tm = 0;

  // Declaration of old tree variables
  nrj_Rawtype  enrj = 0;
  nrj_Rawtype enrj2 = 0;
  nrj_Rawtype enrj3 = 0;
  nrj_Rawtype enrj4 = 0;
  label_type index1 = 0;
  Bool_t pileup1 = 0;
  tm_Rawtype tm1 = 0;

  // maple->SetBranchAddress ("label", &index1);
  // maple->SetBranchAddress ("nrj", &enrj);
  // if(experiment.GetisQDC2()) maple->SetBranchAddress ("nrj2", &enrj2);
  // if(experiment.GetisQDC3()) maple->SetBranchAddress ("nrj3", &enrj3);
  // if(experiment.GetisQDC4()) maple->SetBranchAddress ("nrj4", &enrj4);
  // maple->SetBranchAddress ("time", &tm1);
  // maple->SetBranchAddress ("pileup",&pileup1);

  TFile *outputfile = new TFile(outputfilename,"RECREATE");
  TTree *oak = new TTree("DataTree",outputfilename);

  // Declaration of Branches that will contain the data
  oak->Branch ("label", &index, "label/s");
  oak->Branch ("nrj", &nrj, "nrj/D");
  if(experiment.GetisQDC2()) oak->Branch ("nrj2", &nrj2,"nrj2/D");
  if(experiment.GetisQDC3()) oak->Branch ("nrj3", &nrj3,"nrj3/D");
  if(experiment.GetisQDC4()) oak->Branch ("nrj4", &nrj4,"nrj4/D");
  oak->Branch ("time", &tm, "time/l");
  oak->Branch ("pileup",&pileup,"pileup/O");

  //cout << "New Seed burried " << std::endl;
  //cout << "and in the memory another tree will grow..."<<std::endl ;

  const auto chainentries = chained_oak ->GetEntries();
  int localchainentries = (int) (chainentries);
  int percentage = 10;
  int percent = (int)(percentage/100.*chainentries);
  //cout << "Prepare to recycle the " << chainentries << " entries..." << std::endl;
  BAR[localchainnbr].set_option(option::MaxProgress{localchainentries});

  //mutex stop;
  ROOT::TTreeProcessorMT TP(inputfilename.Data(),experiment.GetDataTreeName());
  int threadnbr = 0;
  auto loop_and_fill = [&] (TTreeReader &myReader){
    TTreeReaderValue<label_Rawtype> labelRV(myReader,"label");
    TTreeReaderValue<nrj_Rawtype> QDC1RV(myReader,"nrj");
    TTreeReaderValue<nrj_Rawtype> QDC2RV(myReader,"nrj2");
    TTreeReaderValue<tm_Rawtype> TMRV(myReader,"time");
    TTreeReaderValue<pu_type> PURV(myReader,"pileup");
    
    ULong64_t hitnumber = 0;
    CHit *hit = new CHit(threadnbr++);
    while(myReader.Next())
      {
        hitnumber++;threadnbr++;
        if(hitnumber%percent==0) // Printing every X percent
        {
          //stop->lock();
          BAR[localchainnbr].set_option(option::PostfixText{
            std::to_string((int)hitnumber) + "/" + std::to_string((int)localchainentries)
          });
          BAR[localchainnbr].set_progress(percentage);//printProgress((double)hitI/(double)chainentries);
          //BAR[localchainnbr].tick();
          //stop->unlock();
        }
        //CHit *hit = new CHit(hitnumber+1000*threadnbr);
        hit->SetHit(*labelRV, *TMRV, *QDC1RV, *QDC2RV, *PURV);
        //hit->PrintHit();
        //cout << "hitnumber " << hitnumber << endl;
        index = *labelRV;
        tm = (Double_t)*TMRV;
        pileup = *PURV;
        
        int detindex = experiment.GetLabel2Detnbrs(static_cast<int>(index));
        
        if(index < 200) // HPGe && BGO
        {
          //cout << "Calibrating " << experiment.GetDetector(detindex)->GetDetectorName() << endl;
          nrj = experiment.GetDetector(detindex)->GetEnergy(*QDC1RV);
          nrj2 = 0;
          //cout << endl << endl << index <<  "\t enrj = " << enrj << "\t" <<"NRJ = " << NRJ << "\t" << "nrj = " << nrj << endl;
          if(nrj > 10. && experiment.GetDetector(detindex)->GetCalibx1()!=1.) {stop->lock();oak->Fill();stop->unlock();}
        }
        
        // Case of PARIS
        if(index1 > 0)
        {
          auto NRJ    = *QDC1RV;  // Short Gate
          auto NRJ2   = *QDC2RV;  // Long Gate
          // Now I treat the energies depending of what part of PARIS has been fired
          //CHit *hit = new CHit(hitnumber);
          //hit->SetHit(index, 0, NRJ, NRJ2, false);
          
          Double_t PSD = hit->PerformPSD();// PSD = (Long-Short)/Long
          Bool_t isLaBr = experiment.GetDetector(detindex)->IsPureLaBr3(PSD,hit->GetHitE1(),hit->GetHitE2());
          Bool_t isBeyondLaBrandNaI = experiment.GetDetector(detindex)->IsBeyondLaBr3andNaI(PSD,hit->GetHitE1(),hit->GetHitE2());


          if(!isBeyondLaBrandNaI)
          {
            if(isLaBr)
            {
              nrj = experiment.GetDetector(detindex)->GetEnergy(NRJ);
              nrj2 = 0.;
              //if(nrj > 10. && experiment.GetDetector(detindex)->GetCalibx1()!=1. && experiment.GetDetector(detindex)->GetCalibx0()!=0.) oak->Fill();
            }
            else if(!isLaBr)
            {
              tie(NRJ2,NRJ) = Rotation((double)experiment.GetDetector(detindex)->GetRotationAngle(),(double)experiment.GetDetector(detindex)->GetRotationAngleTan(),NRJ2,NRJ);
              nrj = experiment.GetDetector(detindex)->GetEnergy(NRJ);
              nrj2 = experiment.GetDetector(detindex)->GetSecondaryEnergy(NRJ2);
              //nrj2 = experiment.GetDetector(detindex)->GetSecondaryEnergy(NRJ2);
              //if(nrj > 10. && experiment.GetDetector(detindex)->GetCalibx1()!=1. && experiment.GetDetector(detindex)->GetCalibx0()!=0.) oak->Fill();
            }
          }
          //cout << endl << endl << index <<   "enrj = " << enrj << "\t" <<"NRJ = " << NRJ << "\t" << "nrj = " << nrj << endl;

          // I place a 50 keV threshold && that the detector has been calibrated.
          
        }
        
        hit->Clear();
      }
      delete hit;
  };
  TP.Process(loop_and_fill);
  cout << "File Converted Saving" << endl;
  outputfile->cd();
  oak->Write("", TObject::kOverwrite);
  //oak->Delete();
  outputfile->Close();
  BAR[localchainnbr].set_progress(100);
  BAR[localchainnbr].mark_as_completed();

  outputfile->Delete();

  //cout << "File " << inputfilename << " energy calibrated "<< endl;
  chained_oak->GetFile()->Close();

  return 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<TH1F*>  DrawTimeShifts_ECal(const CExperiment &experiment, Double_t deltaTinit, Double_t deltaTfin)
{
  //Definition du germe pour le tirage aleatoire
  srand48(time(NULL));

  //Declaration of all the variable used in the function
  ULong64_t chainentries;
  TString detname1;
  TString detname2;
  Double_t deltaT(0);
  //Double_t reso;
  int nperim;
  Int_t nbrchannels = (deltaTfin-deltaTinit)*20; //precision of 50 ps

  //Declaration of all variables for names
  TString title;
  TString spectrumname;

  // Declaration of time spectra
  int nbrofspectra = experiment.GetDetectors().size();
  std::vector<TH1F*> timespectra;
  TH2F *timematrix;
  int highestdetlabel(0);

  for(int sindex = 0; sindex < nbrofspectra; sindex++)
  {
    TH1F* localtimespectrum;
    title = "Time Spectrum of detector ";
    spectrumname = "timespectrum";
    title += experiment.GetDetector(sindex)->GetDetectorName();
    spectrumname +=experiment.GetDetector(sindex)->GetDetectorName();
    if(experiment.GetDetector(sindex)->GetDetectorlabel()>highestdetlabel) highestdetlabel = experiment.GetDetectors().at(sindex)->GetDetectorlabel();
    title += " vs ";
    spectrumname +="vs";
    title += experiment.GetReferenceDetector()->GetDetectorName();
    spectrumname += experiment.GetReferenceDetector()->GetDetectorName();
    localtimespectrum = new TH1F(spectrumname,title,nbrchannels,deltaTinit,deltaTfin);
    timespectra.push_back(localtimespectrum);
    spectrumname.Clear();
    title.Clear();
  }


  // I declare la Time matrix to check alignement later
  timematrix = new TH2F("timealignementmatrix","Time spectra of all detectors",highestdetlabel,1,highestdetlabel,nbrchannels,deltaTinit,deltaTfin);

  // Creation of the file to save all the data
  TString outputfilename = experiment.GetFileDirectory_OUT();
  outputfilename += "Timespectra_all.root";
  TFile *outputfile = new TFile(outputfilename,"RECREATE");

  // Loading the TTree for reading the Data
  std::vector<TChain *> tab_chained_oak = experiment.GettheTChain();
  TChain *chained_oak = tab_chained_oak.at(0);
  label_type index;
  tm_Rawtype tm;
  nrj_type enrj,enrj2,enrj3,enrj4;
  std::cout << FOREBLU << "Loading Tree ..." << std::endl;
  chained_oak -> SetBranchAddress("label",&index);   // Detector number
  chained_oak -> SetBranchAddress("time",&tm);       // Time of the hit
  chained_oak -> SetBranchAddress("nrj",&enrj);
  if(experiment.GetisQDC2()) chained_oak->SetBranchAddress ("nrj2", &enrj2);
  if(experiment.GetisQDC3()) chained_oak->SetBranchAddress ("nrj3", &enrj3);
  if(experiment.GetisQDC4()) chained_oak->SetBranchAddress ("nrj4", &enrj4);
  //chained_oak -> SetBranchAddress("pileup",&pileup); // Time of the hit
  std::cout << "Loaded .." << RESETTEXT << std::endl;

  // Another Timer
  TStopwatch timer2;

  //Starting of the chronometer
  timer2.Reset();
  timer2.Start();

  // Energy windows
  nrj_type E_ref_min = 1300;
  nrj_type E_ref_max = 1360;
  nrj_type E_det_min = 1140;
  nrj_type E_det_max = 1200;

  //---------------------------------------------------------------------------//
  //                                                                           //
  //                    Coincidence reconstruction Algorithm                   //
  //                                                                           //
  //---------------------------------------------------------------------------//
  // Determination of the number of event that have to be treated
  chainentries = chained_oak -> GetEntries();
  int percent = (int)(0.05*chainentries);

  std::cout << "Beginning of coincidences research" << std::endl;
  std::cout << "On " << chainentries << " entries" << std::endl;

  // Creating a ProgressBar which is nice in MT mode
  using namespace indicators;
  ProgressBar bar{
    option::BarWidth{60},
    option::Start{"["},
    option::Fill{"="},
    option::Lead{">"},
    option::Remainder{"-"},
    option::End{"]"},
    option::PostfixText{"Reading"},
    option::ShowElapsedTime{true},
    option::ShowRemainingTime{true},
    option::ForegroundColor{Color::grey},
    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
  };

  tm_type memorytm = 0;

  // I declare a hit collection that will be used for event reconstruction
  CHitCollection *coinc_windows = new CHitCollection();
  coinc_windows->SetCollectionTimeSize(deltaTfin-deltaTinit);
  int ReferenceLabel = experiment.GetReferenceDetector()->GetDetectorlabel();

  for(ULong64_t hitI=0; hitI < chainentries; hitI++)
  {
    if(hitI%percent==0) // Printing every 5 percent
    {
      bar.set_progress((int)((double)hitI/(double)chainentries*100.));
    }

    // Loading of the next entry to compare with the previous one
    int hitexist = chained_oak -> GetEntry(hitI);

    if(hitexist > 0)
    {
      CHit *hit = new CHit(hitI);
      if(!experiment.GetisQDC2()) hit->SetHit(index, tm, enrj, 1);
      if(experiment.GetisQDC2())  hit->SetHit(index, tm, enrj, enrj2, 1);
      if(experiment.GetisQDC3())  hit->SetHit(index, tm, enrj, enrj2, enrj3, 1);
      if(experiment.GetisQDC4())  hit->SetHit(index, tm, enrj, enrj2, enrj4, 1);

      if(hitI > 450 && hitI < 460) hit->PrintHit();

      // If the new hit is in the right time window I add the hit to the collection
      if(coinc_windows->IsHitInside(hit)) {coinc_windows->AddHit(hit);}


      else
      {
        // New hit out of time window
        // I search for the reference detector in the time window
        int Refdetector_pos = coinc_windows->IsReferenceDetectorIn(ReferenceLabel);

        if(Refdetector_pos >= 0 && coinc_windows->GetCollectionSize() > 1)
        {
          cout << endl << endl;
          coinc_windows->PrintHitCollection();
        }

        // If just one or that I have finished working with it I clear the collection and add the last hit that was not added
        coinc_windows->Clear();
        coinc_windows->AddHit(hit);
      }

      delete hit;
    }
  }

  delete coinc_windows;

  std::cout << endl << "Saving in " << outputfilename << std::endl;
  outputfile->cd();
  for (int i=0; i< (int)timespectra.size();i++)
  {
    timespectra.at(i)->Write();
  }
  timematrix->Write();
  outputfile->Close();

  cout << "Time spectra are saved " << endl;


  // Printing of chronometer measurement
  timer2.Stop();
  Double_t rtime2 = timer2.RealTime();
  Double_t ctime2 = timer2.CpuTime();
  std::cout << std::endl;
  std::cout << "End of Time Shifts Calculations" << std::endl;
  std::cout << "# RealTime=" << rtime2 << " seconds, CpuTime="<< ctime2 << " seconds" <<std::endl;
  std::cout << std::endl;

  return timespectra;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<TH1F*>  DrawTimeShifts(const CExperiment &experiment, Double_t deltaTinit, Double_t deltaTfin,  TString filename)
{
  //Definition du germe pour le tirage aleatoire
  srand48(time(NULL));

  //Declaration of all the variable used in the function
  ULong64_t chainentries;
  TString detname1;
  TString detname2;
  Double_t deltaT(0);
  //Double_t reso;
  int nperim;
  Int_t nbrchannels = (deltaTfin-deltaTinit)*20; //precision of 50 ps

  //Declaration of all variables for names
  TString title;
  TString spectrumname;

  // Declaration of time spectra
  int nbrofspectra = experiment.GetDetectors().size();
  std::vector<TH1F*> timespectra;
  TH2F *timematrix;
  int highestdetlabel(0);

  for(int sindex = 0; sindex < nbrofspectra; sindex++)
  {
    TH1F* localtimespectrum;
    title = "Time Spectrum of detector ";
    spectrumname = "timespectrum";
    title += experiment.GetDetector(sindex)->GetDetectorName();
    spectrumname +=experiment.GetDetector(sindex)->GetDetectorName();
    if(experiment.GetDetector(sindex)->GetDetectorlabel()>highestdetlabel) highestdetlabel = experiment.GetDetectors().at(sindex)->GetDetectorlabel();
    title += " vs ";
    spectrumname +="vs";
    title += experiment.GetReferenceDetector()->GetDetectorName();
    spectrumname += experiment.GetReferenceDetector()->GetDetectorName();
    localtimespectrum = new TH1F(spectrumname,title,nbrchannels,deltaTinit,deltaTfin);
    timespectra.push_back(localtimespectrum);
    spectrumname.Clear();
    title.Clear();
  }


  // I declare la Time matrix to check alignement later
  timematrix = new TH2F("timealignementmatrix","Time spectra of all detectors",highestdetlabel,1,highestdetlabel,nbrchannels,deltaTinit,deltaTfin);

  // Creation of the file to save all the data
  TString outputfilename = experiment.GetFileDirectory_OUT();
  outputfilename += "Timespectra_all.root";
  TFile *outputfile = new TFile(outputfilename,"RECREATE");

  // Loading the TTree for reading the Data
  std::vector<TChain *> tab_chained_oak = experiment.GettheTChain();
  TChain *chained_oak = tab_chained_oak.at(0);
  label_type index;
  tm_Rawtype tm;
  nrj_type enrj,enrj2,enrj3,enrj4;
  std::cout << FOREBLU << "Loading Tree ..." << std::endl;
  chained_oak -> SetBranchAddress("label",&index);   // Detector number
  chained_oak -> SetBranchAddress("time",&tm);       // Time of the hit

  //chained_oak -> SetBranchAddress("pileup",&pileup); // Time of the hit
  std::cout << "Loaded .." << RESETTEXT << std::endl;

  // Another Timer
  TStopwatch timer2;

  //Starting of the chronometer
  timer2.Reset();
  timer2.Start();

  //---------------------------------------------------------------------------//
  //                                                                           //
  //                    Coincidence reconstruction Algorithm                   //
  //                                                                           //
  //---------------------------------------------------------------------------//
  // Determination of the number of event that have to be treated
  chainentries = chained_oak -> GetEntries();
  int percent = (int)(0.05*chainentries);

  std::cout << "Beginning of coincidences research" << std::endl;
  std::cout << "On " << chainentries << " entries" << std::endl;

  // Creating a ProgressBar which is nice in MT mode
  using namespace indicators;
  ProgressBar bar{
    option::BarWidth{60},
    option::Start{"["},
    option::Fill{"="},
    option::Lead{">"},
    option::Remainder{"-"},
    option::End{"]"},
    option::PostfixText{"Reading"},
    option::ShowElapsedTime{true},
    option::ShowRemainingTime{true},
    option::ForegroundColor{Color::grey},
    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
  };

  tm_type memorytm = 0;

  // I declare a hit collection that will be used for event reconstruction
  CHitCollection *coinc_windows = new CHitCollection();
  coinc_windows->SetCollectionTimeSize(deltaTfin-deltaTinit);
  int ReferenceLabel = experiment.GetReferenceDetector()->GetDetectorlabel();

  for(ULong64_t hitI=0; hitI < chainentries; hitI++)
  {
    if(hitI%percent==0) // Printing every 5 percent
    {
      bar.set_progress((int)((double)hitI/(double)chainentries*100.));
    }

    // Loading of the next entry to compare with the previous one
    int hitexist = chained_oak -> GetEntry(hitI);

    if(hitexist > 0)
    {
      CHit *hit = new CHit(hitI);
      hit->SetHit(index, tm, 0, 1);

      //if(hitI > 450 && hitI < 460) hit->PrintHit();

      // If the new hit is in the right time window I add the hit to the collection
      if(coinc_windows->IsHitInside(hit)) {coinc_windows->AddHit(hit);}


      else
      {
        // New hit out of time window
        // I search for the reference detector in the time window
        int Refdetector_pos = coinc_windows->IsReferenceDetectorIn(ReferenceLabel);

        if(Refdetector_pos >= 0 && coinc_windows->GetCollectionSize() > 1)
        {
          //cout << endl << endl;
          //coinc_windows->PrintHitCollection();
          // Now I loop on the hit collection to get all the time spectra
          for(int HitC = 0; HitC < coinc_windows->GetCollectionSize(); HitC++)
          {
            if(HitC != Refdetector_pos)
            {
              // I calculate the time difference to Fill the time spectrum with
                deltaT = (double)((Long64_t)coinc_windows->GetHit(HitC).GetHitTime() -(Long64_t)coinc_windows->GetHit(Refdetector_pos).GetHitTime())/1000.;

                // I fill the correct time spectrum
                int spectrumindex = experiment.GetLabel2Detnbrs(coinc_windows->GetHit(HitC).GetHitLabel());
                //cout << coinc_windows->GetHit(HitC).GetHitLabel() << "\t" << spectrumindex << "\t" << deltaT << endl;
                timespectra.at(spectrumindex)->Fill(deltaT);

                // And the time matrix
                timematrix->Fill(coinc_windows->GetHit(HitC).GetHitLabel(),deltaT);
            }
          }

        }

        // If just one or that I have finished working with it I clear the collection and add the last hit that was not added
        coinc_windows->Clear();
        coinc_windows->AddHit(hit);
      }

      delete hit;
    }
  }

  delete coinc_windows;

  std::cout << endl << "Saving in " << outputfilename << std::endl;
  outputfile->cd();
  for (int i=0; i< (int)timespectra.size();i++)
  {
    timespectra.at(i)->Write();
  }
  timematrix->Write();
  outputfile->Close();

  cout << "Time spectra are saved " << endl;


  // Printing of chronometer measurement
  timer2.Stop();
  Double_t rtime2 = timer2.RealTime();
  Double_t ctime2 = timer2.CpuTime();
  std::cout << std::endl;
  std::cout << "End of Time Shifts Calculations" << std::endl;
  std::cout << "# RealTime=" << rtime2 << " seconds, CpuTime="<< ctime2 << " seconds" <<std::endl;
  std::cout << std::endl;

  return timespectra;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int CalculateTimealignementShifts(const CExperiment &experiment, Bool_t isCalibrated)
{
  //Definition du germe pour le tirage aleatoire
  srand48(time(NULL));

  std::vector<TH1F*> TIMESpectra;

  // Define all the TH1F for the analysis
  TString outputfilename = experiment.GetFileDirectory_OUT();
  outputfilename += "Timespectra_all.root";
  TFile *spectrafile = new TFile(outputfilename,"READ");
  if(!spectrafile || spectrafile->IsZombie() || spectrafile->GetNkeys()==0){
    cout << "File" << outputfilename << " with time spectra do not exists; creating it ..." << endl;
    Double_t maxtimewindow = 600.;
    Double_t mintimewindow = -600.;
    cout << "A time window of ["<< mintimewindow << ";" << maxtimewindow << "] ns will be used" << endl;
    if(!isCalibrated) TIMESpectra = DrawTimeShifts(experiment,mintimewindow,maxtimewindow,outputfilename); // Time window in ns
    else TIMESpectra = DrawTimeShifts_ECal(experiment,mintimewindow,maxtimewindow);
    cout << "Number of spectra to analyse = " << TIMESpectra.size() << endl;
  }
  //NRJpectra = DrawAllEnergyUncalibratedSpectra(experiment);
  else
  {
    cout << FORECYN << "File already exists" << endl;
    cout << "Loading spectra from file..." << endl;
    for(auto sindex = 0; sindex < (int)experiment.GetDetectors().size();sindex++)
    {
      TString spectrumname = "timespectrum";
      spectrumname+=experiment.GetDetectors().at(sindex)->GetDetectorName();
      spectrumname+="vs";spectrumname+=experiment.GetReferenceDetector()->GetDetectorName();
      TH1F *temp_h1 = (TH1F*)spectrafile->Get(spectrumname);
      TIMESpectra.push_back(temp_h1);
    }
    cout << "Spectra Loaded ..." << RESETTEXT << endl;
    cout << "Number of spectra to analyse = " << TIMESpectra.size() << endl;
  }

  if(TIMESpectra.size() == 0)
  {
    cout << FORERED << SetBOLD << endl;
    cout << "Problem with Uncalibrated Energy spectra generation" << endl;
    cout << RESETTEXT << endl;
    return 0;
  }

  // I declare two files.
  // One to store the TGraph with the linear fit..
  std::vector<TChain *> tab_chained_oak = experiment.GettheTChain();
  TChain *chained_oak = tab_chained_oak.at(0);
  TString temp = chained_oak->GetName();
  int it1 = temp.Index("Calibration_",12,1,temp.kExact);temp = temp(it1,temp.Length());
  it1 = temp.Index(".",1,temp.kExact);temp = temp(0,it1);
  TString calibfilename=experiment.GetFileDirectory_OUT();
  //calibfilename += temp;//experiment.GetNRJCalibration_filename();
  calibfilename +="deltaT.dat";
  cout << "The Calibration parameters will be saved in " << calibfilename << RESETTEXT<< endl;
  ofstream calibfile(calibfilename,ios::out);
  calibfile.precision(4);

  // Now I proceed to the Analysis of each spectrum
  int bouik = TIMESpectra.size();

  // Getting coincidence peak position & time resolutions for all detectors
  std::vector<Double_t> poz_and_rez;
  for(int sindex = 0; sindex < (int)bouik; sindex++)
  {
    // Checking on LaBr3 who are going to have best time resolution
    Bool_t isLaBr = kFALSE;
    cout << endl << "Analysing spectrum #" << sindex+1 << "/" << bouik << " named " <<  TIMESpectra.at(sindex)->GetName() << endl;
    if(experiment.GetDetectors().at(sindex)->GetDetectorType() == "LaBr") isLaBr = kTRUE;
    cout << isLaBr << " " << TIMESpectra.at(sindex)->GetXaxis()->GetNbins()<< endl;

    poz_and_rez = DeltaTmeasurer(TIMESpectra.at(sindex), isLaBr);
    cout << "Time shift = " << poz_and_rez.at(0) << " & Time Resolution = " << poz_and_rez.at(1) << endl;
    calibfile << experiment.GetDetectors().at(sindex)->GetDetectorlabel() << "\t" << poz_and_rez.at(0) << "\t" << poz_and_rez.at(1) << endl;
  }
  cout << "Time Calibration parameters have been written to " << calibfilename << endl;
  calibfile.close();

  return 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<Double_t> DeltaTmeasurer(TH1F *timespectrum, bool isqdc)
{
  // Preparing the output vector
  std::vector<Double_t> pos_and_rez;

  // I prepare a peak width in advance
  Double_t peakwidth = 20; // in ns
  if(isqdc) peakwidth = 3; // in ns

  // Useful variables
  int i=0;

  // I got the spectrum, I define the TSpectrum
  TSpectrum *s = new TSpectrum();

  // I try to determinate the background with compton edges
  Int_t nbins = timespectrum->GetXaxis()->GetNbins();
  Double_t xmin  = 0;
  Double_t xmax  = nbins;
  Double_t * source = new Double_t[nbins];
  for (int i = 0; i < nbins; i++) source[i]=timespectrum->GetBinContent(i + 1);

  // Now I search for the background
  // I prepare a TH1 to store the background
  // Now I search for the background
  s->Background(source,nbins,10,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder2,kTRUE,TSpectrum::kBackSmoothing15,kTRUE);

  // I prepare a TH1 to store the background
  TH1F *d1 = new TH1F("d1","",nbins,xmin,xmax);
  // I fill the spectrum witht the background
  for (i = 0; i < nbins; i++) d1->SetBinContent(i + 1,source[i]);

  // Now I search for peaks
  //for (i = 0; i < nbins; i++) timespectrum->SetBinContent(i + 1,timespectrum->GetBinContent(i+1)-source[i]);
  //timespectrum->Add(timespectrum,d1,1,-1);
  Int_t nfound = s->Search(timespectrum,peakwidth,"",0.5);
  //nfound = s->Search(timespectrum,peakwidth,"",0.5);
  std::cout <<"Found " << nfound << " candidate peaks to fit" << std::endl;
  nfound = 1;


  // I define the data for the peak
  Double_t *pos      = new Double_t[nfound];
  Double_t *pos_err  = new Double_t[nfound];
  Double_t *reso = new Double_t[nfound];
  Double_t *integral = new Double_t[nfound];
  Double_t *FWTM = new Double_t[nfound];

  // I loop on the peaks
  // I load from TSpectrum the peak list
  Double_t *xpeaks = s->GetPositionX();
  Int_t npeaks = 0;
  Double_t sigma=0.;
  if (nfound < 3)
  {
    for (int p=0;p<nfound;p++)//nfound;p++)
    {
      // First I get bin center
      Double_t xp = xpeaks[p];

      std::cout << "xp " << xp << std::endl;

      // Define my gaussian around my peak
      TF1 *f1 = new TF1("f1", "gaus", xp-20, xp+20);
      f1->SetParLimits(1,xp-5,xp+5);

      // Then I get amplitude
      //Int_t bin = timespectrum->GetXaxis()->FindBin(xp);
      //Double_t yp = timespectrum->GetBinContent(bin);
      //f1->SetParLimits(0,yp-TMath::Sqrt(yp),yp+TMath::Sqrt(yp));

      // I fix the sigma
      f1->SetParLimits(3,10,30);

      // Then I fit my spectrum with a gaussian
      timespectrum->Fit("f1","RIQE");

      // I store and print info
      //std::cout << "Peak amplitude " << f1->GetParameter(0) << std::endl;
      //std::cout << "Peak mean value " << f1->GetParameter(1) << "p/m " << f1->GetParError(1) << std::endl;
      pos[p] = f1->GetParameter(1);
      //std::cout << "Peak sigma " << f1->GetParameter(2) << std::endl;
      sigma=f1->GetParameter(2);
      reso[p] = 2.35482*sigma;
      FWTM[p] = 4.29193*sigma;
      integral[p] = timespectrum->Integral(TMath::Floor(pos[p]-5*sigma),TMath::Floor(pos[p]+5*sigma));
      pos_err[p] =sigma/TMath::Sqrt(integral[p]);

      npeaks++;
    }
  }
  Double_t resof = reso[0];

  if(isqdc) pos_and_rez.push_back(xpeaks[0]);//pos[0];
  else pos_and_rez.push_back(pos[0]);

  //std::cout << "DeltaT sent back " << pos_and_rez.at(0) << std::endl;

  delete[] source;
  delete[] pos;
  delete[] pos_err;
  delete[] reso;
  delete[] integral;
  delete[] FWTM;

  pos_and_rez.push_back(resof);

  // Sending back the vector that contains:
  // 1 - Coinc peak positions
  // 2 - Detector Resolutions
  return pos_and_rez;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int CheckCoincidenceWindow(std::vector<tm_Rawtype> tab_memory_tm, Double_t deltaTfin)
{

  int jmax = tab_memory_tm.size()-1;
  int nperim(0);
  Double_t deltaT(0.);


  for(int ll = 0; ll < jmax;ll++)
  {
    // Calculation of the time difference between detectors
    deltaT = (Double_t)((Long64_t)tab_memory_tm.at(jmax) -(Long64_t)tab_memory_tm.at(ll))/1000.;
    //deltaT = (Double_t)(tab_memory_tm.at(jmax)-tab_memory_tm.at(ll)); // to get in ns
    //deltaT = deltaT/1000.;
    //cout << deltaT << endl;

    // I count the number of hit which are not in the time window anymore
    if((deltaT > deltaTfin))
    {
      //Coucou(1);
      nperim++;
    }
  }
  return nperim;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
