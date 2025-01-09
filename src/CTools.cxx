#include "CTools.h"

// Definition of the I/O function
using namespace std;
using namespace ROOT;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress (double percentage)
{
  float val = (float) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%.3f%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t VectorMean(std::vector<Double_t>& theVector)
{
  double sum = std::accumulate(theVector.begin(), theVector.end(), 0.0);
  double mean = sum / theVector.size();

  return mean;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t VectorWeightedMean(std::vector<Double_t>& theVector, std::vector<Double_t>& theWeight)
{
  double sum(0.);
  double sumweight(0);
  for(auto i = 0; i < (int)theVector.size(); i++)
  {
    sum += theVector.at(i)*theWeight.at(i);
    sumweight += theWeight.at(i);
  }
  double mean = sum / sumweight;

  return mean;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t VectorStdDev(std::vector<Double_t>& theVector)
{
  double mean = VectorMean(theVector);
  std::vector<double> diff(theVector.size());
  std::transform(theVector.begin(), theVector.end(), diff.begin(), [mean](double x) { return x - mean; });
  double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  double stdev = std::sqrt(sq_sum / theVector.size());

  return stdev;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int SimpleRoundTo(Double_t num){return (int)(num < 0 ? (num - 0.5) : (num + 0.5));};
int SimpleRoundTo(Float_t num){return (int)(num < 0 ? (num - 0.5) : (num + 0.5));};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::tuple<double, double> Rotation(Double_t theta, Double_t theta_tan, Double_t x, Double_t y)
{
  Double_t xf = x*TMath::Cos(theta)-y*TMath::Sin(theta)*(TMath::Abs(theta_tan)/theta_tan);
  Double_t yf = x*TMath::Sin(theta)*(TMath::Abs(theta_tan)/theta_tan)+y*TMath::Cos(theta);

  return std::make_tuple(xf,yf);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
