#ifndef CTOOLS_H
#define CTOOLS_H 1

// Call for my libraries
#include "colormod.h"

// Call for standard libraries
#include <iterator>
#include <algorithm>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <numeric>

// Include of ROOT libraries
#include <TROOT.h>
#include "TMath.h"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress (double percentage);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Useful FUNCTIONS
// Template for printing vector
template<typename T>
void printVector(
    const std::vector<T>& t)
{
  for (unsigned int i=0; i<(t).size(); i++)
  {
    std::cout <<  t.at(i) << "\t";
  }
  std::cout << "\n";
}

template<typename T>
void PrintVector(
    const std::vector<T>& t)
{
  for (unsigned int i=0; i<(t).size(); i++)
  {
    std::cout <<  t.at(i) << "\t";
  }
  std::cout << "\n";
}

template <typename firstT, typename... T>
void PrintVectors(const std::vector<firstT>& first, const std::vector<T>& ... vec)
{
     for(int i = 0; i < first.size(); i++)
     {
          std::cout << i << ": " << first.at(i);
          ( (std::cout << ' ' << vec.at(i)), ... );
          std::cout << std::endl;
     }
}

template<typename T>
void printVectorInVector(const T& t) {
    std::for_each(t.cbegin(), t.cend(), PrintVector<typename T::value_type>);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Function to sort a vector and apply the same permutation to other vectors
template <typename T, typename Compare>
void getSortPermutation(
    std::vector<unsigned>& out,
    const std::vector<T>& v,
    Compare compare = std::less<T>())
{
    out.resize(v.size());
    std::iota(out.begin(), out.end(), 0);

    std::sort(out.begin(), out.end(),
        [&](unsigned i, unsigned j){ return compare(v[i], v[j]); });
}

template <typename T>
void applyPermutation(
    const std::vector<unsigned>& order,
    std::vector<T>& t)
{
    assert(order.size() == t.size());
    std::vector<T> st(t.size());
    for(unsigned i=0; i<t.size(); i++)
    {
        st[i] = t[order[i]];
    }
    t = st;
}

template <typename T, typename... S>
void applyPermutation(
    const std::vector<unsigned>& order,
    std::vector<T>& t,
    std::vector<S>&... s)
{
    applyPermutation(order, t);
    applyPermutation(order, s...);
}

// sort multiple vectors using the criteria of the first one
template<typename T, typename Compare, typename... SS>
void sortVectors(
    const std::vector<T>& t,
    Compare comp,
    std::vector<SS>&... ss)
{
    std::vector<unsigned> order;
    getSortPermutation(order, t, comp);
    applyPermutation(order, ss...);
}

// make less verbose for the usual ascending order
template<typename T, typename... SS>
void sortVectorsAscending(
    const std::vector<T>& t,
    std::vector<SS>&... ss)
{
    sortVectors(t, std::less<T>(), ss...);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t VectorMean(std::vector<Double_t>& theVector);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t VectorWeightedMean(std::vector<Double_t>& theVector, std::vector<Double_t>& theWeight);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t VectorStdDev(std::vector<Double_t>& theVector);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Function for debugging
template<typename T>
void Coucou(T i)
{
  std::cout << SetForeMAG << "coucou jsuis la " << i << RESETTEXT << std::endl;
};
template<typename T>
int CoucouF(int i,T toto,int imin)
{
  if (i > imin-1)
  {
    std::cout <<SetForeMAG << toto<<i<<RESETTEXT << std::endl ;
  }
  return 0 ;
};

int SimpleRoundTo(Double_t num);
int SimpleRoundTo(Float_t num);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::tuple<double, double> Rotation(Double_t theta, Double_t theta_tan, Double_t x, Double_t y);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif //MTOBJECT_H
