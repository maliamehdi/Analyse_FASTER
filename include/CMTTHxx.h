/*
 * CMTObject:
 * Base class of CMTTHxx
 *
 * CMTTHxx:
 * Multithreading wrapper for all THxx spectra of root library
 * Author: corentin.hiver@ijclab.in2p3.fr
 *
 * Inspiration :
 * https://root.cern.ch/doc/master/TThreadedObject_8hxx_source.html#l00167
 *
 * MANDATORY : CMTObject::nb_threads must be set to the correct number of threads;
 */
#ifndef CCMTTHxx_H
#define CCMTTHxx_H 1

#include "CMTObject.h"
#include "colormod.h"

template <class THxx>
class CMTTHxx : public CMTObject
{
public:
  CMTTHxx(){};
  ~CMTTHxx()
  {
     for (auto histo: m_collection) if (histo!= nullptr) delete histo;
     if  (m_merged!=nullptr) delete m_merged;
   };

  // General constructor, fully based on reset method
  template <class... ARGS>
  CMTTHxx(std::string name, ARGS &&... args){ this -> reset (name, args...); };

  //Nullptr initialiser, to create a zombie
  template <class... ARGS>
  void reset(std::nullptr_t) {m_exists = false;};

  //General initialiser to construct any root histogram vector
  template <class... ARGS>
  void reset(std::string name, ARGS &&... args)
  {
    m_exists = true;// If SIGSEGV here, did you initialise the container ?
    m_collection.resize(CMTObject::nb_threads, nullptr);
    m_str_name = name;
    for (size_t i = 0; i<m_collection.size(); i++)
    {
      m_collection[i] = new THxx ((name+std::to_string(i)).c_str(), std::forward<ARGS>(args)...);
    }
  };

  operator THxx*() {return m_merged;}

  void Merge();
  void Write(std::string const param = "");

  Bool_t const & exists() {return m_exists;};
  THxx * Merged();
  THxx * Get(UShort_t const & thread_nb) {return m_collection[thread_nb];}
  THxx * operator[](int const & thread_nb) {return m_collection[thread_nb];}
  THxx * operator->() {return m_merged;}

private:
  TFile* m_file = nullptr;
  Bool_t m_exists = false;
  std::string m_str_name = "";
  Bool_t m_is_merged = false;

  TLegend *legend = nullptr;

  std::vector<THxx*> m_collection;
  THxx* m_merged = nullptr;
};

template<class THxx>
void CMTTHxx<THxx>::Merge()
{
  if (!m_exists || m_collection.size() == 0 || m_collection[0] -> IsZombie() || m_collection[0] -> Integral() < 1)
  {
    m_exists = false;
    m_merged = new THxx();
  }
  else
  {
    if (CMTObject::nb_threads == 1)
    {
      m_merged = m_collection[0];
    }
    else
    {
      if ( (m_file = gROOT -> GetFile()) ) gROOT -> cd(); //To get out of scope of any potential TFile
      m_merged = (THxx*) (m_collection[0]->Clone(m_str_name.c_str()));
      for (unsigned int i = 1; i<m_collection.size(); i++)
      {
        m_merged -> Add(m_collection[i]);
      }
      if (m_file) m_file -> cd(); //To return to the scope of any potential TFile
      m_file = nullptr;
    }
    m_is_merged = true;
  }
}

template<class THxx>
void CMTTHxx<THxx>::Write(std::string const param)
{
  if (!m_is_merged) this -> Merge();
  if (m_exists == false) return;
  this -> Merged() -> Write(m_str_name.c_str(), TROOT::kOverwrite);
}

template<class THxx>
THxx* CMTTHxx<THxx>::Merged()
{
  if (!m_is_merged)
  {
    cout << SetBOLD << SetForeRED << endl;
    cout << "No merged histogram yet ! Issues in the code :/" << endl ;
    cout << RESETTEXT << endl;
    return (new THxx());
  }
  return m_merged;
}

template <class THxx>
using Vector_CMTTHxx = std::vector<CMTTHxx<THxx>>;

#endif //CMTTHxx_H
