/*
 * HHEventRecordSingleH.h
 *
 *  Created on: 16.06.2014
 */

#ifndef EVENTRECORD_H_
#define EVENTRECORD_H_

#include "HHKinFit/HHKinFit/interface/HHV4Vector.h"
#include "HHKinFit/HHKinFit/interface/HHParticleList.h"
#include <vector>
#include <TString.h>


class HHEventRecordSingleH{
public:
  enum entry{
    undef = -1,
    MET = 0,   // MET
    htau = 1,  // h --> tau1 tau2
    tau1 = 2,  // tau1
    tau2 = 3,  // tau2
    tauvis1   = 4,   // mu1
    tauinvis1 = 5,  // nu1
    tauvis2   = 6,  // mu2
    tauinvis2 = 7  // nu2
  };

  HHEventRecordSingleH (HHParticleList* particlelist);
  HHEventRecordSingleH (const HHEventRecordSingleH& eventrecord, TString suffix="");
  ~HHEventRecordSingleH ();

  Int_t AddInitialEntry (HHPID::pid id);
  void MakeHEvent ();
  void DecayChain ();
  Int_t GetNEntries () const;
  HHV4Vector* GetEntry (Int_t i) const;
  HHParticleList* GetParticleList () const;
  HHV4Vector* UpdateEntry (Int_t i);
  HHV4Vector* UpdateEntry (entry i) {int temp = i; return UpdateEntry(temp);};
  void Print(TString name="",Int_t mode=1);
  void EventDisplayXY(Int_t style);
  void Recombine();
  void UpdateMothers(Int_t ivDaughter);


private:
  Int_t AddEntry (HHPID::pid id, Int_t origin);
  std::vector< HHV4Vector* > * m_eventrecord;
  HHParticleList* m_particlelist;
};


#endif /* EVENTRECORD_H_ */
