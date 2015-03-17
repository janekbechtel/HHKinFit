/*
 *HHKinFitSingleH.h
 *
 *  Created on: Jun 17, 2014
 */

#ifndef HHKINFIT_H_
#define HHKINFIT_H_

#include "HHKinFit/interface/HHV4Vector.h"
#include "HHKinFit/interface/HHEventRecordSingleH.h"
#include <map>
#include <utility>
#include <Rtypes.h>
#include <TH1F.h>

class HHKinFitSingleH{
public:
  HHKinFitSingleH(HHEventRecordSingleH* recrecord);
  ~HHKinFitSingleH();

  void Fit();
  void ConstrainE2(Int_t iv4, Int_t iv41, Int_t iv42);
  Double_t Chi2Balance();

  //Getters
  Double_t GetChi2();
  Double_t GetConvergence();

  HHEventRecordSingleH* GetFitrecord();
  double GetPullBalance();
  double GetPullBalanceX();
  double GetPullBalanceY();

  TLorentzVector GetFitParticle(Int_t iv4) { 
    HHV4Vector* fitobject = m_fitrecord->GetEntry(iv4);
    return TLorentzVector(fitobject->Px(), fitobject->Py(), fitobject->Pz(), fitobject->E());
  }
  
  //Setters
  void SetPrintLevel(Int_t level);
  void SetGraphicsLevel(Int_t level);
  void SetMaxLoops(Int_t loops);

  void SetAdvancedBalance(Bool_t flag);
  void SetLogLevel(Int_t level);

  bool m_fixedCovMatrix;
private:
  Double_t m_chi2;
  Int_t m_convergence;
  TMatrixD m_covRecoil;
  Int_t m_printlevel;
  Int_t m_graphicslevel;

  Int_t m_maxloops;
  Bool_t m_advancedBalance;
  Int_t m_logLevel;

  HHEventRecordSingleH* m_recrecord;
  HHEventRecordSingleH* m_fitrecord;
  HHParticleList* m_particlelist;

};



#endif /* HHKINFIT_H_ */
