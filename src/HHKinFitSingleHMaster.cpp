#include "../include/HHKinFitSingleHMaster.h"

#include "../include/HHEventRecordSingleH.h"
#include "../include/HHKinFitSingleH.h"
#include "../include/HHParticleList.h"
#include "../include/HHPID.h"
#include "../include/HHV4Vector.h"
#include "TMatrixD.h"

#include "TRandom3.h"

#include <TMath.h>
#include <cmath>
#include <cstdlib>
#include <iterator>

void
HHKinFitSingleHMaster::doFullFit()
{
  //Setup event
  HHParticleList* particlelist = new HHParticleList();
  HHEventRecordSingleH eventrecord_rec(particlelist);

  eventrecord_rec.UpdateEntry(HHEventRecordSingleH::tauvis1)->SetVector(*m_tauvis1);
  eventrecord_rec.UpdateEntry(HHEventRecordSingleH::tauvis2)->SetVector(*m_tauvis2);
  
  if (!m_advancedBalance){
    eventrecord_rec.UpdateEntry(HHEventRecordSingleH::MET)->SetEEtaPhiM(m_simpleBalancePt,0,0,0);
    eventrecord_rec.UpdateEntry(HHEventRecordSingleH::MET)->SetErrors(m_simpleBalanceUncert,0,0);
  }
  else{
    if ((m_MET != NULL) && (m_MET_COV.IsValid())){
      eventrecord_rec.UpdateEntry(HHEventRecordSingleH::MET)->SetVector(*m_MET);
      eventrecord_rec.UpdateEntry(HHEventRecordSingleH::MET)->SetCov(m_MET_COV);
    }
  }

  //loop over all hypotheses
  for(std::vector<Int_t>::iterator mh = m_mh.begin(); mh != m_mh.end(); mh++){
      particlelist->UpdateMass(HHPID::h1, *mh);

      HHKinFitSingleH advancedfitter(&eventrecord_rec);
      advancedfitter.SetPrintLevel(0);
      advancedfitter.SetLogLevel(0);
      advancedfitter.SetAdvancedBalance(m_advancedBalance);
      advancedfitter.Fit();
      
      Double_t chi2_full = advancedfitter.GetChi2();
      Double_t prob_full = TMath::Prob(chi2_full,1);
      std::pair< Int_t, Double_t > entry_chi2_full (*mh, chi2_full);
      std::pair< Int_t, Double_t > entry_fitprob_full (*mh, prob_full);
      std::pair< Int_t, Double_t > entry_pullbalance_full (*mh, advancedfitter.GetPullBalance());
      std::pair< Int_t, Double_t > entry_pullbalance_fullX (*mh, advancedfitter.GetPullBalanceX());
      std::pair< Int_t, Double_t > entry_pullbalance_fullY (*mh, advancedfitter.GetPullBalanceY());
      std::pair< Int_t, Int_t >    entry_convergence_full (*mh, advancedfitter.GetConvergence());
      m_fullFitResultChi2.insert(entry_chi2_full);
      m_fullFitResultFitProb.insert(entry_fitprob_full);
      m_fullFitPullBalance.insert(entry_pullbalance_full);
      m_fullFitPullBalanceX.insert(entry_pullbalance_fullX);
      m_fullFitPullBalanceY.insert(entry_pullbalance_fullY);
      m_fullFitConvergence.insert(entry_convergence_full);
      if (chi2_full<m_bestChi2FullFit) {
        m_bestChi2FullFit = chi2_full;
        m_bestHypoFullFit = *mh;
      }
      
      m_tau1_fitted = advancedfitter.GetFitParticle(HHEventRecordSingleH::tau1);
      m_tau2_fitted = advancedfitter.GetFitParticle(HHEventRecordSingleH::tau2);
      m_fixedCovMatrix = advancedfitter.m_fixedCovMatrix;
      
  }
  
  delete particlelist;
}



HHKinFitSingleHMaster::HHKinFitSingleHMaster(TLorentzVector* tauvis1, TLorentzVector* tauvis2, Bool_t truthinput, TLorentzVector* heavyhiggsgen):
    m_mh(std::vector<Int_t>()),


    m_tauvis1(tauvis1),
    m_tauvis2(tauvis2),

    m_MET(NULL),
    m_MET_COV(TMatrixD(2,2)),

    m_truthInput(truthinput),
    m_advancedBalance(false),
    m_simpleBalancePt(0.0),
    m_simpleBalanceUncert(10.0),
    m_fullFitResultChi2(std::map< Int_t , Double_t>()),
    m_bestChi2FullFit(999),
    m_bestHypoFullFit(-1)
{
  if (m_truthInput){
    TRandom3 r(0);   
    

    TLorentzVector* recoil;
    if(heavyhiggsgen != NULL){
       Double_t pxRecoil = r.Gaus(-(heavyhiggsgen->Px() ), 10.0);
       Double_t pyRecoil = r.Gaus(-(heavyhiggsgen->Py() ), 10.0);
       std::cout << "Higgs Recoil X smeared by: " << pxRecoil + heavyhiggsgen->Px() << std::endl;
       std::cout << "Higgs Recoil Y smeared by: " << pyRecoil + heavyhiggsgen->Py() << std::endl;
       recoil = new TLorentzVector(pxRecoil,pyRecoil,0,sqrt(pxRecoil*pxRecoil+pyRecoil*pyRecoil));
    }
    else{
      recoil = new TLorentzVector(0,0,0,0);
      std::cout << "WARNING! Truthinput mode active but no Heavy Higgs gen-information given! Setting Recoil to Zero!" << std::endl;  
    }
    
    TMatrixD recoilCov(2,2);
    recoilCov(0,0)=100;  recoilCov(0,1)=0;
    recoilCov(1,0)=0;    recoilCov(1,1)=100;

    TLorentzVector* met = new TLorentzVector(-(*tauvis1 + *tauvis2 + *recoil));
    
    TMatrixD metCov(2,2);
    metCov = recoilCov;
    
    setAdvancedBalance(met, metCov);
    m_met_smeared = *met;

    delete recoil;
  }
}

Double_t
HHKinFitSingleHMaster::getBestChi2FullFit()
{
  return m_bestChi2FullFit;
}

std::map< Int_t, Double_t >
HHKinFitSingleHMaster::getChi2FullFit(){
  return m_fullFitResultChi2;
}

std::map< Int_t, Double_t >
HHKinFitSingleHMaster::getFitProbFullFit(){
  return m_fullFitResultFitProb;
}

std::map< Int_t, Double_t >
HHKinFitSingleHMaster::getPullBalanceFullFit(){
  return m_fullFitPullBalance;
}

std::map< Int_t, Double_t >
HHKinFitSingleHMaster::getPullBalanceFullFitX(){
  return m_fullFitPullBalanceX;
}

std::map< Int_t, Double_t >
HHKinFitSingleHMaster::getPullBalanceFullFitY(){
  return m_fullFitPullBalanceY;
}

std::map< Int_t, Int_t >
HHKinFitSingleHMaster::getConvergenceFullFit(){
  return m_fullFitConvergence;
}

Int_t
HHKinFitSingleHMaster::getBestHypoFullFit()
{
  return m_bestHypoFullFit;
}


void
HHKinFitSingleHMaster::addMhHypothesis(Double_t m1, Double_t m2, Double_t m3, Double_t m4, Double_t m5, Double_t m6, Double_t m7, Double_t m8, Double_t m9, Double_t m10)
{
  if (m1 != 0) m_mh.push_back(m1);
  if (m2 != 0) m_mh.push_back(m2);
  if (m3 != 0) m_mh.push_back(m3);
  if (m4 != 0) m_mh.push_back(m4);
  if (m5 != 0) m_mh.push_back(m5);
  if (m6 != 0) m_mh.push_back(m6);
  if (m7 != 0) m_mh.push_back(m7);
  if (m8 != 0) m_mh.push_back(m8);
  if (m9 != 0) m_mh.push_back(m9);
  if (m10 != 0) m_mh.push_back(m10);
}

void
HHKinFitSingleHMaster::addMhHypothesis(std::vector<Int_t> v)
{
  m_mh = v;
}

void
HHKinFitSingleHMaster::setAdvancedBalance(TLorentzVector* met, TMatrixD met_cov)
{
  m_advancedBalance = true;
  m_MET = met;
  m_MET_COV = met_cov;
}

void
HHKinFitSingleHMaster::setSimpleBalance(Double_t balancePt, Double_t balanceUncert)
{
  m_advancedBalance = false;
  m_simpleBalancePt = balancePt;
  m_simpleBalanceUncert = balanceUncert;
}
