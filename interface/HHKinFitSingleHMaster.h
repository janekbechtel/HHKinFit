#ifndef HHKinFitSingleHMaster_H
#define HHKinFitSingleHMaster_H

#include <Rtypes.h>
#include <stdio.h>
#include <TMatrixD.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <map>
#include <utility>
#include <vector>
#include <sstream>

#include "TLorentzVector.h"



class HHKinFitSingleHMaster
{
public:
  HHKinFitSingleHMaster(const TLorentzVector* tauvis1, const TLorentzVector* tauvis2, Bool_t truthinput=0, TLorentzVector* heavyhiggsgen=NULL);

  void doFullFit();
  
  //Setters
  void setAdvancedBalance(const TLorentzVector* met, const TMatrixD& met_cov);
  void setSimpleBalance(Double_t balancePt, Double_t balanceUncert);
  
  //Getters for fit results
  Double_t getBestChi2FullFit() const;
  Int_t getBestHypoFullFit() const;
  std::map< Int_t, Double_t > getChi2FullFit() const;
  Double_t getChi2(Int_t mh) const;
  std::map< Int_t, Double_t > getFitProbFullFit() const;
  Double_t getFitProb(Int_t mh) const;
  std::map< Int_t, Double_t > getPullBalanceFullFit() const;
  std::map< Int_t, Double_t > getPullBalanceFullFitX() const;
  std::map< Int_t, Double_t > getPullBalanceFullFitY() const;
  std::map< Int_t, Int_t > getConvergenceFullFit() const;
  TLorentzVector getTau1Fitted(Int_t mh=-1) const; //mh<0 returns best fit
  TLorentzVector getTau2Fitted(Int_t mh=-1) const;
  TLorentzVector getTau1BestFit() const;
  TLorentzVector getTau2BestFit() const;
  std::map<Int_t,TLorentzVector> getTau1FullFit() const;
  std::map<Int_t,TLorentzVector> getTau2FullFit() const;

  //Hypotheses
  void addMhHypothesis(const std::vector<Int_t>& v);
  void addMhHypothesis(Double_t m1, Double_t m2=0, Double_t m3=0, Double_t m4=0, Double_t m5=0, Double_t m6=0, Double_t m7=0, Double_t m8=0, Double_t m9=0, Double_t m10=0);


  TLorentzVector m_tau1_fitted;
  TLorentzVector m_tau2_fitted;
  TLorentzVector m_met_smeared; 
  bool m_fixedCovMatrix; //returns value for last fit
private:
  //hypotheses
  std::vector< Int_t > m_mh;

  //input vectors
  const TLorentzVector* m_tauvis1;
  const TLorentzVector* m_tauvis2;

  const TLorentzVector* m_MET;
  TMatrixD m_MET_COV;

  //full event fit
  Bool_t m_truthInput;
  Bool_t m_advancedBalance;
  Double_t m_simpleBalancePt;
  Double_t m_simpleBalanceUncert;
  std::map< Int_t, Double_t > m_fullFitResultChi2;
  std::map< Int_t, Double_t > m_fullFitResultFitProb;
  std::map< Int_t, Double_t > m_fullFitPullBalance;
  std::map< Int_t, Double_t > m_fullFitPullBalanceX;
  std::map< Int_t, Double_t > m_fullFitPullBalanceY;
  std::map< Int_t, Int_t > m_fullFitConvergence;
  std::map< Int_t, TLorentzVector > m_tau1_fitted_map;
  std::map< Int_t, TLorentzVector > m_tau2_fitted_map;

  Double_t m_bestChi2FullFit;
  Double_t m_bestMHFullFit;
  Int_t m_bestHypoFullFit;
};

#endif
