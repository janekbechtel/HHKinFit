/*
 * HHKinFitSingleH.cpp
 *
 *  Created on: Jun 17, 2014
 */

#ifdef HHKINFIT_STANDALONE
#include "../interface/HHKinFitSingleH.h"
#include "../interface/PSMath.h"
#include "../interface/PSTools.h"
#else
#include "HHKinFit/HHKinFit/interface/HHKinFitSingleH.h"
#include "HHKinFit/HHKinFit/interface/PSMath.h"
#include "HHKinFit/HHKinFit/interface/PSTools.h"
#endif

#include "TString.h"
#include "TPad.h"
#include "TPolyMarker.h"
#include "TMarker.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TLatex.h"
#include "TH2F.h"
#include <iostream>

HHKinFitSingleH::HHKinFitSingleH(HHEventRecordSingleH* recrecord)
    : m_chi2(-1),
      m_convergence(0),
      m_covRecoil(2,2),
      m_printlevel(1), m_graphicslevel(0),
      m_maxloops(500),
      m_advancedBalance(kTRUE),
      m_logLevel(0),
      m_recrecord(recrecord), m_fitrecord (new HHEventRecordSingleH(*recrecord, "Fit"))
{
  m_particlelist = m_recrecord->GetParticleList();
}


HHKinFitSingleH::~HHKinFitSingleH(){
  delete m_fitrecord;
}


void
HHKinFitSingleH::Fit()
{
  //  ----------  for PSfit -----
  const Int_t np = 1;
  Double_t a[np];
  Double_t astart[np];
  Double_t alimit[np][2];
  Double_t aprec[np];
  Double_t daN[np];
  Double_t h[np];
  Double_t chi2iter[1], aMemory[np][5], g[np], H[np * np], Hinv[np * np];
  Bool_t noNewtonShifts = false;

  Int_t iter = 0;             //  number of iterations
  Int_t method = 1;           //  initial fit method, see PSfit()
  Int_t mode = 1;             //  mode =1 for start of a new fit by PSfit()
//   Int_t icallNewton = -1;     //  init start of Newton Method
//   Int_t iloop = 0;            // counter for falls to fit function

  // calculate htau from tauvis; recombine leaves measured entries in event record untouched
  m_recrecord->Recombine();

  Double_t tau1LowerLimit = 0.7*m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->E();
  Double_t tau2LowerLimit = 0.7*m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->E();
  Double_t mtau  = m_particlelist->GetParticleInfo(HHPID::tau)->M();

  //Calculate Recoil CovMatrix
  TMatrixD Cov_MET     = m_recrecord->GetEntry(HHEventRecordSingleH::MET)->GetCov();
  TMatrixD Cov_tauvis1 = m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->GetCov();
  TMatrixD Cov_tauvis2 = m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->GetCov();
  m_covRecoil = Cov_MET - (Cov_tauvis1 + Cov_tauvis2);
  
  TMatrixDEigen eigenmatrix(m_covRecoil);

  if (eigenmatrix.GetEigenValues()(0,0)<0 || eigenmatrix.GetEigenValues()(1,1)<0){
    m_fixedCovMatrix = true;
    m_covRecoil(0,0) = 100;
    m_covRecoil(1,1) = 100;
    m_covRecoil(1,0) = 0;
    m_covRecoil(0,1) = 0;
  }
  
  // initialise tau1 and tau2 vectors
  if(m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->E() > mtau){
    m_fitrecord->UpdateEntry(HHEventRecordSingleH::tau1)->SetEEtaPhiM(m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->E(),
                                                               m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->Eta(),      
                                                               m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->Phi(),
                                                               mtau);
  }
  else{
    m_fitrecord->UpdateEntry(HHEventRecordSingleH::tau1)->SetEEtaPhiM(pow(m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->P(),2)+pow(mtau,2),
                                                               m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->Eta(),
                                                               m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->Phi(),
                                                               mtau);
  }

  if(tau2LowerLimit > mtau){
    m_fitrecord->UpdateEntry(HHEventRecordSingleH::tau2)->SetEEtaPhiM(tau2LowerLimit,
                                                               m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->Eta(),
                                                               m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->Phi(),
                                                               mtau);
  }
  else{
    m_fitrecord->UpdateEntry(HHEventRecordSingleH::tau2)->SetEEtaPhiM(pow(m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->P(),2)+pow(mtau,2),
                                                               m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->Eta(),
                                                               m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->Phi(),
                                                               mtau);
  }

  //firstly, compute upper limit of E(tau1) by having set E(tau2)=E(tau2)_min and compute E(tau1)
  ConstrainE2(HHEventRecordSingleH::htau, HHEventRecordSingleH::tau2, HHEventRecordSingleH::tau1);
  Double_t maxEtau1 = m_fitrecord->GetEntry(HHEventRecordSingleH::tau1)->E();
  Double_t minEtau1 = tau1LowerLimit;

  //Reset taus
  if(m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->E() > mtau){
    m_fitrecord->UpdateEntry(HHEventRecordSingleH::tau2)->SetEEtaPhiM(m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->E(),
                                                               m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->Eta(),
                                                               m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->Phi(),
                                                               mtau);
  }
  else{
    m_fitrecord->UpdateEntry(HHEventRecordSingleH::tau2)->SetEEtaPhiM(pow(m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->P(),2)+pow(mtau,2),
                                                               m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->Eta(),
                                                               m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->Phi(),
                                                               mtau);
  }

  if(m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->E() > mtau){
    m_fitrecord->UpdateEntry(HHEventRecordSingleH::tau1)->SetEEtaPhiM(m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->E(),
                                                               m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->Eta(),
                                                               m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->Phi(),
                                                               mtau);
  }
  else{
    m_fitrecord->UpdateEntry(HHEventRecordSingleH::tau1)->SetEEtaPhiM(pow(m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->P(),2)+pow(mtau,2),
                                                               m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->Eta(),
                                                               m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->Phi(),
                                                               mtau);
  }

  ConstrainE2(HHEventRecordSingleH::htau, HHEventRecordSingleH::tau1, HHEventRecordSingleH::tau2);
  m_fitrecord->UpdateMothers(HHEventRecordSingleH::tau1);

  if (!(maxEtau1>=m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->E()) ){
    if(m_printlevel > 0)
      std::cout << "tautau mass constraint cannot be fulfilled -> reconstructed visible tau energy greater/smaller than maximal/minimal allowed total tau energy." << std::endl;
    m_convergence=-1;
    m_chi2=9999;
    return;
  }

  // fill initial tau fit parameters
  astart[0] = m_fitrecord->GetEntry(HHEventRecordSingleH::tau1)->E();       // energy of first tau
  aprec[0] = 0.1;               //0.1                 // precision for fit
  // fill initial step width
  h[0] = 0.1*minEtau1;   //                
  daN[0] = 1.0;   //0.0                 // initial search direction

  alimit[0][0] = minEtau1;              // tau: minimum is visible tau1 energy
  alimit[0][1] = maxEtau1;              //      maximum as computed above

  // tau: check initial values against fit range
  if (astart[0] - h[0] < alimit[0][0]) {
    astart[0] = alimit[0][0] + h[0];
  }
  else if (astart[0] + h[0] > alimit[0][1]) {
    astart[0] = alimit[0][1] - h[0];
  }

  for (Int_t ip = 0; ip < np; ip++) {
    a[ip] = astart[ip];
  }
  for (Int_t ip = 0; ip < np; ip++) {
    aMemory[ip][0] = -999.0;
    aMemory[ip][1] = -995.0;
    aMemory[ip][2] = -990.0;
    aMemory[ip][3] = -985.0;
    aMemory[ip][3] = -980.0;
  }

  //  static const Int_t nloopmax = 100;
  //  static Double_t Ya[nloopmax];
  //  static Double_t Ya1[nloopmax];
  //  static Double_t HPx[nloopmax], HPy[nloopmax];
  //  static Double_t HPx1[nloopmax],HPy1[nloopmax];

  ConstrainE2(HHEventRecordSingleH::htau, HHEventRecordSingleH::tau1, HHEventRecordSingleH::tau2);

  if (m_logLevel>1){
    std::cout << "Starting FitLoop! Start-Values: " << std::endl;
    std::cout << "Fit Params: " << std::endl;
    std::cout << "Etau1: " << m_fitrecord->GetEntry(HHEventRecordSingleH::tau1)->E() << " Etau2: " << m_fitrecord->GetEntry(HHEventRecordSingleH::tau2)->E() << std::endl;
    std::cout << "ETau1 Lower Limit: " << alimit[0][0] << "ETau1 Upper Limit: " << alimit[0][1] << std::endl; 
    std::cout << "Etau1 Precision is : " << aprec[0] << std::endl;
    std::cout << "HH Px: " << m_fitrecord->GetEntry(HHEventRecordSingleH::htau)->Px() << " HH Py: " << m_fitrecord->GetEntry(HHEventRecordSingleH::htau)->Py() << std::endl;
  }

  for (Int_t iloop = 0; iloop < m_maxloops * 10 && iter < m_maxloops; iloop++) { // FIT loop
    
    m_fitrecord->UpdateEntry(HHEventRecordSingleH::tau1)->SetEkeepM(a[0]);
    ConstrainE2(HHEventRecordSingleH::htau, HHEventRecordSingleH::tau1, HHEventRecordSingleH::tau2);

    m_chi2 = Chi2Balance();

    if (m_logLevel>1){
      std::cout << std::setprecision(6);
      std::cout << "chi2 : " << m_chi2 << std::endl;
      std::cout << "------------------" << std::endl;
    }

    //    if (iloop >= 0 && iloop < nloopmax) {
    //      Ya1[iloop] = m_fitrecord->GetEntry(HHEventRecordSingleH::tau1)->E();
    //    }
    //    if (iter >= 0 && iter < nloopmax) {
    //      Ya[iter] = m_fitrecord->GetEntry(HHEventRecordSingleH::tau1)->E();
    //    }
    //    if (iloop >= 0 && iloop < nloopmax) {
    //      HPx1[iloop] = m_fitrecord->GetEntry(HHEventRecordSingleH::htau)->Px();
    //      HPy1[iloop] = m_fitrecord->GetEntry(HHEventRecordSingleH::htau)->Py();
    //    }
    //    if (iter >= 0 && iter < nloopmax) {
    //      HPx[iter] = m_fitrecord->GetEntry(HHEventRecordSingleH::htau)->Px();
    //      HPy[iter] = m_fitrecord->GetEntry(HHEventRecordSingleH::htau)->Py();
    //    }

    PSMath::PSfitShow(iloop, m_convergence, iter, method, mode, m_printlevel,
                      m_graphicslevel, np, a, astart, alimit, aprec, daN, h,
                      m_chi2, g, H);

    if (m_convergence != 0) break;
    m_convergence = PSMath::PSfit(iloop, iter, method, mode, noNewtonShifts, m_printlevel,
                                  np, a, astart, alimit, aprec,
                                  daN, h, aMemory, m_chi2, chi2iter, g, H,
                                  Hinv);
  }
  // ------ end of FIT loop
  
  if(m_convergence != 0 && m_convergence != 5){
    if(a[0] < (alimit[0][0] + 2*aprec[0]) ){
      if(m_convergence == 3)
        m_convergence = 5;
      else{
        if (m_logLevel>1) std::cout << "Convergence at lower tau limit!" << std::endl;
        m_convergence = 4;
      }
    }
    if(a[0] > (alimit[0][1] - 2*aprec[0]) ){
      if(m_convergence == 3)
	m_convergence = 5;
      else{
	if (m_logLevel>1)
	  std::cout << "Convergence at upper tau limit!" << std::endl;
	m_convergence = 4;
      }
    }
  }
  if (m_logLevel>1)
    std::cout << "Convergence is " << m_convergence << std::endl;
   
  if (m_logLevel>0){
    Double_t px_MET = m_recrecord->GetEntry(HHEventRecordSingleH::MET)->Px();
    Double_t py_MET = m_recrecord->GetEntry(HHEventRecordSingleH::MET)->Py();

    Double_t px_tauvis1 = m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->Px();
    Double_t py_tauvis1 = m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->Py();

    Double_t px_tauvis2 = m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->Px();
    Double_t py_tauvis2 = m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->Py();

    Double_t px_H_reco = px_MET + px_tauvis1 + px_tauvis2;
    Double_t py_H_reco = py_MET + py_tauvis1 + py_tauvis2;

    std::cout << "MET:" << std::endl;
    PSTools::coutf(9, px_MET);
    PSTools::coutf(9, py_MET); std::cout << std::endl;

    std::cout << "recoil:" << std::endl;
    PSTools::coutf(9, -px_H_reco);
    PSTools::coutf(9, -py_H_reco); std::cout << std::endl;

    std::cout << "CHI2 TOTAL:   " << m_chi2 << std::endl;
  }
}

void
HHKinFitSingleH::ConstrainE2(Int_t iv4, Int_t iv41, Int_t iv42)
{
  m_fitrecord->UpdateMothers(iv41);
  Double_t Mc = m_particlelist->GetParticleInfo( m_fitrecord->GetEntry(iv4)->ID() )->M();
  Double_t M1c = m_particlelist->GetParticleInfo( m_fitrecord->GetEntry(iv41)->ID() )->M();
  Double_t M2c = m_particlelist->GetParticleInfo( m_fitrecord->GetEntry(iv42)->ID() )->M();
  Double_t M = m_fitrecord->GetEntry(iv4)->M();

  Double_t C = 0.5*(pow(Mc,2)-pow(M1c,2)-pow(M2c,2));
  
  int loopCount = 0;
  while(fabs(M-Mc) > 0.000001){
    loopCount++;
    Double_t P1x = m_fitrecord->GetEntry(iv41)->Px();
    Double_t P1y = m_fitrecord->GetEntry(iv41)->Py();
    Double_t P1z = m_fitrecord->GetEntry(iv41)->Pz();
    Double_t P1 = m_fitrecord->GetEntry(iv41)->P();
    Double_t E1 =  m_fitrecord->GetEntry(iv41)->E();
  
    Double_t P2x = m_fitrecord->GetEntry(iv42)->Px();
    Double_t P2y = m_fitrecord->GetEntry(iv42)->Py();
    Double_t P2z = m_fitrecord->GetEntry(iv42)->Pz();
    Double_t P2 = m_fitrecord->GetEntry(iv42)->P();

    Double_t cosa = (P1x*P2x+P1y*P2y+P1z*P2z)/(P1*P2);
    Double_t E2new = -1;

    if(cosa==0){
      E2new=C/E1;
    }
    else{
      Double_t cp=C/(cosa*P1);
      Double_t dp=E1/(cosa*P1);
      Double_t a=pow(dp,2)-1;
      Double_t b=-2*dp*cp;
      Double_t c=pow(cp,2)+pow(M2c,2);

      if (cosa>0) E2new = (1./(2*a))*(-b+sqrt(pow(b,2)-4*a*c));
      if (cosa<0) E2new = (1./(2*a))*(-b-sqrt(pow(b,2)-4*a*c));
    }

    m_fitrecord->UpdateEntry(iv42)->SetEkeepM(E2new);
    m_fitrecord->UpdateMothers(iv42);
    M = m_fitrecord->GetEntry(iv4)->M();

  }

  if(m_logLevel > 1){
    std::cout << "Did Mass Constraint Loop " << loopCount << " times." << std::endl;
  }

}

Double_t
HHKinFitSingleH::Chi2Balance()
{

  //fit objets
  Double_t px_H_fit   = m_fitrecord->GetEntry(HHEventRecordSingleH::htau)->Px();
  Double_t py_H_fit   = m_fitrecord->GetEntry(HHEventRecordSingleH::htau)->Py();


  //reco objects
  Double_t px_MET = m_recrecord->GetEntry(HHEventRecordSingleH::MET)->Px();
  Double_t py_MET = m_recrecord->GetEntry(HHEventRecordSingleH::MET)->Py();

  Double_t px_tauvis1 = m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->Px();
  Double_t py_tauvis1 = m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->Py();

  Double_t px_tauvis2 = m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->Px();
  Double_t py_tauvis2 = m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->Py();
 
  Double_t px_H_reco = px_MET + px_tauvis1 + px_tauvis2;   //+ 1.0 just a test. Remove!
  Double_t py_H_reco = py_MET + py_tauvis1 + py_tauvis2;   //+ 1.0 just a test. Remove!
  Double_t pt_H_reco = sqrt(pow(px_H_reco,2) + pow(py_H_reco,2));

  Double_t res[2];
  res[0] = px_H_fit - px_H_reco;    // residuum in Pt_H
  res[1] = py_H_fit - py_H_reco;    // residuum in Pt_H
  
  Double_t Vxx = m_covRecoil(0,0);
  Double_t Vyy = m_covRecoil(1,1);
  Double_t Vxy = m_covRecoil(0,1);

  Double_t det, Vinvxx, Vinvyy, Vinvxy;
  det = Vxx * Vyy - Vxy * Vxy;
  Vinvxx = Vyy / det;
  Vinvyy = Vxx / det;
  Vinvxy = -Vxy / det;

  Double_t Vinv[2 * 2];
  Vinv[0] = Vinvxx;
  Vinv[1] = Vinvxy;
  Vinv[2] = Vinv[1];
  Vinv[3] = Vinvyy;

  Double_t chi2=-1;
  if (m_advancedBalance){
    chi2 = res[0] * (Vinv[0] * res[0] + Vinv[1] * res[1]) // chi2 = res_transponiert * Vinv * res
         + res[1] * (Vinv[2] * res[0] + Vinv[3] * res[1]);
  }
  else{
    chi2 =  pow(( pt_H_reco - m_fitrecord->GetEntry(HHEventRecordSingleH::htau)->Pt()), 2)
          / pow((m_recrecord->GetEntry(HHEventRecordSingleH::MET)->dE()), 2);
  }

  
  if (m_logLevel>2){
    std::cout << "P_H_Fit - x: " << px_H_fit << "  y: " << py_H_fit << std::endl;
    std::cout << "P_H_reco - x: " << px_H_reco << "  y: " << py_H_reco << std::endl;
    std::cout << "res[0] " << res[0] << "  res[1] " << res[1] << std::endl;
    std::cout << "V[0] " << Vxx << "  V[1] " << Vxy << std::endl;
    std::cout << "V[2] " << Vxy << "  V[3] " << Vyy << std::endl;
    std::cout << "Vinv[0] " << Vinv[0] << "  Vinv[1] " << Vinv[1] << std::endl;
    std::cout << "Vinv[2] " << Vinv[2] << "  Vinv[3] " << Vinv[3] << std::endl;
  }

  if (m_logLevel>1){
    PSTools::coutf(5, TString("pmx"));
    PSTools::coutf(9, 2, px_H_fit);
    PSTools::coutf(9, 2, px_H_reco);
    PSTools::coutf(9, 2, res[0]);
    PSTools::coutf(3, TString("|"));
    PSTools::coutf(9, 2, sqrt(Vxx));
    std::cout << std::endl;

    PSTools::coutf(5, TString("pmy"));
    PSTools::coutf(9, 2, py_H_fit);
    PSTools::coutf(9, 2, py_H_reco);
    PSTools::coutf(9, 2, res[1]);
    PSTools::coutf(3, TString("|"));
    PSTools::coutf(9, 2, sqrt(Vyy));
    PSTools::coutf(3, TString("->"));
    PSTools::coutf(9, 2, chi2);
    std::cout << std::endl;
    std::cout << "--------------" << std::endl;
  }

  return chi2;
}

void
HHKinFitSingleH::SetPrintLevel(Int_t level){
  m_printlevel = level;
}

void
HHKinFitSingleH::SetGraphicsLevel(Int_t level){
  m_graphicslevel=level;
}

void
HHKinFitSingleH::SetMaxLoops(Int_t loops){
  m_maxloops = loops;
}

Double_t
HHKinFitSingleH::GetChi2()
{
  return m_chi2;
}

Double_t
HHKinFitSingleH::GetConvergence()
{
  return m_convergence;
}

HHEventRecordSingleH*
HHKinFitSingleH::GetFitrecord(){
  return m_fitrecord;
}

void
HHKinFitSingleH::SetAdvancedBalance(Bool_t flag)
{
  m_advancedBalance = flag;
}

void
HHKinFitSingleH::SetLogLevel(Int_t level)
{
  m_logLevel = level;
}

double
HHKinFitSingleH::GetPullBalance(){
  Double_t px_MET = m_recrecord->GetEntry(HHEventRecordSingleH::MET)->Px();
  Double_t py_MET = m_recrecord->GetEntry(HHEventRecordSingleH::MET)->Py();

  Double_t px_tauvis1 = m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->Px();
  Double_t py_tauvis1 = m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->Py();

  Double_t px_tauvis2 = m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->Px();
  Double_t py_tauvis2 = m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->Py();

  Double_t px_H_reco = px_MET + px_tauvis1 + px_tauvis2;
  Double_t py_H_reco = py_MET + py_tauvis1 + py_tauvis2;
  Double_t pt_H_reco = sqrt(pow(px_H_reco,2) + pow(py_H_reco,2));


  Double_t pull=-99;
  if (m_advancedBalance){
    if(m_chi2 >=0 ) 
      pull = sqrt(m_chi2);
  }
  else
    pull =  (pt_H_reco - m_fitrecord->GetEntry(HHEventRecordSingleH::htau)->Pt()) / (m_recrecord->GetEntry(HHEventRecordSingleH::MET)->dE());
  
  return pull;
}

double
HHKinFitSingleH::GetPullBalanceX(){

  Double_t px_MET = m_recrecord->GetEntry(HHEventRecordSingleH::MET)->Px();

  Double_t px_tauvis1 = m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->Px();
  Double_t px_tauvis2 = m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->Px();
 
  Double_t px_H_reco = px_MET + px_tauvis1 + px_tauvis2;

  Double_t pull= -999;
  
  if(m_advancedBalance)
    pull =  (px_H_reco - m_fitrecord->GetEntry(HHEventRecordSingleH::htau)->Px())/sqrt(m_covRecoil(0,0) );
  
  return pull;
}

double
HHKinFitSingleH::GetPullBalanceY(){

  Double_t py_MET = m_recrecord->GetEntry(HHEventRecordSingleH::MET)->Py();

  Double_t py_tauvis1 = m_recrecord->GetEntry(HHEventRecordSingleH::tauvis1)->Py();
  Double_t py_tauvis2 = m_recrecord->GetEntry(HHEventRecordSingleH::tauvis2)->Py();
 
  Double_t py_H_reco = py_MET + py_tauvis1 + py_tauvis2;

  Double_t pull= -999;
  
  if(m_advancedBalance)
    pull =  (py_H_reco - m_fitrecord->GetEntry(HHEventRecordSingleH::htau)->Py())/sqrt(m_covRecoil(1,1) );
  
  return pull;
}
