/*
 * HHHHDiJetKinFitNewMini.cpp
 *
 *  Created on: Jun 17, 2014
 *      Author: vormwald
 */

#include "../include/HHDiJetKinFitNewMini.h"
#include "../include/PSMath.h"
#include "../include/PSTools.h"
#include "../include/BVLineSearch.h"
#include "../include/BVNewtonStep.h"

#include "TString.h"
#include "TPad.h"
#include "TPolyMarker.h"
#include "TMarker.h"
#include "TMatrixD.h"
#include "TLatex.h"
#include "TH2F.h"
#include "TFile.h"
#include <iostream>

HHDiJetKinFitNewMini::HHDiJetKinFitNewMini(HHEventRecord* recrecord)
  : m_recrecord(recrecord), 
    m_fitrecord (new HHEventRecord(*recrecord, "Fit")),
    m_logLevel(0),
    m_chi2(-1), 
    m_chi2_b1(-1), 
    m_chi2_b2(-1), 
    m_convergence(0), 
    m_printlevel(1), 
    m_graphicslevel(1),
    m_maxloops(100)
{
  m_particlelist = m_recrecord->GetParticleList();
}


HHDiJetKinFitNewMini::~HHDiJetKinFitNewMini(){
  delete m_fitrecord;
}


Double_t
HHDiJetKinFitNewMini::CalcChi2(TVectorD par){
  m_fitrecord->UpdateEntry(HHEventRecord::b1)->SetEkeepM(par[0]); // update 4-vectors with fit parameters
  ConstrainE2(HHEventRecord::hb, HHEventRecord::b1, HHEventRecord::b2);

  m_chi2_b1 = Chi2V4(HHEventRecord::b1);
  m_chi2_b2 = Chi2V4(HHEventRecord::b2);
  m_chi2 = m_chi2_b1 + m_chi2_b2; // chi2 calculation
  
  return(m_chi2);
}


void
HHDiJetKinFitNewMini::Fit()
{
  m_recrecord->Recombine();
  m_fitrecord->UpdateMothers(HHEventRecord::b1);
  
  int dim = 1;
  TVectorD start(dim);     start[0]=m_fitrecord->GetEntry(HHEventRecord::b1)->E();  // energy of first b-jet
  TVectorD direction(dim); direction[0]=1;
  TVectorD min(dim);       min[0]=start[0]-m_fitrecord->GetEntry(HHEventRecord::b1)->dE() * 5.;
  if (min[0]<0) min[0]=0;
  TVectorD max(dim);       max[0]=start[0]+m_fitrecord->GetEntry(HHEventRecord::b1)->dE() * 5.;
  
  int iter=0;
  bool converged=false;
  TVectorD result(dim);
  
  TFile f("out.root","recreate");
  while (iter<30 && !converged){
    BVLineSearch linesearch(this,start,direction,min,max);
    bool minimumonline = linesearch.run();
    if (!minimumonline) break;
    BVNewtonStep newtonstep(this,linesearch.getEndPoint(),min,max);
    newtonstep.run();
    start=linesearch.getEndPoint();
    direction=newtonstep.getDirection();
    converged = (fabs(newtonstep.getdF())<0.01 && newtonstep.getd()<0.01);

    f.mkdir(Form("%u",iter));
    f.cd(Form("%u",iter));
//     linesearch.getGraphF().Write();
//     linesearch.getGraphLine().Write();
//     linesearch.getGraphStartPoint().Write();
//     linesearch.getGraphEndPoint().Write();
    linesearch.getGraphFOnLine().Write();
    linesearch.getGraphStartPointOnLine().Write();
    linesearch.getGraphEndPointOnLine().Write();
     
    std::cout << "----------------------------------------------------------------------------------" << std::endl;
    result=start;
    iter++;
  }
  f.Close();
 
  if (converged){
    std::cout << "Minimisation converged in " << iter << " global iterations to:" << std::endl;
    result.Print();
  }
  else
    std::cout << "ERROR: Minimisation did not converge in " << iter << " global iterations" <<std::endl;
  
  

  if (m_logLevel>0){
    std::cout << "chi2 b1: " << m_chi2_b1 << std::endl;
    std::cout << "chi2 b2: " << m_chi2_b2 << std::endl;
    std::cout << "------------------" << std::endl;
    std::cout << "chi2 : " << m_chi2 << std::endl;
    std::cout << "------------------" << std::endl;
  }
}

void
HHDiJetKinFitNewMini::ConstrainE2(Int_t iv4, Int_t iv41, Int_t iv42)
{
  m_fitrecord->UpdateMothers(iv41);
  Double_t M1, M2, M, Mc, E1, E2, b2, beta2;
  M = m_fitrecord->GetEntry(iv4)->M();
  Mc = m_particlelist->GetParticleInfo( m_fitrecord->GetEntry(iv4)->ID() )->M();
  E1 = m_fitrecord->GetEntry(iv41)->E();
  M1 = m_fitrecord->GetEntry(iv41)->M();
  E2 = m_fitrecord->GetEntry(iv42)->E();
  M2 = m_fitrecord->GetEntry(iv42)->M();

  if (M2 == 0.) { // massless case
    m_fitrecord->UpdateEntry(iv42)->SetEkeepM(E2 * (Mc / M) * (Mc / M)); // only changes absolute value and keeps eta, phi, m untouched
                                                                         // such that invariant mass of 1 and 2 gives correct mass Mc
    m_fitrecord->UpdateMothers(iv42);
    return;
  }

  Int_t ID2jet = -1;
  Int_t ID2 = m_fitrecord->GetEntry(iv42)->ID();

  if (ID2 == HHPID::q || ID2 == HHPID::c
      || ID2 == HHPID::b || ID2 == HHPID::gluon)
    ID2jet = 1;

  beta2 = sqrt(E2 * E2 - M2 * M2) / E2;
  if (ID2jet < 0) { // is not a jet
    b2 = (M2 / E2) * (M2 / E2); // isn't the no jet case identical to the jet case? 1/gamma**2 = m**2/E**2 = 1-beta**2
  }
  else { // is a jet
    b2 = 1 - beta2 * beta2;
  }

  Double_t d = (M * M - M1 * M1 - M2 * M2) / (2. * E1 * E2);
  Double_t E2lin = (Mc * Mc - M1 * M1) / (2. * E1 * d);
  Double_t E2N = E1 * d / b2;
  Double_t E2new = E2N * (-1. + sqrt(1. + 2. * E2lin / E2N));

  if (ID2jet < 0) { // is not a jet
    m_fitrecord->UpdateEntry(iv42)->SetEkeepM(E2new);
  }
  else { // is a jet
    HHV4Vector* entry = m_fitrecord->GetEntry(iv42);
    m_fitrecord->UpdateEntry(iv42)->SetEEtaPhiM(E2new, entry->Eta(), entry->Phi(), M2 * E2new / E2);
  }

  m_fitrecord->UpdateMothers(iv42);
}

Double_t
HHDiJetKinFitNewMini::Chi2V4(Int_t iv4)
{
  Double_t chi2 = 0;
  Double_t chi2_E = 0;
  Double_t chi2_Eta = 0;
  Double_t chi2_Phi = 0;
  HHV4Vector* fitobject = m_fitrecord->GetEntry(iv4);
  HHV4Vector* recobject = m_recrecord->GetEntry(iv4);

  //  std::cout << "V4chi2 dE " ; V4fit[iv4].Print();
  if (fitobject->dE() > 0.) {
    chi2_E = ((recobject->E() - fitobject->E()) / fitobject->dE())
        * ((recobject->E() - fitobject->E()) / fitobject->dE());
  }

  if (fitobject->dEta() > 0.) {
    chi2_Eta = ((recobject->Eta() - fitobject->Eta()) / fitobject->dEta())
        * ((recobject->Eta() - fitobject->Eta()) / fitobject->dEta());
  }
  if (fitobject->dPhi() > 0.) {
    chi2_Phi = ((recobject->Phi() - fitobject->Phi()) / fitobject->dPhi())
        * ((recobject->Phi() - fitobject->Phi()) / fitobject->dPhi());
  }

  chi2 = chi2_E + chi2_Eta + chi2_Phi;
  if (m_logLevel>1){
    PSTools::coutf(5, TString("E"));
    PSTools::coutf(9, 2, recobject->E());
    PSTools::coutf(9, 2, fitobject->E());
    PSTools::coutf(9, 2, recobject->E()-fitobject->E());
    PSTools::coutf(3, TString("|"));
    PSTools::coutf(9, 2, fitobject->dE());
    PSTools::coutf(3, TString("->"));
    PSTools::coutf(9, 2, chi2_E);
    std::cout << std::endl;

    PSTools::coutf(5, TString("Eta"));
    PSTools::coutf(9, 2, recobject->Eta());
    PSTools::coutf(9, 2, fitobject->Eta());
    PSTools::coutf(9, 2, recobject->Eta()-fitobject->Eta());
    PSTools::coutf(3, TString("|"));
    PSTools::coutf(9, 2, fitobject->dEta());
    PSTools::coutf(3, TString("->"));
    PSTools::coutf(9, 2, chi2_Eta);
    std::cout << std::endl;

    PSTools::coutf(5, TString("Phi"));
    PSTools::coutf(9, 2, recobject->Phi());
    PSTools::coutf(9, 2, fitobject->Phi());
    PSTools::coutf(9, 2, recobject->Phi()-fitobject->Phi());
    PSTools::coutf(3, TString("|"));
    PSTools::coutf(9, 2, fitobject->dPhi());
    PSTools::coutf(3, TString("->"));
    PSTools::coutf(9, 2, chi2_Phi);
    std::cout << std::endl;
  }

  return chi2;
}

void
HHDiJetKinFitNewMini::SetPrintLevel(Int_t level){
  m_printlevel = level;
}

void
HHDiJetKinFitNewMini::SetGraphicsLevel(Int_t level){
  m_graphicslevel=level;
}

void
HHDiJetKinFitNewMini::SetMaxLoops(Int_t loops){
  m_maxloops = loops;
}

Double_t
HHDiJetKinFitNewMini::GetChi2()
{
  return m_chi2;
}

Double_t
HHDiJetKinFitNewMini::GetChi2_b1()
{
  return m_chi2_b1;
}

Double_t
HHDiJetKinFitNewMini::GetChi2_b2()
{
  return m_chi2_b2;
}

Double_t
HHDiJetKinFitNewMini::GetConvergence()
{
  return m_convergence;
}

HHEventRecord*
HHDiJetKinFitNewMini::GetFitrecord(){
  return m_fitrecord;
}

void
HHDiJetKinFitNewMini::SetLogLevel(Int_t level)
{
  m_logLevel = level;
}

double
HHDiJetKinFitNewMini::GetPullE(Int_t iv4){
  HHV4Vector* fv = m_fitrecord->GetEntry(iv4);
  HHV4Vector* rv = m_recrecord->GetEntry(iv4);
  return (fv->E()-rv->E())/rv->dE();
}

