/*
 * HHDiJetKinFit.cpp
 *
 *  Created on: Jun 17, 2014
 *      Author: vormwald
 */

#include "HHKinFit/HHKinFit/interface/HHDiJetKinFit.h"
#include "HHKinFit/HHKinFit/interface/PSMath.h"
#include "HHKinFit/HHKinFit/interface/PSTools.h"

#include "TString.h"
#include "TPad.h"
#include "TPolyMarker.h"
#include "TMarker.h"
#include "TMatrixD.h"
#include "TLatex.h"
#include "TH2F.h"
#include <iostream>

HHDiJetKinFit::HHDiJetKinFit(HHEventRecord* recrecord)
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

void
HHDiJetKinFit::FitNew()
{
  //  ----------  for PSfit -----
  const Int_t np = 1;
  //Double_t a[np];
  Double_t astart[np];
  Double_t alimit[np][2];

  // calculate htau from tauvis; recombine leaves measured entries in event record untouched
  m_recrecord->Recombine();

  m_fitrecord->UpdateMothers(HHEventRecord::b1);
  
  // fill initial fit parameters
  astart[0] = m_fitrecord->GetEntry(HHEventRecord::b1)->E();         // energy of first b-jet

  // fit range
  alimit[0][0] = astart[0] - m_fitrecord->GetEntry(HHEventRecord::b1)->dE() * 5.; // b-jets: 5 sigma
  if (alimit[0][0]<0) alimit[0][0]=0;
  alimit[0][1] = astart[0] + m_fitrecord->GetEntry(HHEventRecord::b1)->dE() * 5.;
  
  //first grid scan
  Int_t steps1 = 10;
  Double_t stepsize1=(alimit[0][1]-alimit[0][0])/(steps1*1.0);
  Double_t minE1=0;
  Double_t chi2_min=999999;
  
  for (Int_t i=0; i<steps1; i++) { // FIT loop
    Double_t Eb1test=alimit[0][0]+i*stepsize1;
    
    m_fitrecord->UpdateEntry(HHEventRecord::b1)->SetEkeepM(Eb1test); // update 4-vectors with fit parameters
    ConstrainE2(HHEventRecord::hb, HHEventRecord::b1, HHEventRecord::b2);

    m_chi2_b1 = Chi2V4(HHEventRecord::b1);
    m_chi2_b2 = Chi2V4(HHEventRecord::b2);
    m_chi2 = m_chi2_b1 + m_chi2_b2; // chi2 calculation

    if (m_chi2 < chi2_min){
      chi2_min=m_chi2;
      minE1=Eb1test;
    }
    
    if (m_logLevel>1){
      std::cout << "chi2 b1: " << m_chi2_b1 << std::endl;
      std::cout << "chi2 b2: " << m_chi2_b2 << std::endl;
      std::cout << "------------------" << std::endl;
      std::cout << "chi2 : " << m_chi2 << std::endl;
      std::cout << "------------------" << std::endl;
    }
  }
    
  //second grid scan
  steps1 = 20;
  Double_t limit1_up   = ((minE1+2*stepsize1)>alimit[0][1])?alimit[0][1]:(minE1+1*stepsize1);
  Double_t limit1_down = ((minE1-2*stepsize1)<alimit[0][0])?alimit[0][0]:(minE1-1*stepsize1);
  stepsize1=(limit1_up-limit1_down)/(steps1*1.0);
  
  for (Int_t i=0; i<steps1; i++) { // FIT loop
    Double_t Eb1test=limit1_down+i*stepsize1;
    
    m_fitrecord->UpdateEntry(HHEventRecord::b1)->SetEkeepM(Eb1test); // update 4-vectors with fit parameters
    ConstrainE2(HHEventRecord::hb, HHEventRecord::b1, HHEventRecord::b2);

    m_chi2_b1 = Chi2V4(HHEventRecord::b1);
    m_chi2_b2 = Chi2V4(HHEventRecord::b2);
    m_chi2 = m_chi2_b1 + m_chi2_b2; // chi2 calculation

    if (m_chi2 < chi2_min){
      chi2_min=m_chi2;
      minE1=Eb1test;
    }
    
    if (m_logLevel>1){
      std::cout << "chi2 b1: " << m_chi2_b1 << std::endl;
      std::cout << "chi2 b2: " << m_chi2_b2 << std::endl;
      std::cout << "------------------" << std::endl;
      std::cout << "chi2 : " << m_chi2 << std::endl;
      std::cout << "------------------" << std::endl;
    }
  }
  
  //calculate final best SetNextPoint
  m_fitrecord->UpdateEntry(HHEventRecord::b1)->SetEkeepM(minE1); // update 4-vectors with fit parameters
  ConstrainE2(HHEventRecord::hb, HHEventRecord::b1, HHEventRecord::b2);
  m_chi2_b1 = Chi2V4(HHEventRecord::b1);
  m_chi2_b2 = Chi2V4(HHEventRecord::b2);
  m_chi2 = m_chi2_b1 + m_chi2_b2; // chi2 calculation
   
  if (m_logLevel>0){
    std::cout << "chi2 b1: " << m_chi2_b1 << std::endl;
    std::cout << "chi2 b2: " << m_chi2_b2 << std::endl;
    std::cout << "------------------" << std::endl;
    std::cout << "chi2 : " << m_chi2 << std::endl;
    std::cout << "------------------" << std::endl;
  }
}

void
HHDiJetKinFit::Fit()
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
  Int_t icallNewton = -1;     //  init start of Newton Method
  Int_t iloop = 0;            // counter for falls to fit function

  // calculate htau from tauvis; recombine leaves measured entries in event record untouched
  m_recrecord->Recombine();

//  recrecord->GetEntry(V4EventRecord::tauvis1)->PrintEPxPyPzM(V4EventRecord::tauvis1);
//  recrecord->GetEntry(V4EventRecord::tau1)->PrintEPxPyPzM(V4EventRecord::tau1);
//  recrecord->GetEntry(V4EventRecord::tauvis2)->PrintEPxPyPzM(V4EventRecord::tauvis2);
//  recrecord->GetEntry(V4EventRecord::tau2)->PrintEPxPyPzM(V4EventRecord::tau2);
//  recrecord->GetEntry(V4EventRecord::htau)->PrintEPxPyPzM(V4EventRecord::htau);

  m_fitrecord->UpdateMothers(HHEventRecord::b1);

  // fill initial fit parameters
  astart[0] = m_fitrecord->GetEntry(HHEventRecord::b1)->E();         // energy of first b-jet
  aprec[0] = 0.01;   //0.1                 // precision for fit

  // fill initial step width
  h[0] = 0.5*m_fitrecord->GetEntry(HHEventRecord::b1)->dE();        // step width = bjet uncertainty

  daN[0] = 1.0;   //0.0                 // initial search direction in Eb-Etau diagonal

  // fit range
  alimit[0][0] = astart[0] - m_fitrecord->GetEntry(HHEventRecord::b1)->dE() * 5.; // b-jets: 5 sigma
  if (alimit[0][0]<0) alimit[0][0]=0;
  alimit[0][1] = astart[0] + m_fitrecord->GetEntry(HHEventRecord::b1)->dE() * 5.;

  for (Int_t ip = 0; ip < np; ip++) {
    a[ip] = astart[ip];
  }

  static const Int_t nloopmax = 100;
  static Double_t Xa[nloopmax];
  static Double_t Xa1[nloopmax];
  static Double_t HPx[nloopmax], HPy[nloopmax];
  static Double_t HPx1[nloopmax], HPy1[nloopmax];

  for (Int_t iloop = 0; iloop < m_maxloops * 10 && iter < m_maxloops; iloop++) { // FIT loop
    //std::cout << "iloop = " << iloop << std::endl;

    m_fitrecord->UpdateEntry(HHEventRecord::b1)->SetEkeepM(a[0]); // update 4-vectors with fit parameters
    //std::cout << "break-point 1 reached" << std::endl;
    //m_fitrecord->GetEntry(HHEventRecord::b1)->PrintEPxPyPzM(HHEventRecord::b1);
    //m_fitrecord->GetEntry(HHEventRecord::b2)->PrintEPxPyPzM(HHEventRecord::b2);
    ConstrainE2(HHEventRecord::hb, HHEventRecord::b1, HHEventRecord::b2);
    //std::cout << "break-point 2 reached" << std::endl;
    //m_fitrecord->GetEntry(HHEventRecord::b1)->PrintEPxPyPzM(HHEventRecord::b1);
    //m_fitrecord->GetEntry(HHEventRecord::b2)->PrintEPxPyPzM(HHEventRecord::b2);
    //std::cout << "break-point 3 reached" << std::endl;

    m_chi2_b1 = Chi2V4(HHEventRecord::b1);
    m_chi2_b2 = Chi2V4(HHEventRecord::b2);
    m_chi2 = m_chi2_b1 + m_chi2_b2; // chi2 calculation

    if (m_logLevel>1){
      std::cout << "chi2 b1: " << m_chi2_b1 << std::endl;
      std::cout << "chi2 b2: " << m_chi2_b2 << std::endl;
      std::cout << "------------------" << std::endl;
      std::cout << "chi2 : " << m_chi2 << std::endl;
      std::cout << "------------------" << std::endl;
    }

//    chi2 = Chi2V4(recrecord, fitrecord, V4EventRecord::b1)
//        + Chi2V4(recrecord, fitrecord, V4EventRecord::b2)
//        + Chi2PTmiss(*ptmissrec, fitrecord);  // chi2 calculation

    if (iloop >= 0 && iloop < nloopmax) {
      Xa1[iloop] = m_fitrecord->GetEntry(HHEventRecord::b1)->E();
    }
    if (iter >= 0 && iter < nloopmax) {
      Xa[iter] = m_fitrecord->GetEntry(HHEventRecord::b1)->E();
    }
    if (iloop >= 0 && iloop < nloopmax) {
      HPx1[iloop] = m_fitrecord->GetEntry(HHEventRecord::H)->Px();
      HPy1[iloop] = m_fitrecord->GetEntry(HHEventRecord::H)->Py();
    }
    if (iter >= 0 && iter < nloopmax) {
      HPx[iter] = m_fitrecord->GetEntry(HHEventRecord::H)->Px();
      HPy[iter] = m_fitrecord->GetEntry(HHEventRecord::H)->Py();
    }

    //    if (iloop>=0 && iloop<nloopmax) {Xloop[iloop] = a[0] ;  Yloop[iloop] = a[1] ; }
    //    if (iter>=0) { nIterations = Iterations->SetNextPoint( V4fit[V4EventRecord::b1].E(), V4fit[V4EventRecord::tau1].E() ) ;   }
    //    if (iter>=0) { nIterations = Iterations->SetNextPoint( a[0], a[1] ) ;   }

    PSMath::PSfitShow(iloop, m_convergence, iter, method, mode, m_printlevel,
                      m_graphicslevel, np, a, astart, alimit, aprec, daN, h,
                      m_chi2, g, H);

    if (m_convergence != 0) {
      break;
    }
    m_convergence = PSMath::PSfit(iloop, iter, method, mode, noNewtonShifts, m_printlevel,
                                  np, a, astart, alimit, aprec,
                                  daN, h, aMemory, m_chi2, chi2iter, g, H,
                                  Hinv);
  }

  if (m_logLevel>0){
    std::cout << "chi2 b1: " << m_chi2_b1 << std::endl;
    std::cout << "chi2 b2: " << m_chi2_b2 << std::endl;
    std::cout << "------------------" << std::endl;
    std::cout << "chi2 : " << m_chi2 << std::endl;
    std::cout << "------------------" << std::endl;
  }
}

void
HHDiJetKinFit::ConstrainE2(Int_t iv4, Int_t iv41, Int_t iv42)
{
  //std::cout << "<HHKinFit::ConstrainE2>:" << std::endl;

  m_fitrecord->UpdateMothers(iv41);

  Double_t M = m_fitrecord->GetEntry(iv4)->M();
  Double_t Mc = m_particlelist->GetParticleInfo( m_fitrecord->GetEntry(iv4)->ID() )->M();
  //std::cout << "Mc = " << Mc << ", M = " << M << std::endl;

  Double_t E1 = m_fitrecord->GetEntry(iv41)->E();
  Double_t Px1 = m_fitrecord->GetEntry(iv41)->Px();
  Double_t Py1 = m_fitrecord->GetEntry(iv41)->Py();
  Double_t Pz1 = m_fitrecord->GetEntry(iv41)->Pz();
  Double_t M1 = m_fitrecord->GetEntry(iv41)->M();
  //std::cout << "E1 = " << E1 << ", Px1 = " << Px1 << ", Py1 = " << Py1 << ", Pz1 = " << Pz1 << ", M1 = " << M1 << std::endl;

  Double_t E2 = m_fitrecord->GetEntry(iv42)->E();
  Double_t Px2 = m_fitrecord->GetEntry(iv42)->Px();
  Double_t Py2 = m_fitrecord->GetEntry(iv42)->Py();
  Double_t Pz2 = m_fitrecord->GetEntry(iv42)->Pz();
  Double_t M2 = m_fitrecord->GetEntry(iv42)->M();
  //std::cout << "E2 = " << E2 << ", Px2 = " << Px2 << ", Py2 = " << Py2 << ", Pz2 = " << Pz2 << ", M2 = " << M2 << std::endl;

  if ( M2 < (1.e-3*E2) ) { // massless case
    m_fitrecord->UpdateEntry(iv42)->SetEkeepM(E2 * (Mc / M) * (Mc / M)); // only changes absolute value and keeps eta, phi, m untouched
                                                                         // such that invariant mass of 1 and 2 gives correct mass Mc
    m_fitrecord->UpdateMothers(iv42);
    //std::cout << "Mnew = " << m_fitrecord->GetEntry(iv4)->M() << std::endl;
    return;
  }

  //Double_t p1_dot_p2 = E1*E2 - (Px1*Px2 + Py1*Py2 + Pz1*Pz2);
  Double_t p1_dot_p2 = 0.5*(M*M - (M1*M1 + M2*M2));
  Double_t sf = (TMath::Sqrt(p1_dot_p2*p1_dot_p2 + (Mc*Mc - M1*M1)*M2*M2) - p1_dot_p2)/(M2*M2);
  assert(sf >= 0.);
  Double_t E2new_sf = sf*E2;
  //std::cout << "sf = " << sf << " --> E2new = " << E2new_sf << std::endl;

  HHV4Vector* entry = m_fitrecord->GetEntry(iv42);
  m_fitrecord->UpdateEntry(iv42)->SetEEtaPhiM(E2new_sf, entry->Eta(), entry->Phi(), sf*M2);
  m_fitrecord->UpdateMothers(iv42);
  //std::cout << "Mnew = " << m_fitrecord->GetEntry(iv4)->M() << std::endl;
}

Double_t
HHDiJetKinFit::Chi2V4(Int_t iv4)
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
HHDiJetKinFit::SetPrintLevel(Int_t level){
  m_printlevel = level;
}

void
HHDiJetKinFit::SetGraphicsLevel(Int_t level){
  m_graphicslevel=level;
}

void
HHDiJetKinFit::SetMaxLoops(Int_t loops){
  m_maxloops = loops;
}

Double_t
HHDiJetKinFit::GetChi2()
{
  return m_chi2;
}

Double_t
HHDiJetKinFit::GetChi2_b1()
{
  return m_chi2_b1;
}

Double_t
HHDiJetKinFit::GetChi2_b2()
{
  return m_chi2_b2;
}

Double_t
HHDiJetKinFit::GetConvergence()
{
  return m_convergence;
}

HHEventRecord*
HHDiJetKinFit::GetFitrecord(){
  return m_fitrecord;
}

void
HHDiJetKinFit::SetLogLevel(Int_t level)
{
  m_logLevel = level;
}

double
HHDiJetKinFit::GetPullE(Int_t iv4){
  HHV4Vector* fv = m_fitrecord->GetEntry(iv4);
  HHV4Vector* rv = m_recrecord->GetEntry(iv4);
  return (fv->E()-rv->E())/rv->dE();
}

