/*
 * HHParticleList.h
 *
 *  Created on: 16.06.2014
 */

#ifndef HHPARTICLELIST_H_
#define HHPARTICLELIST_H_

#include <Rtypes.h>
#include <TString.h>

#ifdef HHKINFIT_STANDALONE
#include "../interface/HHPID.h"
#include "../interface/HHParticle.h"
#include "../interface/HHV4Vector.h"
#else
#include "HHKinFit/HHKinFit/interface/HHPID.h"
#include "HHKinFit/HHKinFit/interface/HHParticle.h"
#include "HHKinFit/HHKinFit/interface/HHV4Vector.h"
#endif

class HHParticleList
{
public:
  HHParticleList();

  HHParticle*  GetParticleInfo(HHPID::pid id);
  HHV4Vector*  GenerateV4(HHPID::pid id);
  void  UpdateMass(HHPID::pid id, Double_t m);
  void  Print(TString name = "");

private:
  // largest number for array dimensions
  const static Int_t m_v4maxPartList = 22;
  HHParticle m_particlelist[m_v4maxPartList];
};

#endif /* HHPARTICLELIST_H_ */
