#ifndef BVLINESEARCH_H_
#define BVLINESEARCH_H_

#include "../include/HHDiJetKinFitNewMini.h"
#include "TVectorD.h"
#include "TGraph2D.h"
#include "TGraph.h"

class BVLineSearch{
public:
  BVLineSearch(HHDiJetKinFitNewMini* fit, TVectorD start, TVectorD direction, TVectorD min, TVectorD max);
  bool run();
   
  TVectorD getStartPoint()  {return(m_start);};
  TVectorD getEndPoint()    {return(m_end);};
  TVectorD getDirection()   {return(m_direction);};
  
  TGraph2D getGraphF(int resolution=30);
  TGraph   getGraphLine(int resolution=30);
  TGraph   getGraphFOnLine(int resolution=30);
  TGraph   getGraphStartPoint();
  TGraph   getGraphEndPoint();
  TGraph   getGraphEndPointOnLine();
  TGraph   getGraphStartPointOnLine();
  
private:
  HHDiJetKinFitNewMini* m_fit;
  
  TVectorD m_start;
  TVectorD m_end;  // contains later the result
  TVectorD m_direction;
  
  TVectorD m_min;
  TVectorD m_max;
  
  int m_Npar;
  bool m_initialised;
  
  TVectorD getVector(double x);
  TVectorD calcLimits();
};

#endif
