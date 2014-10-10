#ifndef BVNEWTONSTEP_H_
#define BVNEWTONSTEP_H_

#include "../include/HHDiJetKinFitNewMini.h"
#include<vector>

#include"TVectorD.h"
#include"TMatrixD.h"
#include "TGraph2D.h"
#include "TGraph.h"


class BVNewtonStep{
public:
  BVNewtonStep(HHDiJetKinFitNewMini* fit, TVectorD start, TVectorD min, TVectorD max);
  bool run();

  TVectorD getStartPoint()  {return(m_start);};
  TVectorD getEndPoint()    {return(m_start+m_direction);};
  TVectorD getDirection()   {return(m_direction);};

  double getdF()            {return(m_dF);}; 
  double getd()             {return(m_d);};
  
  TGraph2D getGraphF(int resolution=30);
  TGraph   getGraphStartPoint();
  TGraph   getGraphDirectionDerivative();
  TGraph   getGraphDirectionNewtonStep();
  TGraph   getGraphTestedPoints();
  
private:
  HHDiJetKinFitNewMini* m_fit;
  
  TVectorD m_start;
  TVectorD m_end;
  TVectorD m_direction;  // contains later the result
  TVectorD m_derivative;
  TMatrixD m_hesse;
  
  TVectorD m_min;
  TVectorD m_max;
  
  int m_Npar;
  bool m_initialised;

  double m_d;
  double m_dF;

  std::vector<TVectorD> m_testedPoints;
};

#endif