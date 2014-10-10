#include "../include/BVNewtonStep.h"

#include<iostream>
#include <algorithm>
#include <cmath>

#include"TMatrixDEigen.h"

BVNewtonStep::BVNewtonStep(HHDiJetKinFitNewMini* fit, TVectorD start, TVectorD min, TVectorD max)
    : m_fit(fit),
      m_start(start),
      m_end(start),
      m_min(min),
      m_max(max),
      m_Npar(m_start.GetNoElements()),
      m_initialised(true),
      m_d(100000),
      m_dF(100000)
{

  for (int i=0; i<m_Npar; i++){
    if ((m_start[i]<m_min[i]) || (m_start[i]>m_max[i])) {
      std::cout << "ERROR: start value of parameter " << i << " is not within the given range" << std::endl;
      m_initialised=false;
    }
  }
  
  if ((m_start.GetNoElements()!=m_min.GetNoElements()) || (m_start.GetNoElements()!=m_max.GetNoElements())){
    std::cout << "ERROR: input vectors have different dimensions" << std::endl;
    m_initialised=false;
  }
  
  m_direction.ResizeTo(m_Npar);
  m_derivative.ResizeTo(m_Npar);
  m_hesse.ResizeTo(m_Npar,m_Npar);
  
  m_testedPoints.clear();
  
  std::cout << "#########################################" << std::endl;
  std::cout << "## NEWTON METHOD ########################" << std::endl;
  std::cout << "#########################################" << std::endl;
  std::cout << "calculating Newton step in " << m_Npar << "-dimensional parameter space" << std::endl;
  std::cout << "-----------------------------------------" << std::endl;
}

bool 
BVNewtonStep::run(){
  if(!m_initialised) return(false);
  
  ////////////////////////////////////////////
  // config parameters
  ////////////////////////////////////////////  
  double epsF = 0.01;   //arbitrary value
  double epsx = 0.004;  //arbitrary value
  double h = 0.1;       //arbitrary value, todo: the optimal h can be estimated (cf. Blobel ebook, p.199, (8.8))
    
  ////////////////////////////////////////////
  //define offset vector for making the shifts in the different parameters
  //the sign of the shift alternates
  ////////////////////////////////////////////
  std::vector<TVectorD> dx;
  for (int i=0; i<m_Npar; i++){
    TVectorD dxi(m_Npar);
    dxi.Zero();
    dxi[i]=((i%2==0)?h:-h);
//     dxi[i]=h;
    dx.push_back(dxi);
  }
  
  ////////////////////////////////////////////
  //calculate all the different ingredients for the numerical calculation of the gradient vector and the Hesse matrix
  //the code is organised like this in order to minimise the number of calls of function F
  ////////////////////////////////////////////
  double Fnoshift = m_fit->CalcChi2(m_start);
  
  std::vector<double> Fsingleshiftpos(m_Npar); Fsingleshiftpos.clear();
  std::vector<double> Fsingleshiftneg(m_Npar); Fsingleshiftneg.clear();
  
  //caculate first all function values for the gradient vector and the diagonal entries of the Hesse matrix
  for (int i=0; i<m_Npar; i++){
    Fsingleshiftpos[i]=m_fit->CalcChi2(m_start+dx[i]);
    Fsingleshiftneg[i]=m_fit->CalcChi2(m_start-dx[i]);
    m_testedPoints.push_back(m_start+dx[i]);
    m_testedPoints.push_back(m_start-dx[i]);
    
    m_derivative[i] = (Fsingleshiftpos[i]-Fsingleshiftneg[i])/(2*dx[i][i]);  // symmetrised derivative by up and down fluctuation ("zentrale Differenzenformel")
    m_hesse[i][i] = (Fsingleshiftpos[i] - 2*Fnoshift + Fsingleshiftneg[i])/(dx[i][i]*dx[i][i]);
  }
  
  //and now all function values for the off-diagonal entries of the Hesse matrix
  for (int i=0; i<m_Npar; i++){
    for (int j=i+1; j<m_Npar; j++){
      double Fdoubleshift=m_fit->CalcChi2(m_start+dx[i]+dx[j]);
      m_testedPoints.push_back(m_start+dx[i]+dx[j]);
      m_hesse[i][j]=(Fdoubleshift+Fnoshift-Fsingleshiftpos[i]-Fsingleshiftpos[j])/(dx[i][i]*dx[j][j]);
      m_hesse[j][i]=m_hesse[i][j];
    }
  }
  

  ////////////////////////////////////////////
  //check positive-definitness of hesse matrix
  ////////////////////////////////////////////
  TMatrixDEigen eigenmatrix(m_hesse);
  bool ispositivedefinite=true;
  double smallesteigenvalue=1000;
 
  for (int i=0; i<m_Npar; i++){
    double eigenvalue = eigenmatrix.GetEigenValues()(i,i);
    if (eigenvalue<0) {
      ispositivedefinite=false;
      if (eigenvalue<smallesteigenvalue){
        smallesteigenvalue=eigenvalue;
      }
    }
  }
  
  std::cout << "gradient vector: " << std::endl;
  m_derivative.Print();
  
  ////////////////////////////////////////////
  //correct Hesse matrix in order to make it positive definite
  //the larger the diagonal entries on the diagonalmatrix are, 
  //the more becomes the Newton direction the pure gradient
  ////////////////////////////////////////////
  if (!ispositivedefinite){
    TMatrixD diagonalmatrix(m_Npar,m_Npar);
    for (int i=0; i<m_Npar; i++){
      diagonalmatrix[i][i]=-smallesteigenvalue;
    }
    std::cout << "Hesse matrix is NOT positive-definite: " << std::endl;
    m_hesse.Print();
    m_hesse += diagonalmatrix;
    std::cout << "WARNING: Hesse matrix was modified in order to make it positive-definite" <<std::endl;
    m_hesse.Print();    
  }
  else{
    std::cout << "Hesse matrix is positive-definite: " << std::endl;
    m_hesse.Print();    
  }
  
  
  ////////////////////////////////////////////
  // calculate Newton step: H.direction = -derivative
  ////////////////////////////////////////////
  TMatrixD invhesse = m_hesse.Invert();
  m_direction=m_derivative;
  m_direction*=invhesse;
  m_direction*=-1;
  


  ////////////////////////////////////////////
  //calculate convergence quantities
  ////////////////////////////////////////////
  double Fnew = m_fit->CalcChi2(m_start+m_direction);
  m_dF = Fnoshift - Fnew;
  m_d = m_direction*m_derivative;
  m_d *= -1;

  std::cout << "#########################################" << std::endl;
  std::cout << "convergence parameters: dF=" << getdF() << " d=" << getd() << std::endl;
  std::cout << "found Newton direction: " << std::endl;
  m_direction.Print();
  std::cout << "#########################################" << std::endl;
  
  return(true);
}


TGraph2D 
BVNewtonStep::getGraphF(int resolution){
  TGraph2D gr((resolution+1)*(resolution+1));
  gr.SetName("F");
  double stepx0=(m_max[0]-m_min[0])/resolution;
  double stepx1=(m_max[1]-m_min[1])/resolution;
  
  for (int i=0; i<=resolution; i++){
    for (int j=0; j<=resolution; j++){
      double x0=m_min[0]+i*stepx0;
      double x1=m_min[1]+j*stepx1;
      TVectorD x(2); x[0]=x0; x[1]=x1;
      double F=m_fit->CalcChi2(x);
      gr.SetPoint(i*(resolution+1)+j,x0,x1,F);
    }
  }
  return(gr);
}


TGraph 
BVNewtonStep::getGraphStartPoint(){
  TGraph gr(1);
  gr.SetName("StartPoint");
  gr.SetPoint(0,m_start[0],m_start[1]);
  gr.SetMarkerSize(2);
  return(gr);
}

TGraph 
BVNewtonStep::getGraphDirectionDerivative(){
  TGraph gr(2);
  gr.SetName("DirectionDerivative");
  gr.SetPoint(0,m_start[0],m_start[1]);
  gr.SetPoint(1,m_start[0]+(m_derivative[0]/sqrt(m_derivative.Norm2Sqr())),m_start[1]+(m_derivative[1]/sqrt(m_derivative.Norm2Sqr())));
  return(gr);
}

TGraph 
BVNewtonStep::getGraphDirectionNewtonStep(){
  TGraph gr(2);
  gr.SetName("DirectionNewtonStep");
  gr.SetPoint(0,m_start[0],m_start[1]);
  gr.SetPoint(1,m_start[0]+(m_direction[0]/sqrt(m_direction.Norm2Sqr())),m_start[1]+(m_direction[1]/sqrt(m_direction.Norm2Sqr())));
  return(gr);
}

TGraph 
BVNewtonStep::getGraphTestedPoints(){
  TGraph gr(m_testedPoints.size());
  gr.SetName("TestedPoints");
  for (int i=0; i<m_testedPoints.size(); i++){
    gr.SetPoint(i,m_testedPoints[i][0],m_testedPoints[i][1]);
  }
  return(gr);
}
