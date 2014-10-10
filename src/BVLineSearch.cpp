#include "../include/BVLineSearch.h"

#include<vector>
#include<iostream>
#include <algorithm>
#include <cmath>


struct absvalsort {
  bool operator() (double i, double j) { return (fabs(i)<fabs(j));}
};

BVLineSearch::BVLineSearch(HHDiJetKinFitNewMini* fit, TVectorD start, TVectorD direction, TVectorD min, TVectorD max)
    : m_fit(fit),
      m_start(start),
      m_end(start),
      m_direction(direction),
      m_min(min),
      m_max(max),
      m_Npar(m_start.GetNoElements()),
      m_initialised(true)
{

  for (int i=0; i<m_Npar; i++){
    if ((m_start[i]<m_min[i]) || (m_start[i]>m_max[i])) {
      std::cout << "ERROR: start value of parameter " << i << " is not within the given range" << std::endl;
      m_initialised=false;
    }
  }
  
  if ((m_start.GetNoElements()!=m_direction.GetNoElements()) || (m_start.GetNoElements()!=m_min.GetNoElements()) || (m_start.GetNoElements()!=m_max.GetNoElements())){
    std::cout << "ERROR: input vectors have different dimensions" << std::endl;
    m_initialised=false;
  }
  
  std::cout << "#########################################" << std::endl;
  std::cout << "## 1D LINE SEARCH #######################" << std::endl;
  std::cout << "#########################################" << std::endl;
  
}


TVectorD
BVLineSearch::getVector(double x){
  return( m_start + (x*m_direction) );
}

TVectorD
BVLineSearch::calcLimits(){
  std::vector<double> lambda; lambda.clear(); //stores all crossing points with axes
  for (int i=0; i<m_Npar; i++){
    lambda.push_back((m_min[i]-m_start[i])/m_direction[i]);
    lambda.push_back((m_max[i]-m_start[i])/m_direction[i]);
  }
  
  std::sort (lambda.begin(), lambda.end(), absvalsort());
  
  TVectorD result(2); //result vector
  result[0]=lambda[0];
  if (lambda.size()>2) // if more than 2 limits (parameter space > 1) some degeneracy of two limits could appear
    result[1]=(fabs(lambda[1]-lambda[0])<0.00001?lambda[2]:lambda[1]);
  else
    result[1]=lambda[1];
  
  if (result[0]>result[1]){
    double temp = result[0];
    result[0]=result[1];
    result[1]=temp;
  }
  
  return(result);
}

bool 
BVLineSearch::run(){
  if(!m_initialised) return(false);
  
  TVectorD limit=calcLimits();
  double xmin=limit[0];
  double xmax=limit[1];

  for (int i=0; i<m_Npar; i++){
    std::cout << "min" <<i<< ": " << m_min[i] << " "<< "max" <<i<< ": " << m_max[i] << " "<< "start" <<i<< ": " << m_start[i] << " " << "direction" <<i<< ": "<< m_direction[i] << std::endl;
  }
  std::cout << "limits: (" << xmin << ";" << xmax << ")" << std::endl;
  
  ////////////////////////////////////////////
  // config parameters
  ////////////////////////////////////////////  
  double tau = 0.618034; // Golden Ratio
  double alpha = 1;      // step weighting for interval search - for golden ratio user alpha= tau - tau*tau
  double epsF = 0.01;    // arbitrary value
  double epsx = 0.004;   // arbitrary value
  double minsteps = 10;  // by default the norm of the m_direction is used as initial step width, but...
  double stepx = std::min(sqrt(m_direction.Norm2Sqr()),(xmax-xmin)/minsteps); //... if this is too large devide the whole scan interval in minsteps steps
  
  std::cout << "initial step size: " << stepx << std::endl;
  std::cout << "-----------------------------------------" << std::endl;

  ////////////////////////////////////////////
  // 1.) search for minimum interval
  ////////////////////////////////////////////
  std::vector<double> xsearch;   xsearch.clear();  //xvalues
  std::vector<double> Fsearch;   Fsearch.clear();  //functionvalues
  
  //calculate two start values
  xsearch.push_back(0);        Fsearch.push_back(m_fit->CalcChi2(getVector(xsearch[0])));
  xsearch.push_back(0+stepx);  Fsearch.push_back(m_fit->CalcChi2(getVector(xsearch[1])));
  
  //and sort them in right order
  if (Fsearch[0]<Fsearch[1]) { //sort x values such that Fsearch[x0]>Fsearch[x1]
    int temp=0;
    temp = xsearch[0];
    xsearch[0] = xsearch[1];
    xsearch[1] = temp;
    temp = Fsearch[0];
    Fsearch[0] = Fsearch[1];
    Fsearch[1] = temp;
  }
 
  //calculate iteratively new points and check if between three points a minimum can be localised
  bool foundminimuminterval=false;
  bool reachedlimit=false;
  std::vector<double> x(3); x.clear(); //xvaluesinterval
  std::vector<double> F(3); F.clear(); //functionvaluesinterval
  for (int i = 2; !foundminimuminterval ;i++){
    double xnew = xsearch[i-1] + alpha*(xsearch[i-1]-xsearch[i-2]); // the way how new x values are generated
    if (xnew<xmin) { // if xnew crosses the limits, set it to the crossed limit
      xnew=xmin;
      reachedlimit=true; // if this is true the series cannot go on
      std::cout << "WARNING: touched the limit" << std::endl;
    }
    if (xnew>xmax){
      xnew=xmax;
      reachedlimit=true;
      std::cout << "WARNING: touched the limit" << std::endl;
    }
    xsearch.push_back(xnew);
    Fsearch.push_back(m_fit->CalcChi2(getVector(xsearch[i])));
    std::cout << "checked: x    = (" << xsearch[i-2] << ";" << xsearch[i-1] << ";" << xsearch[i] << ")" <<std::endl;
    std::cout << "checked: F(x) = (" << Fsearch[i-2] << ";" << Fsearch[i-1] << ";" << Fsearch[i] << ")" <<std::endl;
    if (Fsearch[i]>Fsearch[i-1]) {
      foundminimuminterval=true;
      if (xsearch[i-2]>xsearch[i]) {      //sort x values such that x[0]<x[1]<x[2] - this sorting is needed for 2.)
        x[2]=xsearch[i-2];  F[2]=Fsearch[i-2];
        x[1]=xsearch[i-1];  F[1]=Fsearch[i-1];
        x[0]=xsearch[i];    F[0]=Fsearch[i];
      }
      else {
        x[0]=xsearch[i-2];  F[0]=Fsearch[i-2];
        x[1]=xsearch[i-1];  F[1]=Fsearch[i-1];
        x[2]=xsearch[i];    F[2]=Fsearch[i];
      }
    }
    if (foundminimuminterval)
      std::cout << "  ->found minimum" << std::endl;
    else{
       std::cout << "  ->no minimum" << std::endl;
    }
    
    if (!foundminimuminterval && reachedlimit){  //if in an iteration where a limit was reached no minimum has been found, the whole method stops and returns false
      std::cout << "ERROR: found no minimum within the limits" << std::endl;
      return(foundminimuminterval);
    }
    std::cout << "-----------------------------------------" << std::endl;
  }

  if (!foundminimuminterval) return(foundminimuminterval);
  
  ////////////////////////////////////////////
  // 2.) reduce minimum interval
  ////////////////////////////////////////////
  
  bool convergence = false;
  int iterations = 0;
  std::vector<double> xresult(3); xresult.clear(); //xresult
  std::vector<double> Fresult(3); Fresult.clear(); //functionvaluesresult
  while (!convergence){
    //todo: in this while loop, one could also additionally use a maximum number of iterations as break condition
    //however, if a minimum interval has been found, convergence should always be possible for smooth functions F
    //a non-convergence is a hint for bugs
    
    // monitoring of the shrinking interval
    std::cout << "interval after " << iterations << " minimisation iterations:" << std::endl;
    std::cout << "x0=" << x[0] << "    " << "F[x0]=" << F[0] << std::endl;
    std::cout << "x1=" << x[1] << "    " << "F[x1]=" << F[1] << std::endl;
    std::cout << "x2=" << x[2] << "    " << "F[x2]=" << F[2] << std::endl;
  
    // sanity check of the interval
    if (!(F[0]>F[1] && F[2]>F[1] && x[0]<x[1] && x[1]<x[2])){
      std::cout << "ERROR: found interval is wrong" << std::endl;
      return(false);
    }
    
    // iterative minimisation
    //try first quadratic interpolation (\approx Newton method in 1D)
    double dF  = 1.0/(x[2]-x[0])*((F[2]-F[1])*(x[1]-x[0])/(x[2]-x[1]) + (F[1]-F[0])*(x[2]-x[1])/(x[1]-x[0]));
    double ddF = 2.0/(x[2]-x[0])*((F[2]-F[1])/(x[2]-x[1]) - (F[1]-F[0])/(x[1]-x[0]));
    double xt  = x[1]-dF/ddF;
  
    // if the found test point is too close to x[1] use the Golden Ratio
    if (fabs(dF/ddF) < (x[2]-x[0])*0.1){ // xt is very close to x[1] (10% of total interval)
      xt = x[0] + tau * (x[2]-x[0]);     // in this case use better Golden Ratio to minimise the interval
      std::cout << "using Golden Ratio:" << std::endl;
    }
    else{
      std::cout << "using quadratic interpolation:" << std::endl;
    }
    
    double Ft = m_fit->CalcChi2(getVector(xt));
    std::cout << "xt=" << xt << "    " << "F[xt]=" << Ft << std::endl;
    
    // calculate the change of the minimum point
    double dMinx = xt-x[1];
    double dMinF = Ft-F[1];
    std::cout << "dx=" << dMinx << std::endl;
    std::cout << "dF=" << dMinF << std::endl;
    
    //test for convergence
    if(fabs(dMinx)<epsx || fabs(dMinF)<epsF) {
      convergence=true;
      xresult[0]=x[0]; Fresult[0]=F[0];
      xresult[1]=x[1]; Fresult[1]=F[1];
      xresult[2]=x[2]; Fresult[2]=F[2];
      break;
    }
    
    //if not converged yet, reduce the interval further depending on the value of Ft
    if (Ft<F[1]){ // xt becomes the new x1 and x1 becomes the new x3
      if (xt<x[1]) x[2]=x[1]; F[2]=F[1];
      if (xt>x[1]) x[0]=x[1]; F[0]=F[1];
      x[1]=xt;   F[1]=Ft;
    }   
    else{
      if (xt<x[1]) x[0]=xt;  F[0]=Ft;
      if (xt>x[1]) x[2]=xt;  F[2]=Ft;
    }
    std::cout << "-----------------------------------------" << std::endl;
    iterations++; 
  }
  
  std::cout << "#########################################" << std::endl;
  std::cout << "Found a good minimum after " << iterations << " iterations:" << std::endl;
  std::cout << "xmin=" << xresult[1] << " (" << xresult[0] << ";" << xresult[2] << ")" << std::endl;
  std::cout << "Fmin=" << Fresult[1] << " (" << Fresult[0] << ";" << Fresult[2] << ")" << std::endl;
  std::cout << "The corresponding n-dimensional parameter vector becomes:" << std::endl;
  m_end=getVector(xresult[1]);
  m_end.Print();
  std::cout << "#########################################" << std::endl;
  
  return(convergence);
}


TGraph2D 
BVLineSearch::getGraphF(int resolution){
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
BVLineSearch::getGraphLine(int resolution){
  TGraph gr(resolution+1);
  gr.SetName("Line");
  TVectorD limit=calcLimits();
  double xmin=limit[0];
  double xmax=limit[1];
  double xstep = (xmax-xmin)/resolution;
  
  for (int i=0; i<=resolution; i++){
    TVectorD x = getVector(xmin+i*xstep);
    gr.SetPoint(i,x[0],x[1]);
  }
  return(gr);
}



TGraph 
BVLineSearch::getGraphFOnLine(int resolution){
  TGraph gr(resolution+1);
  gr.SetName("FOnLine");
  TVectorD limit=calcLimits();
  double xmin=limit[0];
  double xmax=limit[1];
  double xstep = (xmax-xmin)/resolution;
  
  for (int i=0; i<=resolution; i++){
    TVectorD x = getVector(xmin+i*xstep);
    double F=m_fit->CalcChi2(x);
    gr.SetPoint(i,xmin+i*xstep,F);
  }
  return(gr);
  
}

TGraph 
BVLineSearch::getGraphStartPoint(){
  TGraph gr(1);
  gr.SetName("StartPoint");
  gr.SetPoint(0,m_start[0],m_start[1]);
  gr.SetMarkerSize(2);
  return(gr);
}


TGraph 
BVLineSearch::getGraphEndPoint(){
  TGraph gr(1);
  gr.SetName("EndPoint");
  gr.SetPoint(0,m_end[0],m_end[1]);
  gr.SetMarkerSize(2);
  return(gr);
}

TGraph 
BVLineSearch::getGraphEndPointOnLine(){
  TGraph gr(1);
  gr.SetName("EndPointOnLine");
  double F=m_fit->CalcChi2(m_end);
  double x= (m_end[0]-m_start[0])/m_direction[0];  
  gr.SetPoint(0,x,F);
  gr.SetMarkerSize(2);
  return(gr);
}

TGraph 
BVLineSearch::getGraphStartPointOnLine(){
  TGraph gr(1);
  gr.SetName("StartPointOnLine");
  double F=m_fit->CalcChi2(m_start);
  gr.SetPoint(0,0,F);
  gr.SetMarkerSize(2);
  return(gr);
}
