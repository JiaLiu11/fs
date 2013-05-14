#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cmath>
#include "arsenal.h"
#include "stdlib.h"

static double T00 =     0.01359279277226;
static double T01 =    0.008201951912154;
static double T02 =     -0.0081977832369;
static double T11 =    0.006806438809342;
static double T12 =   -0.003553110944959;
static double T22 =    0.006786353962921;


using namespace std;

double func(double phiv)
{
    double A,B;

  	double func1;  


    A = T01*cos(phiv)+T02*sin(phiv);
    B = T00+T11*cos(phiv)*cos(phiv)
       +T22*sin(phiv)*sin(phiv)
       +T12*sin(2.0*phiv);

    if((B*B-4.0*A*A)<0)
    {
      cout<<"angle="<<phiv<<"; sqrt negative"<<endl;
      exit(-1);
    }
    
    func1 = sin(2.0*phiv)/(2.0*A)*(B*B-2*A*A-B*sqrt(B*B-4.0*A*A))
           -1/(2.0*A)*((T00+0.5*T11+0.5*T22)*sin(2.0*phiv)+T12)*(B-sqrt(B*B-4.0*A*A))            
           +T01*sin(phiv)+T02*cos(phiv);



    return func1;

}

double funcDev(double phi)
{

  double result;  


  result=T01*cos(phi) - T02*sin(phi) - ((T12 + (T00 + T11/2. + T22/2.)*sin(2*phi))*
      (2*T12*cos(2*phi) - 2*T11*cos(phi)*sin(phi) + 2*T22*cos(phi)*sin(phi) - 
        (-8*(T02*cos(phi) - T01*sin(phi))*(T01*cos(phi) + T02*sin(phi)) + 
           2*(2*T12*cos(2*phi) - 2*T11*cos(phi)*sin(phi) + 2*T22*cos(phi)*sin(phi))*(T00 + T11*pow(cos(phi),2) + T22*pow(sin(phi),2) + T12*sin(2*phi)))
          /(2.*sqrt(-4*pow(T01*cos(phi) + T02*sin(phi),2) + pow(T00 + T11*pow(cos(phi),2) + T22*pow(sin(phi),2) + T12*sin(2*phi),2)))))/
    (2.*(T01*cos(phi) + T02*sin(phi))) - ((T00 + T11/2. + T22/2.)*cos(2*phi)*
      (T00 + T11*pow(cos(phi),2) + T22*pow(sin(phi),2) + T12*sin(2*phi) - 
        sqrt(-4*pow(T01*cos(phi) + T02*sin(phi),2) + pow(T00 + T11*pow(cos(phi),2) + T22*pow(sin(phi),2) + T12*sin(2*phi),2))))/
    (T01*cos(phi) + T02*sin(phi)) + ((T02*cos(phi) - T01*sin(phi))*(T12 + (T00 + T11/2. + T22/2.)*sin(2*phi))*
      (T00 + T11*pow(cos(phi),2) + T22*pow(sin(phi),2) + T12*sin(2*phi) - 
        sqrt(-4*pow(T01*cos(phi) + T02*sin(phi),2) + pow(T00 + T11*pow(cos(phi),2) + T22*pow(sin(phi),2) + T12*sin(2*phi),2))))/
    (2.*pow(T01*cos(phi) + T02*sin(phi),2)) + (sin(2*phi)*(-4*(T02*cos(phi) - T01*sin(phi))*(T01*cos(phi) + T02*sin(phi)) + 
        2*(2*T12*cos(2*phi) - 2*T11*cos(phi)*sin(phi) + 2*T22*cos(phi)*sin(phi))*(T00 + T11*pow(cos(phi),2) + T22*pow(sin(phi),2) + T12*sin(2*phi)) - 
        ((T00 + T11*pow(cos(phi),2) + T22*pow(sin(phi),2) + T12*sin(2*phi))*
           (-8*(T02*cos(phi) - T01*sin(phi))*(T01*cos(phi) + T02*sin(phi)) + 
             2*(2*T12*cos(2*phi) - 2*T11*cos(phi)*sin(phi) + 2*T22*cos(phi)*sin(phi))*
              (T00 + T11*pow(cos(phi),2) + T22*pow(sin(phi),2) + T12*sin(2*phi))))/
         (2.*sqrt(-4*pow(T01*cos(phi) + T02*sin(phi),2) + pow(T00 + T11*pow(cos(phi),2) + T22*pow(sin(phi),2) + T12*sin(2*phi),2))) - 
        (2*T12*cos(2*phi) - 2*T11*cos(phi)*sin(phi) + 2*T22*cos(phi)*sin(phi))*
         sqrt(-4*pow(T01*cos(phi) + T02*sin(phi),2) + pow(T00 + T11*pow(cos(phi),2) + T22*pow(sin(phi),2) + T12*sin(2*phi),2))))/
    (2.*(T01*cos(phi) + T02*sin(phi))) + (cos(2*phi)*(-2*pow(T01*cos(phi) + T02*sin(phi),2) + 
        pow(T00 + T11*pow(cos(phi),2) + T22*pow(sin(phi),2) + T12*sin(2*phi),2) - 
        (T00 + T11*pow(cos(phi),2) + T22*pow(sin(phi),2) + T12*sin(2*phi))*
         sqrt(-4*pow(T01*cos(phi) + T02*sin(phi),2) + pow(T00 + T11*pow(cos(phi),2) + T22*pow(sin(phi),2) + T12*sin(2*phi),2))))/
    (T01*cos(phi) + T02*sin(phi)) - ((T02*cos(phi) - T01*sin(phi))*sin(2*phi)*
      (-2*pow(T01*cos(phi) + T02*sin(phi),2) + pow(T00 + T11*pow(cos(phi),2) + T22*pow(sin(phi),2) + T12*sin(2*phi),2) - 
        (T00 + T11*pow(cos(phi),2) + T22*pow(sin(phi),2) + T12*sin(2*phi))*
         sqrt(-4*pow(T01*cos(phi) + T02*sin(phi),2) + pow(T00 + T11*pow(cos(phi),2) + T22*pow(sin(phi),2) + T12*sin(2*phi),2))))/
    (2.*pow(T01*cos(phi) + T02*sin(phi),2));
}

double func2(double phiv)
{


    return cos(phiv);

}

double BinaryRootFinder(double (*func)(double), double value, double xL, double xR, double precision)
// Binary search for a root func(x)=value, i.e. find root for func(x)-value=0
// 
{   double x1=xL;
    double x2=xR;
    double xx=(x1+x2)/2.0;
    double f1,f2,f0;   //f(x)-value at the left ending, right ending and mid-point of the region

    while(abs((*func)(xx)-value)>precision)
    {
      xx=(x1+x2)/2.0;
      f1=(*func)(x1)-value;
      f2=(*func)(x2)-value;
      f0=(*func)(xx)-value;

      if(f1*f0<0)
        x2=xx;
      else if(f0*f2<0)
        x1=xx;
      else if(f1*f0>0&&f0*f2>0)
      {
        cout<<"no root in this region"<<endl;
        exit (-1);
      }
    }
    return xx;
}

double CalVelocity(double phiv)
{
    double A,B;

    A = T01*cos(phiv)+T02*sin(phiv);
    B = T00+T11*cos(phiv)*cos(phiv)
       +T22*sin(phiv)*sin(phiv)
       +T12*sin(2.0*phiv);

    if((B*B-4.0*A*A)<0)
    {
      cout<<"angle="<<phiv<<"; sqrt negative"<<endl;
      exit(-1);
    }
    return (B-sqrt(B*B-4.0*A*A))/(2.0*A);
}

bool MatchingCheck(double v, double phiv, double tolerance)
{
  double ed, u0, u1, u2;
  double eq0_LHS, eq0_RHS, eq1_LHS, eq1_RHS, eq2_LHS, eq2_RHS;
  bool eq1=false, eq2=false, eq0=false;


  ed=T00-T01*v*cos(phiv)-T02*v*sin(phiv);

//safty check
  if(1-v*v*cos(phiv)*cos(phiv)-v*v*sin(phiv)*sin(phiv)<0)
    exit(-1);

  double gamma= 1/sqrt(1-v*v*cos(phiv)*cos(phiv)-v*v*sin(phiv)*sin(phiv));

  u0=gamma;
  u1=gamma*v*cos(phiv);
  u2=gamma*v*sin(phiv);

// Matching Condition:
  eq0_LHS=T00*u0-T01*u1-T02*u2;
  eq0_RHS=ed*u0;
  eq1_LHS=T01*u0-T11*u1-T12*u2;
  eq1_RHS=ed*u1;
  eq2_LHS=T02*u0-T12*u1-T22*u2;
  eq2_RHS=ed*u2;


//debug 
  cout<<"root is: phi="<<phiv<<endl;
  cout<<"energy density ="<<ed<<" velocity="<<v<<endl;
  cout<<"Eq0, LHS= "<<eq0_LHS<<"; RHS="<<eq0_RHS<<endl;
  cout<<"Eq1, LHS= "<<eq1_LHS<<"; RHS="<<eq1_RHS<<endl;
  cout<<"Eq2, LHS= "<<eq2_LHS<<"; RHS="<<eq2_RHS<<endl<<endl;


  if(abs((eq0_LHS-eq0_RHS)/(eq0_LHS+eq0_RHS))<tolerance)
    eq0=true;
  else return 0;

  if(abs((eq1_LHS-eq1_RHS)/(eq1_LHS+eq1_RHS))<tolerance)
    eq1=true;
  else return 0;

  if(abs((eq2_LHS-eq2_RHS)/(eq2_LHS+eq2_RHS))<tolerance)
    {eq2=true; }
  else return 0;

  return eq0*eq1*eq2;

}

bool RootFilter(double root)
//return "true" if the root is in accordance with the rule
{   
  double vel=CalVelocity(root);
  bool eq;

  eq=MatchingCheck(vel, root, 1e-2);

  if(vel>=0&&eq==true)
    return true;
  else return false;
}



int MultiRootFinder(double (*func)(double), double value, double *roots, int division, double xL, double xR, double precision)
// find multiple roots using binary search method. Divide the region into several smaller ones
{
  int tolerance=100;
  int idx=0, idx_roots=0;
  double x1,x2;
  double ix=(xR-xL)/division;

  for(int i=0;i<division;i++)
  {
    x1=xL+i*ix;
    x2=x1+ix;
    double f1=(*func)(x1)-value;
    double f2=(*func)(x2)-value;

    if(f1*f2<0)    //there is at least one root in [x1,x2]
    {
      roots[idx]=BinaryRootFinder(func,value,x1,x2,precision);
      idx++;
    }

    else if(f1*f2>0)   //there is no root within [x1,x2], or two root within this region
    {
      for(int j=0;j<tolerance;j++)
      {
        double xx=drand(x1,x2);
        double fx=(*func)(xx)-value;
        if(f1*fx<=0)
        { 
          roots[idx]=BinaryRootFinder(func,value,x1,xx,precision);
          idx++;
          roots[idx]=BinaryRootFinder(func,value,xx,x2,precision);
          idx++;
//if roots are found, this loop should not continue, which causes double counting          
          break;   
        }
      }
    }
  }//division loop

  if(idx>0)
    {
      return idx;    
    }
    
  else
  {
    cout<<"Multifinder: no roots in this region"<<endl;
    return 0;
  }
}

int EhMultiRootFinder(double (*func)(double), double (*funcDev)(double), double value, 
  double *roots, int division, double xL, double xR, double precision, bool filter)
//this function is designed to find the points of a curve which are tangential to the x-axis. This kind
//of roots cannot be found by a regular binary search routine, which needs the function values at 
//boundary be a positive one and a negative one.
//1. f(x)=(*func)(), the function needs to be invert
//2. f'(x)=(*funcDev)(), the first derivative of f(x)
//3. y(x0)=f(x0)-value==0, this routine finds x0;
//4. *roots, the x0 will be stored here;
//5. division: specify the number of sub-regions, make sure the function is monotonous in each sub-region
//6. xL and xR are the left and right boundary of the function;
//7. precision: the precision of the binary search;
//8. boolean variable filter: decide if the roots need to be in accordance with certain rule. rule must 
//  be contained in function: bool RootFilter(double);
{
  double temp[20];
  double temp1[20];
  double fx;
  int idx_roots=0;
  int idx=MultiRootFinder(func, value, temp, division, xL, xR, precision);
  int idx1=MultiRootFinder(funcDev, value, temp1, division, xL, xR, precision);
  
//test the zero points in 
  for(int i=0;i<idx1;i++)
  {
    fx=(*func)(temp1[i])-value;
    cout<<" "<<fx<<"precision is "<<precision<<endl;
    if(abs(fx)<=precision)
    {
      temp[idx]=temp1[i];
      idx++;
    }
    else continue;
  }
//test if the roots are in accordance with the rule

  if(filter==true)
  {
    for(int i=0;i<idx;i++)
    { 

      if(RootFilter(temp[i])==true)
      {
        roots[idx_roots]=temp[i];
        idx_roots++;
      }
      else continue;
    }
  }

  else
  {
    for(int i=0;i<idx;i++)
    {
        roots[idx_roots]=temp[i];
        idx_roots++;
    }
  }

  if(idx_roots>0)
  {
// debug
    cout<<"There are "<<idx_roots<<" "<<"roots in this region"<<endl;   
    for(int i=0;i<idx_roots;i++)
    {
      cout<<"root "<<i+1<<": "<<setw(24)<<setprecision(8)<<roots[i]<<"   f(x)="<<CalVelocity(roots[i])<<endl;
 
    } 
    return idx_roots;    
  }
    
  else
  {
    cout<<"no roots in this region"<<endl;
    return 0;
  }
    
}

int main()
{
	  double phiv;
    double x1,x2,xp;
    int sol;
    double roots[20];

    x1=0;
    x2=2*M_PI;
     sol=EhMultiRootFinder(func, funcDev,0.0, roots, 10, x1,x2,1e-8,true);
    // sol=MultiRootFinder(func,0.0, roots, 10, x1,x2,1e-9);

    return 0;
}



