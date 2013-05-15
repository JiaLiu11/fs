#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include "LdMatching.h"
#include "Freestreaming.h"
#include "gauss_quadrature.h"

using namespace std;

int main()
{
// LdMatching(double xmax, double ymax, double ptmin, double ptmax, double dx0,double dy0,
// 		    double dpt, int nrap, double rmin, double rmax, string filename);

//----------------resume after debug    
  double tau_min=0.0, dtau=0.6;
  double tau_max=0.6;
  int nevents=20;

  //processing events
  for(int event_num=1;event_num<=nevents;event_num++)
  {
    //prepare readin filename
    ostringstream filename_stream;
    filename_stream.str("");
    filename_stream << "data/sd_event_"
                    << event_num  <<"_5col.dat";

    LdMatching *Matching;
    Matching = new LdMatching(13, 13, 0.1 , 0.1, 1, 0, 0, 1);
    Matching->MultiMatching(filename_stream.str(), 
        tau_min, tau_max, dtau);
    delete Matching;
  }
  return 0;
}
