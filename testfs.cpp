#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include "Freestreaming.h"
#include "gauss_quadrature.h"

using namespace std;

int main()
{
  FreeStrm *Streaming;
  ofstream of;
  of.open("data/Epxn_mckln_new.dat", std::ios_base::app);
  of<<"#Proper time                   Epx_2 ..."<<endl;
//FreeStrm(double xmax, double ymax, double ptmin, double ptmax, double dx0,double dy0,
//		    double dpt, int nrap, double rmin, double rmax, double taui, double tauf);

//   ostringstream infilename;
// for(int i_n=3;i_n<4;i_n++)  //deal with 5 events
// { 
//   of<<"#Event No."<<i_n<<endl;
// 
//   infilename.str("");
//   infilename<<"sd_event_"<<i_n<<"_5col.dat";
//   infilename.str("GaussianProfile.dat");
  double taui=0.0, tauf = 0.0;

  for(int i=1;i<18;i++)
  {
    tauf = 0.1*i + 0.2;
	  Streaming= new FreeStrm(13.0, 13.0, 0.1, 12, 0.1, 0.1,
		    0.1, 1, 0.0, 0.0, taui, tauf);
    Streaming->ReadTable("sd_event_3_5col.dat");

    Streaming->generateEdTable();

    //store epx table for various order
    of << setw(10) << setprecision(5) << tauf;
    cout << "Calculating Spatial Eccentricities------" << endl;
    for(int order = 2; order <= 3; order ++)
    {
       of << setw(16) << setprecision(8) << Streaming->getEpx(order);
    }
    of << endl;

    cout << "Loop "<<i << " complete!" <<endl <<endl;
    delete Streaming;
  }    
// }
    of.close();
    return 0;

}


