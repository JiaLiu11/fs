#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <ctime>
#include "LdMatching.h"
#include "Freestreaming.h"
#include "gauss_quadrature.h"

using namespace std;

int main()
{
  double tau_min=0.0, dtau=0.1;
  double tau_max=3.;
  int nevents=20;

  //Timing the current run
  time_t start, end;
  double cpu_time_used;
  start = clock();

  //processing events
  for(int event_num=1;event_num<=nevents;event_num++)
  {
    //prepare readin filename
    ostringstream filename_stream;
    filename_stream.str("");
    filename_stream << "data/sd_event_"
                    << event_num  <<"_block.dat";

    LdMatching *Matching;
    //LdMatching(double xmax, double ymax, double dx0,double dy0,
        // int ny, double rapmin, double rapmax,
        //     int iEOS, bool outputdata)
    Matching = new LdMatching(13, 13, 0.1 , 0.1, 1, 0, 0, 1, false);
    Matching->MultiMatching(filename_stream.str(), 
        tau_min, tau_max, dtau);
    delete Matching;
  }

  end = clock();
  cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

  cout << "Time elapsed (in seconds): " << cpu_time_used << endl;
  return 0;
}
