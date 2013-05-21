#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <ctime>
#include <stdlib.h>
#include "LdMatching.h"
#include "Freestreaming.h"
#include "gauss_quadrature.h"

using namespace std;

int main(int argc, char *argv[])
{
  if(argc != 2)  //check if there is only one argument
  {
    cout << "Usage: " << argv[0] << " <total events number>" << endl;
    cout << "Try again!" << endl;
    exit(-1); 
  }

  int nevents=atoi(argv[1]);  //read in events # from command-line
  double tau_min=0.0, dtau=0.2;
  double tau_max=1.2;

  //Timing the current run
  time_t start, end;
  double cpu_time_used;
  start = clock();

  //processing events
  for(int event_num=nevents;event_num<=nevents;event_num++)
  {
    //prepare readin filename for event-by-event eccentricity fluctuation
    ostringstream filename_stream;
    filename_stream.str("");
    filename_stream << "data/sd_event_"
                    << event_num  <<"_block.dat";

    //prepare data directory for final profiles of different events and 
    //different matching time
    ostringstream result_dir_stream;
    result_dir_stream.str("");
    result_dir_stream << "./data/result/event_" << event_num;
    string result_directory = result_dir_stream.str();
    system(("rm -rf " + result_directory).c_str());
    system(("mkdir " + result_directory).c_str());

    LdMatching *Matching;
    //LdMatching(double xmax, double ymax, double dx0,double dy0,
        // int ny, double rapmin, double rapmax,
        //     int iEOS, bool outputdata, string result_location)
    Matching = new LdMatching(13, 13, 0.1 , 0.1, 1, 0, 0, 2, true, result_directory);
    Matching->MultiMatching(filename_stream.str(), 
        tau_min, tau_max, dtau);
    delete Matching;
  }

  end = clock();
  cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

  cout << "Time elapsed (in seconds): " << cpu_time_used << endl;
  return 0;
}
