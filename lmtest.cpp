#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include "LdMatching.h"
#include "Freestreaming.h"
#include "gauss_quadrature.h"
#include <ctime>

using namespace std;

int main()
{
   LdMatching *Matching;
// LdMatching(double xmax, double ymax, double ptmin, double ptmax, double dx0,double dy0,
// 		    double dpt, int nrap, double rmin, double rmax, int iEOS, string filename);

//----------------resume after debug    
    double  dtau=0.2;
    double tau0, tau1;
    int norder=1;
    
    ofstream epx_of;
    epx_of.open("data/Epx_matched_Avg_time_evolve.dat",std::ios_base::app);
    epx_of << "%  Order      Matching Time   Epx_order" <<endl;
    
    //Timing this run
    time_t start, end;
    double cpu_time_used;
    start = clock();

    //processing events
    for(int order=1;order<=norder;order++)
    {
        int loops=1;
        //prepare readin filename
        ostringstream filename_stream;
        filename_stream.str("");
        filename_stream << "data/sd_event_"
                        << order  <<"_5col.dat";
        //processing free-streaming and matching                
        for(int i=1;i<2;i++)
        {
            tau0=0.0;   //start the matching from tau0
            tau1=tau0+i*dtau;

            Matching= new LdMatching(13.0, 13.0, 0.1, 12, 0.1, 0.1,
                    0.1, 1, 0.0, 0.0, tau0, tau1, 2, filename_stream.str());
            Matching->CalTmunu(0);
            Matching->Matching_eig(1);
        
            Matching->CalPresTable();  //only can be done if ed table is complete
           Matching->GenerateSdTable(); //only can be done if ed table is complete
            Matching->CalBulkVis();
            Matching->CalShearVis();
    //         Matching->CalVis2Ideal_Ratio();
            // Matching->Test_piSq_part();
            // Matching->CalPiSquare();

            ostringstream filename_stream_ux;
            filename_stream_ux.str("");
            filename_stream_ux << "data/ux_profile_kln_tauf_" << tau1 << ".dat";
            ostringstream filename_stream_uy;
            filename_stream_uy.str("");
            filename_stream_uy << "data/uy_profile_kln_tauf_" << tau1 << ".dat";
            Matching->OutputTable_ux(filename_stream_ux.str().c_str());
            Matching->OutputTable_uy(filename_stream_uy.str().c_str());
            Matching->OutputTable_Sd("data/sd_profile.dat");

            Matching->OutputTable_BulkPi("data/bulk_pi.dat");
            epx_of << setw(8)  << setprecision(5) << order
                   << setw(12) << setprecision(5) << tau1
                   << setw(20) << setprecision(10)<< Matching->getEpx(2)<<endl;
            cout<<"Loop "<<loops<<" complete!"<<endl<<endl; 
            loops++;
            delete Matching;
        }
    }
    epx_of.close();

    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

    cout << "Time elapsed (in seconds): " << cpu_time_used << endl;

    return 0;

}
