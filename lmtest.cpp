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
    LdMatching *Matching;
// LdMatching(double xmax, double ymax, double ptmin, double ptmax, double dx0,double dy0,
// 		    double dpt, int nrap, double rmin, double rmax, string filename);

//----------------resume after debug    
    double tau_min=0.0, dtau=0.59;
    double tau0, tau1;
    int nevents=20;
    
    ofstream epx_of;
    epx_of.open("data/Epx_matched_diff_events.dat",std::ios_base::app);
    epx_of << "%  Event #      Matching Time     Epx_2             Epx_3" <<endl;

    //processing events
    for(int event_num=1;event_num<=nevents;event_num++)
    {
        int loops=1;
        //prepare readin filename
        ostringstream filename_stream;
        filename_stream.str("");
        filename_stream << "data/sd_event_"
                        << event_num  <<"_5col.dat";
        //processing free-streaming and matching                
        for(int i=0;i<2;i++)
        {
            tau0=0.01;   //start the matching from tau0
            tau1=tau0+i*dtau;

            Matching= new LdMatching(13.0, 13.0, 0.1, 12, 0.1, 0.1,
                    0.1, 1, 0.0, 0.0, tau_min, tau1, 1, filename_stream.str());
            Matching->CalTmunu(0);
            Matching->Matching_eig(1);
        
            Matching->CalPresTable();  //only can be done if ed table is complete
    //        Matching->GenerateSdTable(); //only can be done if ed table is complete
            Matching->CalBulkVis();
            Matching->CalShearVis();
    //         Matching->CalVis2Ideal_Ratio();
            // Matching->Test_piSq_part();
            // Matching->CalPiSquare();

            // ostringstream filename_stream_ux;
            // filename_stream_ux.str("");
            // filename_stream_ux << "data/ux_profile_kln_tauf_" << tau1 << ".dat";
            // ostringstream filename_stream_uy;
            // filename_stream_uy.str("");
            // filename_stream_uy << "data/uy_profile_kln_tauf_" << tau1 << ".dat";
            // Matching->OutputTable_ux(filename_stream_ux.str().c_str());
            // Matching->OutputTable_uy(filename_stream_uy.str().c_str());

    // //         Matching->OutputTable_Bulkpi("data/bulk_pi.dat");
    // //         Matching->OutputTmnTable("data/T00_kln.dat");
            epx_of << setw(8)  << setprecision(5) << event_num
                   << setw(12) << setprecision(5) << tau1
                   << setw(20) << setprecision(10)<< Matching->getEpx(2,0)
                   << setw(20) << setprecision(10)<< Matching->getEpx(3,0)<<endl;
            cout<<"Loop "<<loops<<" complete!"<<endl<<endl; 
            loops++;
            delete Matching;
        }
    }

    
    epx_of.close();

    return 0;

}
