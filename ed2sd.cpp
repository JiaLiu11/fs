#include "EOS.h"
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

int main(void)
{
  EOS eos;
  int Maxx = 261;
  int Maxy = 261;
  ifstream readin;
  readin.open("InitialSd.dat");

  ofstream of;
  of.open("test_blockTable_ed.dat",std::ios_base::out);

  eos.loadEOSFromFile("s95p-PCE/EOS_converted.dat", "s95p-PCE/coeff.dat");  //load S95 EOS table
  double sd_temp = 0;
  double ed_temp = 0;
  cout<<"Generating Entropy density table from energy density------------"<<endl;

  while(!readin.eof())
  {
  for(int i=0;i<Maxx;i++)
    {
      for(int j=0;j<Maxy;j++)
      {
        readin>>sd_temp;
        ed_temp = eos.edFromSd(sd_temp);
        //ed_temp=eos.sd(ed_temp); // Return the entropy density from the energy density ed0.
        of<<setw(22)<<setprecision(10)<<ed_temp;
      }
      of<<endl;
    }
  readin>>sd_temp;   //reach the end of the file
  cout<<"loop"<<endl;
  }


  of.close();
  readin.close();
  cout<<"Entropy density table generated!"<<endl;  
}
