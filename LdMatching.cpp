#include <cmath>
#include <iomanip>
#include <vector>
#include "LdMatching.h"
#include "Freestreaming.h"
#include "gauss_quadrature.h"
#include "mistools.h"
#include "arsenal.h"

using namespace std;


LdMatching::LdMatching(double xmax, double ymax, double dx0,double dy0,
        int ny, double rapmin, double rapmax,
            int iEOS)
{
  Xmax=xmax;
  Ymax=ymax;
  Xmin=-xmax;
  Ymin=-ymax;

  dx=dx0;
  dy=dy0;
  Maxx=(int)((Xmax-Xmin)/dx+0.1)+1;
  Maxy=(int)((Ymax-Ymin)/dy+0.1)+1;
  nRap=ny;
  rapMin=rapmin;
  rapMax=rapmax;

  //center of the profile, max energy initialization
  Xcm = 0.;
  Ycm = 0.;
  edMax = -1.;
  dNd2rdyTable = 0;
  echo();  //output basic information of the current run

  EOS_type = iEOS;
  if(EOS_type==2)
  eos.loadEOSFromFile((char*)"s95p-PCE/EOS_converted.dat", (char*)"s95p-PCE/coeff.dat");  //load S95 EOS table
    
// // Read from streamed profile
//     Streaming->CreateDataTable("GaussianProfile.dat"); 
}


LdMatching::~LdMatching()
{
  //Clean data table
  if(dNd2rdyTable!=0)
  {
    for(int iy=0;iy<nRap;iy++)  {
        {
          for(int i=0;i<Maxx;i++) 
            delete [] dNd2rdyTable[iy][i];
        }
      delete [] dNd2rdyTable[iy];
      }
    delete [] dNd2rdyTable;      
  }  
}

void LdMatching::echo()
{
  cout<<"************************************************"<<endl<<endl;
  cout<<"Free streaming and Landau matching program"<<endl<<endl;
  cout<<"************************************************"<<endl<<endl;
  cout<<"Parameters for free-streaming and Landau matching:"<<endl
      <<"x range: "<<Xmin<<", "<<Xmax<<"     "<<endl
      <<"y range: "<<Ymin<<", "<<Ymax<<"     "<<endl
      <<"Rapidity range: "<<rapMin<<", "<<rapMax<<"     "<<endl
      <<"slicing: dx="<<dx<<", dy="<<dy<<endl;
  if(EOS_type == 1)
    cout << "EOS: Ideal gas"<<endl<<endl;
  else if(EOS_type == 2)
  {
    cout<<"EOS: S95p-PCE"<<endl<<endl;
  }
}




void LdMatching::MultiMatching(string filename, double taui, double tauf, double dtau)
{
/* Read in one matter distribution and do multiple times free-streaming and 
   Landau matching
*/
  //use KLN output file
  ReadTable(filename);
  //
  Taui = taui;
  Tauf = tauf;
  Dtau = dtau;
  int maxT = (int)((Tauf-Taui)/Dtau+0.1)+1;

  for(int t_i=1; t_i < maxT; t_i++)
  {
    double tau0 = Taui;
    double tau1 = tau0 + t_i * Dtau;
    double delta_tau = t_i*Dtau;

    cout << "Start to Free-stream the profile from tau0=" << tau0
         << " to tau=" << tau1
         << ", then do the Landau Matching." << endl;

    Streaming=new FreeStrm(Xmax, Ymax, dx, dy,
          nRap, rapMin, rapMax, tau0, tau1);
    Streaming->CopyTable(dNd2rdyTable);  //assign data table to FreeStrm class
                                         //to speed up interpolation
    //data table for LdMatching result et.al.
    DataTable=new CellData(Xmax, Ymax, dx, dy, nRap, rapMin, rapMax);
    
    CalTmunu(0, delta_tau);
    Matching_eig(1);
        
    CalPresTable();  //only can be done if ed table is ready
    //  Matching->GenerateSdTable(); //only can be done if ed table is ready
    CalBulkVis();
    CalShearVis();

    // ostringstream filename_stream_ux;
    // filename_stream_ux.str("");
    // filename_stream_ux << "data/ux_profile_kln_tauf_" << tau1 << ".dat";
    // ostringstream filename_stream_uy;
    // filename_stream_uy.str("");
    // filename_stream_uy << "data/uy_profile_kln_tauf_" << tau1 << ".dat";
    // Matching->OutputTable_ux(filename_stream_ux.str().c_str());
    // Matching->OutputTable_uy(filename_stream_uy.str().c_str());

    cout << setw(8)  << setprecision(5) << tau0
       << setw(12) << setprecision(5) << tau1
       << setw(20) << setprecision(10)<< getEpx(2,0)
       << setw(20) << setprecision(10)<< getEpx(3,0)<<endl;

    cout << "Matching is done at tau=" << tau1 
         << "!" << endl;
    //clean up before leaving
    delete DataTable;
    delete Streaming;
  }
}


void LdMatching::ReadTable(string filename)
{
  //initialize dN/dyd^2r table
  dNd2rdyTable  = new double** [nRap];
  for(int iy=0;iy<nRap;iy++) {
  dNd2rdyTable[iy] =  new double* [Maxx];
  for(int i=0;i<Maxx;i++) {
      dNd2rdyTable[iy][i] = new double [Maxy];
      for(int j=0;j<Maxy;j++) {
          dNd2rdyTable[iy][i][j]=0.0;
      }
    }
  }

  double a,b;    //temp variables help to reach the end of the data file
  ifstream DataFile(filename.c_str());
  if (!DataFile)
  {
     cout << "ReadTable::readFromFile error: file " << filename << " does not exist." << endl;
     exit(-1);
  }

  cout << "Start to read in data table " << filename << endl;
  while (!DataFile.eof())
  {
    string comment_line;  //store comment line
    getline(DataFile, comment_line);  //read in to skip the comment line
    cout << comment_line << endl;
    
    for(int iy=0;iy<nRap;iy++)
      for(int i=0;i<Maxx;i++)
          for(int j=0;j<Maxy;j++)
          {
            DataFile>>dNd2rdyTable[iy][i][j];
          }
    DataFile >> a >> b;  //safty procedure, reach the end of the file
    cout << "Read in complete!" << endl;  //this line should only be output once
  }
  DataFile.close();
  cout << "Data table has been read into memory!" << endl;
}



void LdMatching::CalTmunu(const int iRap, double delta_tau)
{
  //safty check
  if(dNd2rdyTable==0)
  {
    cout << "Calculating Tmn needs dN/dyd^2r table!" << endl;
    exit(-1);
  }
  if(delta_tau==0)
  {
    cout << "Matching at inital time! Scale the energy density afterwards!" << endl;
    delta_tau = 1;
  }
  //block for gaussian integration
  int kind=1;
  const int order=100;   //debugging
  double alpha=0.0, beta=0.0;
  double xphip[order],wphi[order];

  double phipmin=0.0, phipmax=2.0*M_PI;
  //temp variables
  double T00i=0,T01i=0,T02i=0,T11i=0,T12i=0,T22i=0;

  cout<<"Start to calculate T_mn matrix----------------------"<<endl;

  gauss_quadrature(order, kind, alpha, beta, phipmin, phipmax, xphip, wphi); 

  for(int i=0;i<Maxx;i++)  
    for(int j=0;j<Maxy;j++)  //loop over the transverse plane
    {
      for(int iphi=0;iphi<order;iphi++)  //loop over Gaussian points for phi
      {
        T00i=0,T01i=0,T02i=0,T11i=0,T12i=0,T22i=0;
        double rotate_angle = xphip[iphi];
        double rotatedDens=Streaming->GetDensity(iRap, i, j, rotate_angle);

        T00i=rotatedDens*wphi[iphi];
        T01i=rotatedDens*wphi[iphi]*cos(rotate_angle);
        T02i=rotatedDens*wphi[iphi]*sin(rotate_angle);
        T11i=rotatedDens*wphi[iphi]*cos(rotate_angle)*cos(rotate_angle);
        T12i=rotatedDens*wphi[iphi]*cos(rotate_angle)*sin(rotate_angle);
        T22i=rotatedDens*wphi[iphi]*sin(rotate_angle)*sin(rotate_angle);

        DataTable->SetTmn(iRap ,i, j, 0, 0, T00i+DataTable->GetTmn(iRap,i,j,0,0));
        DataTable->SetTmn(iRap ,i, j, 0, 1, T01i+DataTable->GetTmn(iRap,i,j,0,1));
        DataTable->SetTmn(iRap ,i, j, 0, 2, T02i+DataTable->GetTmn(iRap,i,j,0,2));
        DataTable->SetTmn(iRap ,i, j, 1, 1, T11i+DataTable->GetTmn(iRap,i,j,1,1));
        DataTable->SetTmn(iRap ,i, j, 1, 2, T12i+DataTable->GetTmn(iRap,i,j,1,2));
        DataTable->SetTmn(iRap ,i, j, 2, 2, T22i+DataTable->GetTmn(iRap,i,j,2,2));
      }
    } //<->for int i=0;i<Maxx;i++

//fill in the symmetric part of T\mu\nu, and scale the Tmn table by Tau
  for(int iy=0;iy<nRap;iy++)
    for(int i=0;i<Maxx;i++)
      for(int j=0;j<Maxy;j++)
      {
        DataTable->SetTmn(iy ,i, j, 1, 0, DataTable->GetTmn(iy,i,j,0,1)/delta_tau);
        DataTable->SetTmn(iy ,i, j, 2, 0, DataTable->GetTmn(iy,i,j,0,2)/delta_tau);
        DataTable->SetTmn(iy ,i, j, 2, 1, DataTable->GetTmn(iy,i,j,1,2)/delta_tau);

        DataTable->SetTmn(iy ,i, j, 0, 0, DataTable->GetTmn(iy,i,j,0,0)/delta_tau);
        DataTable->SetTmn(iy ,i, j, 0, 1, DataTable->GetTmn(iy,i,j,0,1)/delta_tau);
        DataTable->SetTmn(iy ,i, j, 0, 2, DataTable->GetTmn(iy,i,j,0,2)/delta_tau);
        DataTable->SetTmn(iy ,i, j, 1, 1, DataTable->GetTmn(iy,i,j,1,1)/delta_tau);
        DataTable->SetTmn(iy ,i, j, 1, 2, DataTable->GetTmn(iy,i,j,1,2)/delta_tau);
        DataTable->SetTmn(iy ,i, j, 2, 2, DataTable->GetTmn(iy,i,j,2,2)/delta_tau);
      } 
  OutputTmnTable("data/T00.dat");
  cout<<"Tmn table complete!"<<endl;
}


void LdMatching::Matching_eig(const int nRap)
{
  gsl_complex v0;
  double factor;
  double Tmn_data[16];
  double tolerance = 1e-18;  //regard quantities below this value are zero
  int count=0;

  cout<<"Start Matching-----------------------"<<endl;

  for(int iy=0;iy<nRap;iy++)
    for(int i=0;i<Maxx;i++)
      for(int j=0;j<Maxy;j++) 
      {
          int itemp=0;
          for(int ir=0;ir<4;ir++)   //copy Tmn values out
            for(int ic=0;ic<4;ic++)
            {
              Tmn_data[itemp]=DataTable->GetTmn(iy,i,j,ir,ic);
              itemp++;
            }

        gsl_matrix_view Tmn_temp 
          = gsl_matrix_view_array (Tmn_data, 4, 4);

        double data[]={1, -1, -1, -1,
                       1, -1, -1, -1,
                       1, -1, -1, -1,
                       1, -1, -1, -1};
        gsl_matrix_view gmn 
         = gsl_matrix_view_array (data, 4, 4);
        gsl_matrix_mul_elements(&Tmn_temp.matrix, &gmn.matrix);  //Equivalent to Tmn*gmn
//debug
        // cout<<"At point x= "<<Xmin+i*dx<<" y= "<<Ymin+j*dy<<" needs further check!"<<endl
        //          <<"Tmn at this point:"<<endl;
        // for(int testi=0; testi<4; testi++)
        // {
        //   for(int testj=0; testj<4; testj++)
        //     cout<<setw(21)<<setprecision(12)<<gsl_matrix_get(&Tmn_temp.matrix, testi, testj);
        //     cout<<endl;
        // }
//debug        
  // double data[]=   {64.6206632, -446.206632996633, 0, 0,
        //                446.206632996633, -4442.06632996633 ,0, 0,
        //                0,                0, -9.66666666666667, 0,
        //                0,                0, 0, 4378.11233333333
        //               };  //testing T_mn*g_mn

//convert array date[] to a matrix
        gsl_vector_complex *eval = gsl_vector_complex_alloc (4);
        gsl_matrix_complex *evec = gsl_matrix_complex_alloc (4, 4);
        gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (4);  

//Solve for eigen-values and eigen-vecotrs
        gsl_eigen_nonsymmv (&Tmn_temp.matrix, eval, evec, w);
     
        gsl_eigen_nonsymmv_free (w);

//Check if there is any non-real eigen-value is found
        // bool flag=false;
        // for (int i1 = 0; i1 < 4; i1++)
        // {
        //   gsl_complex eval_ii 
        //      = gsl_vector_complex_get (eval, i1);
        //   gsl_vector_complex_view evec_ii
        //      = gsl_matrix_complex_column (evec, i1);
        //   if(GSL_IMAG(eval_ii)!=0)
        //   {
        //   flag=true;
        //   printf("\nNon-zero eigenvalue found!");
        //   printf ("eigenvalue = %g + %gi\n",
        //           GSL_REAL(eval_ii), GSL_IMAG(eval_ii));
        //   printf ("eigenvector = \n");
        //   for (int j1 = 0; j1 < 4; ++j1)
        //     {
        //       gsl_complex z = 
        //         gsl_vector_complex_get(&evec_ii.vector, j1);
        //       printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));
        //     }
        //    }
        // }
        // gsl_vector_complex_free(eval);
        // gsl_matrix_complex_free(evec);
        // if (flag==true)
        // {
        //    for(int ir=0;ir<4;ir++)
        //      {
        //       for(int int ic=0;ic<4;ic++)
        //       cout<<setw(13)<<setprecision(10)DataTable->GetTmn(iy,i,j,ir,ic);
        //       cout<<endl;
        //      } 
        //  }
        

        for (int k = 0; k < 4; k++)
           {
             gsl_complex eval_k 
                = gsl_vector_complex_get (eval, k);

             if(GSL_REAL(eval_k)>0&&GSL_IMAG(eval_k)==0)   //select eigen-value
             {
              gsl_vector_complex_view  evec_k = 
                        gsl_matrix_complex_column (evec, k);
              v0= gsl_vector_complex_get(&evec_k.vector, 0);
              if(GSL_IMAG(v0)==0&&(2.0*GSL_REAL(v0)*GSL_REAL(v0)-1)>0)  //select eigen-vector
              {
                factor=sqrt(1.0/(2.0*GSL_REAL(v0)*GSL_REAL(v0)-1));

                if(GSL_REAL(v0)<0)
                  factor=-factor;

                if(GSL_REAL(eval_k)>tolerance)
                {
                  DataTable->SetEd(iy, i, j, GSL_REAL(eval_k));

                  for(int ivec=0;ivec<4;ivec++)
                    DataTable->SetUm(iy, i, j, ivec, factor*GSL_REAL(gsl_vector_complex_get(&evec_k.vector, ivec)));
                }

                else if(GSL_REAL(eval_k)<=tolerance)
                {
                  DataTable->SetEd(iy, i, j, 0);

                  for(int ivec=0;ivec<4;ivec++)
                    DataTable->SetUm(iy, i, j, ivec, 0);
                }

                if(count>1)
                {
                  cout<<"At point x= "<<Xmin+i*dx<<" y= "<<Ymin+j*dy<<" needs further check!"<<endl;
                  cout<<"eigen value selected: "<<scientific<<setprecision(13)<<DataTable->GetEd(iy,i,j)<<endl;
                  exit(0);
                
                }
                // count++;
              }
             }

             else continue;
           }//for k=0;k<4;k++
      }//for j=0:Maxy
  cout<<"Matching Complete!"<<endl;

  //find the weighted center of the energy density profile
  findCM_ed();
  //regulate u^\mu in the dilute region
  regulateUm();

  //prepare file name
  ostringstream filename_stream;
  filename_stream.str("");
  filename_stream << "data/ed_profile_kln_tauf_" << Tauf << ".dat";
  OutputTable_ed(filename_stream.str().c_str(), 0);  
}


void LdMatching::regulateUm(const int iRap)
/* set u1=u2=1 in dilute region, where energy density is less than max energy density by 10^10.
   Modify u0 accordingly.
*/
{
  double ed_max = edMax;
  if(ed_max<0)
  {
    cout << "Max of the energy density has not been assigned!" << endl;
    findCM_ed();
    ed_max = edMax;
    cout << "Max of the energy density is found to be: " << ed_max 
         << " Now proceed." << endl;
  }

  for(int i = 0; i < Maxx; i++)
  {
    for(int j = 0; j< Maxy; j++)
    {
      double ed_temp = DataTable->GetEd(iRap,i,j);
      if(ed_temp*1e10 < ed_max)  //rule for regulation
      {
        for(int k=0;k<3;k++)
          DataTable->SetUm(iRap, i, j, k, 1.);
      }
    }
  } //<-> for i=0:Maxx
}



void LdMatching::findCM_ed( int iRap)
{
/*Find the principal axis of the elliptic profile and align it with the x-axis
// {x} = \int{dx * dy * e(x,y) *gamma(x,u) * x} / \int{dx * dy * e(x, y) * gamma(x,y)}
// {y} = \int{dx * dy * e(x,y) *gamma(x,u) * y} / \int{dx * dy * e(x, y) * gamma(x,y)}
This is excuted after the matching
*/
  double x_ave = 0.;
  double y_ave = 0.;
  double weight_total = 0.;
  double weight_max=0.;

  for(int i = 0; i < Maxx; i++)
  {
    for(int j = 0; j< Maxy; j++)
    {
      double ed = DataTable->GetEd(iRap,i,j);
      double u0_temp = DataTable->GetUm(iRap, i, j, 0);
      double gamma = u0_temp;

      double weight = ed * gamma;   //energy density in lab frame
      if(weight > weight_max)
        weight_max = weight;

      double x_current = Xmin + dx * i;
      double y_current = Ymin + dy * j;

      x_ave += weight * x_current * dx * dy;
      y_ave += weight * y_current * dx * dy;
      weight_total += weight * dx * dy;
    }
  } //<-> for i=0:Maxx

  x_ave = x_ave/(weight_total + 1e-18);
  y_ave = y_ave/(weight_total + 1e-18);

  Xcm = x_ave;
  Ycm = y_ave;
  edMax = weight_max;

  cout << "Weighted center of profile is: " << endl
       << "("
       << Xcm << ", " << Ycm 
       << ")" <<endl;  //debug
}



double LdMatching::GetPressure(double edens)
{
  return edens/3.0;   //ideal gas
}

void LdMatching::CalPresTable(const int nRap)
{
  double pres_temp=0;

  cout<<"Start generating Pressure Table---------------------------"<<endl;
  if(EOS_type ==1)
  {
    cout<<"Generating pressure table use ideal EOS"<<endl;
    for(int iy=0;iy<nRap;iy++)
      for(int i=0;i<Maxx;i++)
        for(int j=0;j<Maxy;j++)
        {
          pres_temp=GetPressure(DataTable->GetEd(iy, i, j));
          DataTable->SetPres(iy, i, j, pres_temp);
        }
  }
  else if(EOS_type == 2)
  {
    cout<<"Generating pressure table use s95p-PCE"<<endl;
    for(int iy=0;iy<nRap;iy++)
      for(int i=0;i<Maxx;i++)
        for(int j=0;j<Maxy;j++)
        {
          pres_temp=eos.p(DataTable->GetEd(iy, i, j)); // Return the pressure from the energy density ed0.
          DataTable->SetPres(iy, i, j, pres_temp);
        }
  }

  cout<<"Pressure table complete!----------------------------"<<endl<<endl;
}

void LdMatching::GenerateSdTable(const int nRap)
{
  double sd_temp = 0;
  cout<<"Generating Entropy density table from energy density------------"<<endl;
  for(int iy=0;iy<nRap;iy++)
    for(int i=0;i<Maxx;i++)
      for(int j=0;j<Maxy;j++)
      {
        sd_temp=eos.sd(DataTable->GetEd(iy, i, j)); // Return the entropy density from the energy density ed0.
        DataTable->SetSd(iy, i, j, sd_temp);
      }
  cout<<"Entropy density table generated!"<<endl;  
}


void LdMatching::CalBulkVis(const int nRap)
{
  double T00i, T01i, T02i, T03i, T11i, T12i,  T13i,
         T22i, T23i, T33i;
  double u0i, u1i, u2i, u3i;
  double tr, result;
  double tolerance=1e-10;

  cout<<"Calculating Bulk Viscosity--------------------------"<<endl;

  for(int iy=0;iy<nRap;iy++)
    for(int i=0;i<Maxx;i++)
      for(int j=0;j<Maxy;j++)
      {
        tr=0;

        T00i = DataTable->GetTmn(iy, i, j, 0, 0);
        T01i = DataTable->GetTmn(iy, i, j, 0, 1);
        T02i = DataTable->GetTmn(iy, i, j, 0, 2);
        T03i = DataTable->GetTmn(iy, i, j, 0, 3);
        T11i = DataTable->GetTmn(iy, i, j, 1, 1);
        T12i = DataTable->GetTmn(iy, i, j, 1, 2);
        T13i = DataTable->GetTmn(iy, i, j, 1, 3);
        T22i = DataTable->GetTmn(iy, i, j, 2, 2);
        T23i = DataTable->GetTmn(iy, i, j, 2, 3);
        T33i = DataTable->GetTmn(iy, i, j, 3, 3);

        u0i  = DataTable->GetUm(iy, i, j, 0);
        u1i  = DataTable->GetUm(iy, i, j, 1);
        u2i  = DataTable->GetUm(iy, i, j, 2);
        u3i  = DataTable->GetUm(iy, i, j, 3);

        tr = ((1.-u0i*u0i)*T00i - (1.+u1i*u1i)*T11i -(1.+u2i*u2i)*T22i - (1.+u3i*u3i)*T33i
              + 2.*u0i*u1i*T01i + 2.*u0i*u2i*T02i   +2.*u0i*u3i*T03i
              - 2.*u1i*u2i*T12i  - 2.*u1i*u3i*T13i   -2.*u2i*u3i*T23i);            
        result = -1./3.*tr- DataTable->GetPres(iy, i, j);
        
       //  if(i==116 && j==93)  //debug
       //  cout << "bulk pi trace" << tr << endl
       //       << "bulk pi pressure" << DataTable->GetPres(iy, i, j) << endl
						 // << "bulk pi result" << result << endl;
//         if(i==57&&j==129) 
// 				if(result >0.1&&DataTable->GetPres(iy, i, j)>1)  //debugging
//         {
//           cout<<"u_m at this point: "<<endl
//               <<scientific<<setprecision(10)<<u0i<<" "<<u1i<<" "<<u2i<<" "<<u3i<<endl;
//           cout<<"energy density"<<DataTable->GetEd(iy, i, j)<<endl;
//           cout<<"Tmn at this point: "<<endl
//               <<scientific<<setprecision(10)<< T00i<<" "<<T01i<<" "<<T02i<<" "<<T03i<<endl
//                                             << T01i<<" "<<T11i<<" "<<T12i<<" "<<T13i<<endl
//                                             << T02i<<" "<<T12i<<" "<<T22i<<" "<<T23i<<endl
//                                             << T03i<<" "<<T13i<<" "<<T23i<<" "<<T33i<<endl;
//           cout<< tr<<endl
//             << result<<endl; 
//             exit(0);
//         }//debug

        if(fabs(result)<tolerance) 
          result = 0;

        DataTable->SetBulk_Pi(iy, i, j, result);

      }
  OutputTable_BulkPi("data/bulk_pi.dat");

  cout<<"Bulk viscosity table complete!"<<endl<<endl; 

}



void LdMatching::CalShearVis(const int nRap)
{
  double Tmn_data[16];  //temp Tmn table
  double tolerance=1e-15;
  double u0i, u1i, u2i, u3i;  //temp u_mu table

  cout<<"Calculating shear viscosity----------------"<<endl;

  ofstream logfile;
  logfile.open("logfile_pimu_trace.dat", std::ios_base::out);  //keep a 
  logfile<<"%Trace of shear viscosity largger than 1e-8"<<endl
         <<"% i     j      Trace"<<endl;

  for(int iy=0;iy<nRap;iy++)
    for(int i=0;i<Maxx;i++)
      for(int j=0;j<Maxy;j++)
      { 
        int itemp=0; 
        double factor_um = 0;
        double factor_gmn = 0;

        for(int ir=0;ir<4;ir++)   //copy Tmn values out
          for(int ic=0;ic<4;ic++)
            {
              Tmn_data[itemp]=DataTable->GetTmn(iy,i,j,ir,ic);
              itemp++;
            }

        gsl_matrix_view Tmn_temp 
          = gsl_matrix_view_array (Tmn_data, 4, 4);

        u0i  = DataTable->GetUm(iy, i, j, 0);
        u1i  = DataTable->GetUm(iy, i, j, 1);
        u2i  = DataTable->GetUm(iy, i, j, 2);
        u3i  = DataTable->GetUm(iy, i, j, 3);

        double b[]={u0i*u0i,  u0i*u1i,  u0i*u2i,  u0i*u3i,
                    u0i*u1i,  u1i*u1i,  u1i*u2i,  u1i*u3i,
                    u0i*u2i,  u1i*u2i,  u2i*u2i,  u2i*u3i,
                    u0i*u3i,  u1i*u3i,  u2i*u3i,  u3i*u3i};
        gsl_matrix_view Umn_temp 
          = gsl_matrix_view_array (b, 4, 4);        
        factor_um = DataTable->GetEd(iy, i, j)+ DataTable->GetPres(iy, i, j)       
                   +DataTable->GetBulk_Pi(iy, i, j) ;
        gsl_matrix_scale(&Umn_temp.matrix, factor_um);   //(e+p+Bulk_pi)*u'u

        double data[]={1,  0,  0,  0,
                       0, -1,  0,  0, 
                       0,  0, -1,  0,
                       0,  0,  0, -1};
        gsl_matrix_view g_mn 
          = gsl_matrix_view_array (data, 4, 4);
        factor_gmn = DataTable->GetPres(iy, i, j)       
                          +DataTable->GetBulk_Pi(iy, i, j) ;
        gsl_matrix_scale(&g_mn.matrix, factor_gmn);  //(p+Bulk_pi)*g_mn


//Pi_mu=T_mn-(e+p+Bulk_pi)*u'u+(e+Bulk_pi)*g_mn
        gsl_matrix_sub(&Tmn_temp.matrix, &Umn_temp.matrix);
        gsl_matrix_add(&Tmn_temp.matrix, &g_mn.matrix);  

// //trace of Pi_mn for debug
        // double tr_pi=0;
        // gsl_matrix *temp0 = gsl_matrix_alloc(4,4);
        // gsl_matrix_memcpy (temp0, &Tmn_temp.matrix);
        // Lower_matrix_single(temp0);

        // for(int k=0;k<4;k++)
        //   tr_pi+=gsl_matrix_get(temp0, k, k);
        // gsl_matrix_free(temp0);

        // if(tr_pi>1e-8) 
        //   {
        //     logfile<<i<<" "<<j<<" "<<scientific<<tr_pi<<endl;
        //     cout<<"Trace of shear viscosity matrix out of range: "
        //                    <<setw(22)<<setprecision(10)<<tr_pi<<endl;
        //     cout<<"At this point: "<<Xmin+i*dx<<" "<<Ymin+j*dy<<endl;
        //     cout<<"i= "<<i<<" j="<<j<<endl;

        //     cout<<"Tmn at this point: "<<endl;
        //     gsl_matrix *Tmn_current = gsl_matrix_alloc(4,4);
        //     DataTable->Tmn_matrix_get(Tmn_current, iy, i, j);

        //     for(int ip=0;ip<4;ip++)
        //     {
        //       for(int jp=0;jp<4;jp++)
        //         cout<<scientific<<setw(18)<<setprecision(10)<<gsl_matrix_get(Tmn_current, ip, jp);
        //       cout<<endl;
        //     }

        //     cout<<"U_mu at this point: "<<endl;
        //     for(int ip=0;ip<4;ip++)
        //       cout<<scientific<<setw(18)<<setprecision(10)<<DataTable->GetUm(iy, i, j, ip);
        //     cout<<endl;

        //     cout<<"Energy density at this point: "<<endl
        //         <<scientific<<setw(18)<<setprecision(10)<<DataTable->GetEd(iy, i, j)<<endl;

        //     cout<<"Bulk viscosity at this point: "<<endl
        //         <<scientific<<setw(18)<<setprecision(10)<<DataTable->GetBulk_Pi(iy, i, j)<<endl; 

        //     cout<<"Pi_mn at this point:"<<endl;    
        //     for(int ir=0;ir<4;ir++)
        //     {
        //       for(int jr=0;jr<4;jr++)
        //         cout<<scientific<<setw(18)<<setprecision(10)<<gsl_matrix_get(&Tmn_temp.matrix, ir, jr);
        //       cout<<endl;
        //     }
            
        //     gsl_matrix_free(Tmn_current);
        //     exit(0);
        //     }
//Store Pi_mn table
        for(int ir=0;ir<4;ir++)
          for(int ic=0;ic<4;ic++)
          {
            if(fabs(gsl_matrix_get(&Tmn_temp.matrix, ir, ic))>tolerance)
              DataTable->SetPi_mn(iy, i, j, ir, ic, gsl_matrix_get(&Tmn_temp.matrix, ir, ic));
            else 
              DataTable->SetPi_mn(iy, i, j, ir, ic, 0);
            // cout<<DataTable->GetPi_mn(iy, i, j, ir, ic)<<endl;
          }

//         if(i==116 && j==93)  //debug
//         Diagnostic(iy, i, j);
      }
  logfile.close();
  cout<<"Shear viscosity table Pi_mu nu complete!"<<endl<<endl;
  OutputTables_Pimn(0);
}



void LdMatching::CalVis2Ideal_Ratio(const int iRap)
{
  double tr;
  // double temp;
  double tolerance=1e-7;

  ofstream of;
  of.open("pi_directCal.dat",std::ios_base::out); //debug
  of<<" %Pi*Pi directed calculated for comparison"<<endl; //debug

  cout<<"Calculating Viscosity vs. Energy momentum tensor ratio---------"<<endl;
  for(int i=0;i<Maxx;i++)
  {
    for(int j=0;j<Maxy;j++)
    {
      tr = 0.0;
      double data[]={1,  0,  0,  0,
                     0, -1,  0,  0, 
                     0,  0, -1,  0,
                     0,  0,  0, -1};
      gsl_matrix_view g_mn = gsl_matrix_view_array (data, 4, 4);

      double gmn_convert[]={ 1, -1, -1, -1,
                            -1,  1,  1,  1,
                            -1,  1,  1,  1,
                            -1,  1,  1,  1};
      double c[] = { 0.00, 0.00, 0.00, 0.00,
               0.00, 0.00, 0.00, 0.00,
               0.00, 0.00, 0.00, 0.00,
               0.00, 0.00, 0.00, 0.00 }; 

      gsl_matrix_view gmn_converter
          =gsl_matrix_view_array (gmn_convert, 4, 4);

      gsl_matrix *pi_mn_up=gsl_matrix_alloc(4,4);
      DataTable->Pimn_matrix_get(pi_mn_up, iRap, i, j);  //Pi^{\mu\nu}

      gsl_matrix *pi_mn_low=gsl_matrix_alloc(4,4);
      gsl_matrix_memcpy(pi_mn_low, pi_mn_up); 
      gsl_matrix_mul_elements(pi_mn_low, &gmn_converter.matrix);  //Pi^{\mu\nu}--->Pi_{\mu\nu}

      gsl_matrix_view result = gsl_matrix_view_array(c, 4, 4);

//Matrix multiplication: Pi^{\mu\nu}*Pi_{\mu\nu}
      gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                 1.0, pi_mn_up, pi_mn_low,
                 0.0, &result.matrix);

      for(int k=0;k<4;k++)
        tr+=gsl_matrix_get(&result.matrix,k,k);
      if(fabs(tr)<tolerance)
        tr=0.0;
      else if(fabs(tr)>tolerance&&tr<0)
      {
        cout<<"Pi_mn square is smaller than zero: sqrt(pi^mn * pi_mn)="<<tr<<endl;
//        exit(0);
       tr = 0.0;//debug
      }

      of << setw(22) <<setprecision(10) << tr; //debug

      double a0=sqrt(tr);

//calculate T^mn*T_mn for ideal case
      double b0=0, u0i, u1i, u2i, u3i;
      u0i  = DataTable->GetUm(iRap, i, j, 0);
      u1i  = DataTable->GetUm(iRap, i, j, 1);
      u2i  = DataTable->GetUm(iRap, i, j, 2);
      u3i  = DataTable->GetUm(iRap, i, j, 3);

      double b[]={u0i*u0i,  u0i*u1i,  u0i*u2i,  u0i*u3i,
                  u0i*u1i,  u1i*u1i,  u1i*u2i,  u1i*u3i,
                  u0i*u2i,  u1i*u2i,  u2i*u2i,  u2i*u3i,
                  u0i*u3i,  u1i*u3i,  u2i*u3i,  u3i*u3i};
      gsl_matrix_view Umn_temp 
          = gsl_matrix_view_array (b, 4, 4);        
      double factor_um = DataTable->GetEd(iRap, i, j)+ DataTable->GetPres(iRap, i, j);
      gsl_matrix_scale(&Umn_temp.matrix, factor_um);   //(e+p)*u'u

      double factor_gmn = DataTable->GetPres(iRap, i, j);
      gsl_matrix_scale(&g_mn.matrix, factor_gmn);  //p*g_mn

      gsl_matrix_sub(&Umn_temp.matrix, &g_mn.matrix);// T^mn ideal is stored in Umn_temp.matrix
      gsl_matrix *Tmn_ideal_up=gsl_matrix_alloc(4,4);  //Store T^mn
      gsl_matrix_memcpy(Tmn_ideal_up, &Umn_temp.matrix);

      gsl_matrix_mul_elements(&Umn_temp.matrix, &gmn_converter.matrix);  //T^mn--->T_mn, then store in Umn_temp

//Matrix multiplication: T^{\mu\nu}*T_{\mu\nu}
      gsl_matrix_view result2 = gsl_matrix_view_array(c, 4, 4);
      gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                 1.0, Tmn_ideal_up, &Umn_temp.matrix,
                 0.0, &result2.matrix);
      for(int k=0;k<4;k++)
        b0+=gsl_matrix_get(&result2.matrix,k,k);


      if(fabs(b0)<tolerance)
        b0=0.0;
      else if(fabs(b0)>tolerance&&b0<0)
      {
        cout<<"T_mn square is smaller than zero: sqrt(T^mn * T_mn)="<<b0<<endl;
        exit(0);
      }
      b0=sqrt(b0);

      double ratio=a0/(b0+1e-16);
      DataTable->SetVisRatio(iRap, i, j, ratio);
                    
// //debug
      // temp=0.0;
      // for(int ir=0;ir<4;ir++)
      //   for(int ic=0;ic<4;ic++)
      //   {
      //     temp+=gsl_matrix_get(pi_mn_up, ir, ic)*gsl_matrix_get(pi_mn_low, ir, ic);
      //   }
      // cout<<"Trace sum gives: "<<tr<<endl
      //     <<"Summation gives: "<<temp<<endl;
//debug
      gsl_matrix_free(pi_mn_up);
      gsl_matrix_free(pi_mn_low);
      gsl_matrix_free(Tmn_ideal_up);
    }
    of << endl;  //debug
  }
  of.close();
  cout<<"Ratio complete!"<<endl;
  Output4colTable_visratio("data/visideal_ratio_kln.dat", 0);
}


// void LdMatching::Test_piSq_part(const int iRap)
// {
// //Calculating integration over p and phi_p to get one part of test function
// //block for gaussian integration
//   int kind=1;
//   const int order=50;   //debugging
//   double alpha=0.0, beta=0.0;
//   double xphip[order],wphi[order];
//   double xpt[order],wpt[order];
//   double elem;

//   double phipmin=0.0, phipmax=2.0*M_PI;
// //temp
//   double result=0.;

//   cout<<"Start to calculate test t_mn------------------------"<<endl;
// //free-streamed gluon density table for a certain phi_p
//   double ****dNd2rdyPhipTable;
//     dNd2rdyPhipTable  = new double*** [nRap];
//     for(int iy=0;iy<nRap;iy++) 
//     {
//       dNd2rdyPhipTable[iy] =  new double** [Maxx];
//       for(int i=0;i<Maxx;i++) 
//       {
//         dNd2rdyPhipTable[iy][i] = new double* [Maxy];
//         for(int j=0;j<Maxy;j++) 
//         {
//             dNd2rdyPhipTable[iy][i][j]=new double[order];
//               for(int iphip=0;iphip<order; iphip++)
//                     dNd2rdyPhipTable[iy][i][j][iphip]=0.0;
//         }
//       }
//     } 

//   gauss_quadrature(order, kind, alpha, beta, phipmin, phipmax,xphip, wphi); 
//   gauss_quadrature(order, kind, alpha, beta, PTmin, PTmax, xpt, wpt); 

//   for(int iphi=0;iphi<order;iphi++)
//     {
//      Streaming->ShiftDensity(0, xphip[iphi]); 


// //Initialize vectors for cubic interpolation
//       vector<double>* dens1=new vector<double>(MaxPT,0.0);
//       vector<double>* pt0=new vector<double>(MaxPT,0.0);    
//       for(int ipt0=0;ipt0<MaxPT;ipt0++)
//         (*pt0)[ipt0]=(PTmin+dpt*ipt0);

// //inner integration: integrate over pt for a certain phi_p.
//       for(int i=0;i<Maxx;i++)  {
//         for(int j=0;j<Maxy;j++) {
//           for(int ipt=0;ipt<order;ipt++)
//           {

//             for(int k=0;k<MaxPT;k++)
//               (*dens1)[k]=Streaming->GetShiftdeDensity(iRap, i, j, k); //debug
// //Use 1D cubic interpolation interpCubicDirect(vector<double>* x, vector<double>* y, double x0)
// // to interpolate 
//           elem=interpCubicDirect(pt0, dens1, xpt[ipt]);

//           dNd2rdyPhipTable[iRap][i][j][iphi]+=elem*wpt[ipt]*xpt[ipt]*xpt[ipt];
//          }
//        }

//       } //for int i=0;i<Maxx;i++
// //temp variables:

//     for(int iy=0;iy<nRap;iy++)
//       for(int i=0;i<Maxx;i++)
//         for(int j=0;j<Maxy;j++)
//         {
//             result=0;
//             double um[4]={0, 0, 0, 0};
//             for(int ic=0;ic<4;ic++)
//               um[ic]=DataTable->GetUm(iy, i, j, ic);

//             double factor = (um[0]-um[1]*cos(xphip[iphi])-um[2]*sin(xphip[iphi]))
//                            *(um[0]-um[1]*cos(xphip[iphi])-um[2]*sin(xphip[iphi]))
//                            /Tauf ;

//             result=dNd2rdyPhipTable[iRap][i][j][iphi]*wphi[iphi]
//                      *factor;
//             DataTable->SetupVal(iy ,i, j, result+DataTable->GetupVal(iy,i,j));
//         } //for int j...   
//   } //for int iphi=0;iphi<order;iphi++

// //fill in the symmetric part of T\mu\nu, and scale the Tmn table by Tau
 

//   //Clean intermediate table
//     for(int iy=0;iy<nRap;iy++) 
//     {
//       for(int i=0;i<Maxx;i++) 
//       {
//         for(int j=0;j<Maxy;j++) 
//           delete [] dNd2rdyPhipTable[iy][i][j];
//         delete [] dNd2rdyPhipTable[iy][i];
//         }
//     delete [] dNd2rdyPhipTable[iy];
//     }
//     delete [] dNd2rdyPhipTable; 
//     cout<<"Inte[U*P*f] table complete!"<<endl;
// }

// void LdMatching::CalPiSquare(const int iRap)
// {
// // Calculate Pi^mn * Pi_mn use another way, do a cross check to the previous code
//   double result=0.;
//   for(int i=0;i<Maxx;i++)
//     for(int j=0;j<Maxy;j++)
//     {
//       gsl_matrix *Tmn_upper=gsl_matrix_alloc(4,4);
//       DataTable->Tmn_matrix_get(Tmn_upper, iRap, i, j);  //Get T^mn for the current point

//       double TSq = Contract_matrix(Tmn_upper);  //TSq = T^mn * T_mn, contraction
//       double upinte = DataTable->GetupVal(iRap, i, j);
//       double pressure = DataTable->GetPres(iRap, i, j);
//       result = TSq - upinte*upinte -3*pressure*pressure;
//       DataTable->SetPiSq(iRap, i, j, result);
//     }  

//   Output_picontract_comp("pi_contraction_comp.dat");
//   cout<<"Pi*Pi completed!"<<endl;

// }

double LdMatching::getShiftedProfile(double ***data, int i, int j, double x0, double y0, double phi, bool limit, const int iRap)
{
    double x,y;     //unshifted coordinate
    double shfedx, shfedy;    //shifed coordinate
    double Phi;
    double it, jt;             //indices corresponding to shifted coordinates

    const int method=2;  //method=1 : Bilinear interpolation, no worry about negative value of ed;
											   //method=2 : Cubic interpolation, better precision, but ed may be negative. 
                         //           Set limit=true to avoid negative ed.
    Phi=phi;

    x= Xmin + i * dx;
    y= Ymin + j * dy;

    //Shift the coordinate
    if( Phi ==0 )
    {
      shfedx = x + x0;
      shfedy = y + y0;
    }
    else
    {
      shfedx= (x - x0) * cos(Phi) + (y - y0) * sin(Phi);
      shfedy=-(x - x0) * sin(Phi) + (y - y0) * cos(Phi);
    }

    // finite grid size
    if( (shfedx>Xmax||shfedx<Xmin) ||(shfedy>Ymax||shfedy<Ymin) )  //coordinates locate outside of the grid
        return 0;

    else
    {
      it=(shfedx-Xmin)/dx;        
      jt=(shfedy-Ymin)/dy;

      // boundary safty control
      if (jt<=0) jt += 1e-10;
      if (jt>=Maxy-1) jt -= 1e-10;
      if (it<=0) it += 1e-10;
      if (it>=Maxx-1) it -= 1e-10;
      // get integer parts:
      long int jti = (long int)floor(jt);
      long int iti = (long int)floor(it);

      if(method==1)  //bilinear interpolation
      {
        if (jti<0) jti=0;
        if (jti>=Maxy-2) 
          {
            jti=Maxy-2; // need 2 points
            // cout<<"out of y boundary"<<endl;
          }

        if (iti<0) iti=0;       
        if (iti>=Maxx-2) 
          { 
            iti=Maxx-2; // need 2 points
           // cout<<"out of x boundary"<<endl;
          }
        return Bilinear2dInterp(iti, jti, 1, 1, data[iRap][iti][jti], data[iRap][iti][jti+1], 
            data[iRap][iti+1][jti+1], data[iRap][iti+1][jti]);
      } 

      if(method==2)  //cubic interpolation
      {
        // Cubic interpolation
        if (jti<0) jti=0;
        if (jti>=Maxy-4) 
          {
            jti=Maxy-4; // need 4 points
            // cout<<"out of y boundary"<<endl;
          }

        if (iti<0) iti=0;       
        if (iti>=Maxx-4) 
          { 
            iti=Maxx-4; // need 4 points
           // cout<<"out of x boundary"<<endl;
          }
        //interpolation is done on the x-y grid, 

        double xfraction = it-iti;
        double yfraction = jt-jti;
         // row interpolation + extrapolation first
        if(limit==true)
        {
          double A0 = interpCubic4Points(data[iRap][iti][jti], data[iRap][iti][jti+1], 
            data[iRap][iti][jti+2], data[iRap][iti][jti+3], 1, yfraction, limit);

          double A1 = interpCubic4Points(data[iRap][iti+1][jti], data[iRap][iti+1][jti+1], 
            data[iRap][iti+1][jti+2], data[iRap][iti+1][jti+3], 1, yfraction, limit);

          double A2 = interpCubic4Points(data[iRap][iti+2][jti], data[iRap][iti+2][jti+1], 
            data[iRap][iti+2][jti+2], data[iRap][iti+2][jti+3], 1, yfraction, limit);

          double A3 = interpCubic4Points(data[iRap][iti+3][jti], data[iRap][iti+3][jti+1], 
            data[iRap][iti+3][jti+2], data[iRap][iti+3][jti+3], 1, yfraction, limit);

          return interpCubic4Points(A0,A1,A2,A3,1, xfraction, limit);
        }//<-> if method ==1
        else
        {
          double A0 = interpCubic4Points(data[iRap][iti][jti], data[iRap][iti][jti+1], 
            data[iRap][iti][jti+2], data[iRap][iti][jti+3], 1, yfraction);

          double A1 = interpCubic4Points(data[iRap][iti+1][jti], data[iRap][iti+1][jti+1], 
            data[iRap][iti+1][jti+2], data[iRap][iti+1][jti+3], 1, yfraction);

          double A2 = interpCubic4Points(data[iRap][iti+2][jti], data[iRap][iti+2][jti+1], 
            data[iRap][iti+2][jti+2], data[iRap][iti+2][jti+3], 1, yfraction);

          double A3 = interpCubic4Points(data[iRap][iti+3][jti], data[iRap][iti+3][jti+1], 
            data[iRap][iti+3][jti+2], data[iRap][iti+3][jti+3], 1, yfraction);

          return interpCubic4Points(A0,A1,A2,A3,1, xfraction);
        }
      }//<->if method ==2
      else 
      {
        cout << "Method:" << method << "  is not supported!" << endl;
        exit(0);
      }
    }
}

double LdMatching::getEpx_test(const int iRap)
// Epx = \int{dx*dy*e(x,y)*(y^2-x^2)*gamma(ux,uy)}/\int{dx*dy*e(x,y)*(y^2+x^2)*gamma(ux,uy)}
// gamma factor is inlucded, since Landau matching generates flow, in order to transform energy
// density to the lab frame, gamma factor should be included.
{
  cout<<"Calculating spatial eccentricity Epx----------------"<<endl;
  double Epx=0.0;

  double epx_nu_real = 0.0;  //numerator of epx
  double epx_nu_img = 0.0;
  double Epx_angle = 0.0;
  double epx_dn = 0.0;  //denominator of epx

  //load in ed, ux, uy table
  double ***ed_temp  = new double** [nRap];
  double ***ux_temp  = new double** [nRap];
  double ***uy_temp  = new double** [nRap];
  ux_temp  = new double** [nRap];
  for(int iy=0;iy<nRap;iy++) 
  {
    ed_temp[iy] =  new double* [Maxx];
    ux_temp[iy] =  new double* [Maxx];
    uy_temp[iy] =  new double* [Maxx];
    
    for(int i=0;i<Maxx;i++) 
    {
      ed_temp[iy][i] = new double [Maxy];
      ux_temp[iy][i] = new double [Maxy];  
      uy_temp[iy][i] = new double [Maxy];
     
      for(int j=0;j<Maxy;j++) 
      {
        ed_temp[iy][i][j]=DataTable->GetEd(iy, i, j);
        ux_temp[iy][i][j]=DataTable->GetUm(iy, i, j, 1);
        uy_temp[iy][i][j]=DataTable->GetUm(iy, i, j, 2);
      }
    }
  }//<-> for iy=0:nRap

//printout shifted profile, debug
  cout << "output shifted profile for check!" << endl;
	ofstream ed_of, ux_of, uy_of, gamma_of;
  ed_of.open("data/ed_kln_shifted.dat",std::ios_base::out);
  ux_of.open("data/ux_kln_shifted.dat",std::ios_base::out);
  uy_of.open("data/uy_kln_shifted.dat",std::ios_base::out);
  gamma_of.open("data/gamma_kln_shifted.dat",std::ios_base::out);
  for(int i=0; i<Maxx; i++)
  {
      for(int j=0; j<Maxy; j++)
     {
       double x = Xmin + dx*i;
       double y = Ymin + dy*j;

       double ed_shifted = getShiftedProfile(ed_temp, i, j, Xcm, Ycm, 0, true, 0);  //forbid ed=negative

       double ux_shifted = getShiftedProfile(ux_temp, i, j, Xcm, Ycm, 0, false, 0);
       double uy_shifted = getShiftedProfile(uy_temp, i, j, Xcm, Ycm, 0, false, 0);

       double gamma = sqrt(1. + ux_shifted * ux_shifted + uy_shifted * uy_shifted);

       //output shifted profiles , debug
			 ed_of << setw(20) << setprecision(10) << ed_shifted;
       ux_of << setw(20) << setprecision(10) << ux_shifted;
	     uy_of << setw(20) << setprecision(10) << uy_shifted;
			 gamma_of << setw(20) << setprecision(10) << gamma;

       epx_nu_real += ed_shifted * (y*y - x*x)* gamma * dx*dy;
       epx_nu_img += ed_shifted * ( y * x) * gamma * dx * dy;  //imaginary part
       epx_dn += ed_shifted * (y*y + x*x)* gamma * dx*dy;
     }       
			 ed_of << endl;   //debug
       ux_of << endl;
			 uy_of << endl;
			 gamma_of << endl;
  }//<-> for i=0:Maxx
			ed_of.close();
      ux_of.close();
      uy_of.close();
			gamma_of.close();
  
  Epx = sqrt( epx_nu_real * epx_nu_real + 4.0 * epx_nu_img * epx_nu_img )
       /(epx_dn + 1e-18);
  Epx_angle = atan2(2.0 * epx_nu_img, (epx_nu_real + 1e-18));

  cout<<"Spatial Eccentricity complete!"<<endl
      <<"epx real=" << epx_nu_real << ", epx img=" << -2*epx_nu_img << endl;
  cout << "Epx_angle =" <<Epx_angle << endl;

  //clean up before leave
  for(int iy=0;iy<nRap;iy++)  
  {
    for(int i=0;i<Maxx;i++) 
    {
      delete [] ed_temp[iy][i];
      delete [] ux_temp[iy][i];
      delete [] uy_temp[iy][i];
    }
  
    delete [] ed_temp[iy];
    delete [] ux_temp[iy];
    delete [] uy_temp[iy];
  }
  return Epx;
}


double LdMatching::getEpx(int nth_order, const int iRap)
// Epx = \int{dx*dy*e(x,y)*(y^2-x^2)*gamma(ux,uy)}/\int{dx*dy*e(x,y)*(y^2+x^2)*gamma(ux,uy)}
// gamma factor is inlucded, since Landau matching generates flow, in order to transform energy
// density to the lab frame, gamma factor should be included.
{
  cout<<"Calculating spatial eccentricity Epx----------------"<<endl;
  double Epx=0.0;

  double epx_nu_real = 0.0;  //numerator of epx
  double epx_nu_img = 0.0;
  double Epx_angle = 0.0;
  double epx_dn = 0.0;  //denominator of epx
  
  for(int i=0; i<Maxx; i++)
  {
      for(int j=0; j<Maxy; j++)
     {
       double x = Xmin + dx*i - Xcm;  //center the profile
       double y = Ymin + dy*j - Ycm;
       double phi = atan2(y,x);

       double ed_temp = DataTable->GetEd(iRap, i, j);
       double ux_temp = DataTable->GetUm(iRap, i, j, 1);
       double uy_temp = DataTable->GetUm(iRap, i, j, 2);
       double gamma = sqrt(1 + ux_temp * ux_temp + uy_temp * uy_temp);
       double weight = ed_temp * gamma;

       epx_nu_real += pow(x*x + y*y, double(nth_order)/2.) 
                  * cos(nth_order * phi) * weight * dx*dy;
       epx_nu_img += pow(x*x + y*y, double(nth_order)/2.) 
                  * sin(nth_order * phi) * weight * dx*dy;
       epx_dn += pow(y*y + x*x, double(nth_order)/2.)* weight * dx*dy;
     }       
  }//<-> for i=0:Maxx
  
  Epx = sqrt( epx_nu_real * epx_nu_real + epx_nu_img * epx_nu_img )
       /(epx_dn + 1e-18);
  Epx_angle = -atan2( epx_nu_img, (epx_nu_real + 1e-18))/nth_order;

  cout<<"Spatial Eccentricity complete!"<<endl
      <<"epx real=" << epx_nu_real <<", epx imaginary=" << epx_nu_img << endl;
  cout << "Epx_angle =" <<Epx_angle << endl;

  return Epx;
}

void LdMatching::OutputTable_ux(const char *filename, const int iRap)
{
  ofstream of;
  // of.open(filename, std::ios_base::out | std::ios_base::app);

  //   for(int i=0;i<Maxx;i++)      
  //   {
  //     for(int j=0;j<Maxy;j++)
  //     of  << Tauf-Taui<<"   "
  //         << Xmin+i*dx<<"   "
  //         << Ymin+j*dy<<"   "
  //         << setprecision(10) << setw(14) << DataTable->GetUm(iRap,i,j,1)/(DataTable->GetUm(iRap,i,j,0)+1e-16)<<endl;  
//vx = u1 / u0
  //   }
  // cout<<"vx has been printed out for time Delta_Tau= "<<Tauf-Taui<<" fm/c"<<endl;
//output u instead
    of.open(filename, std::ios_base::out);

    for(int i=0;i<Maxx;i++)      
    {
      for(int j=0;j<Maxy;j++)
      of << scientific << setprecision(10) << setw(20) << DataTable->GetUm(iRap,i,j,1);  
      of << endl;
    }
  cout<<"ux has been printed out for time Delta_Tau= "<<Tauf-Taui<<" fm/c"<<endl;

  of.close();
}



void LdMatching::OutputTable_uy(const char *filename, const int iRap)
{
  ofstream of;
//   of.open(filename, std::ios_base::out | std::ios_base::app);

// //  Output Tmn table

//     for(int i=0;i<Maxx;i++)      
//     {
//       for(int j=0;j<Maxy;j++)
//       of  << Tauf-Taui<<"   "
//           << Xmin+i*dx<<"   "
//           << Ymin+j*dy<<"   "
//           << setprecision(10) << setw(14) << DataTable->GetUm(iRap,i,j,2)/(DataTable->GetUm(iRap,i,j,0)+1e-16)<<endl;  
//vy = u1 / u0
//     }
  //output u instead
    of.open(filename, std::ios_base::out);

    for(int i=0;i<Maxx;i++)      
    {
      for(int j=0;j<Maxy;j++)
      of << scientific << setprecision(10) << setw(20) << DataTable->GetUm(iRap,i,j,2);
      of << endl; 
    }
  cout<<"uy has been printed out for time Delta_Tau= "<<Tauf-Taui<<" fm/c"<<endl;

  of.close();
}
//


// void LdMatching::OutputTable_vz(const char *filename, const int nRap)
// {
//   ofstream of;
//   of.open(filename, std::ios_base::out);

// //  Output Tmn table
//   for(int iRap=0;iRap<nRap;iRap++)
//     for(int i=0;i<Maxx;i++)
//     {
//       for(int j=0;j<Maxy;j++)
//       of  << setprecision(12) << setw(22) << uz[iRap][i][j]/(u0[iRap][i][j]+1e-16);  //need revise
//       of  << endl;
//     }


//   cout<<"vz has been printed out!"<<endl;
//   of.close();
// }

void LdMatching::Output4colTable_ed(const char *filename, const int iRap)
{
  ofstream of;
  of.open(filename, std::ios_base::out | std::ios_base::app);

  // of<<"% Energy Density profile for x=(" << Xmin <<", "<<Xmax
  //                                        << Ymin <<", "<<Ymax
  //                                        <<"Rapidity "<<rapMin+ iRap* nRap
  //                                        <<"Matching time "<<Tauf-Taui<<" fm/c"<<endl;                                   


    for(int i=0;i<Maxx;i++)      
    {
      for(int j=0;j<Maxy;j++)
      of  << Tauf-Taui<<"   "
          << Xmin+i*dx<<"   "
          << Ymin+j*dy<<"   "
          << scientific << setprecision(10) << setw(19) << DataTable->GetEd(iRap, i, j)<<endl;  //need revise
      
    }

  cout<<"Energy density table has been printed out for time Delta_Tau= "<<Tauf-Taui<<" fm/c"<<endl<<endl;
  of.close();
}


void LdMatching::Output4colTable_visratio(const char *filename, const int iRap)
{
  ofstream of;
  of.open(filename, std::ios_base::out | std::ios_base::app);                                  

    for(int i=0;i<Maxx;i++)      
    {
      for(int j=0;j<Maxy;j++)
      of  << Tauf-Taui<<"   "
          << Xmin+i*dx<<"   "
          << Ymin+j*dy<<"   "
          << scientific << setprecision(10) << setw(19) << DataTable->GetVisRatio(iRap, i, j)<<endl;  //need revise
      
    }

  cout<<"vis/ideal ratio table has been printed out for time Delta_Tau= "<<Tauf-Taui<<" fm/c"<<endl<<endl;
  of.close();
}

void LdMatching::OutputTables_Pimn(const int iRap)
{
  cout <<"Start to dump Pi_mn table"<<endl;
  ofstream of_pi00, of_pi01, of_pi02;
  ofstream of_pi11, of_pi12, of_pi22;
  ofstream of_pi33;

  of_pi00.open("data/Pi00_kln_stred.dat", std::ios_base::out);
  of_pi01.open("data/Pi01_kln_stred.dat", std::ios_base::out);
  of_pi02.open("data/Pi02_kln_stred.dat", std::ios_base::out);
  of_pi11.open("data/Pi11_kln_stred.dat", std::ios_base::out);
  of_pi12.open("data/Pi12_kln_stred.dat", std::ios_base::out);
  of_pi22.open("data/Pi22_kln_stred.dat", std::ios_base::out);
  of_pi33.open("data/Pi33_kln_stred.dat", std::ios_base::out);

    for(int i=0;i<Maxx;i++)      
    {
      for(int j=0;j<Maxy;j++)
      {
        of_pi00 << setprecision(12) << setw(22) << DataTable->GetPi_mn(iRap,i, j, 0, 0); 
        of_pi01 << setprecision(12) << setw(22) << DataTable->GetPi_mn(iRap,i, j, 0, 1); 
        of_pi02 << setprecision(12) << setw(22) << DataTable->GetPi_mn(iRap,i, j, 0, 2); 
        of_pi11 << setprecision(12) << setw(22) << DataTable->GetPi_mn(iRap,i, j, 1, 1); 
        of_pi12 << setprecision(12) << setw(22) << DataTable->GetPi_mn(iRap,i, j, 1, 2); 
        of_pi22 << setprecision(12) << setw(22) << DataTable->GetPi_mn(iRap,i, j, 2, 2); 
        of_pi33 << setprecision(12) << setw(22) << DataTable->GetPi_mn(iRap,i, j, 3, 3); 
      }
        of_pi00 << endl;
        of_pi01 << endl;
        of_pi02 << endl;
        of_pi11 << endl;
        of_pi12 << endl;
        of_pi22 << endl;
        of_pi33 << endl;
    }
  of_pi00.close();
  of_pi01.close();
  of_pi02.close();
  of_pi11.close();
  of_pi12.close();
  of_pi22.close();
  of_pi33.close();

  cout << "Pi_mn table dumped!"<<endl;
}


void LdMatching::OutputTable_Sd(const char *filename, const int iRap)
{
  ofstream of;
  of.open(filename, std::ios_base::out);

  of<<"% Entropy Density profile for x=(" << Xmin <<", "<<Xmax
                                         << Ymin <<", "<<Ymax
                                         <<"Rapidity "<<rapMin+ iRap* nRap
                                         <<"Matching time "<<Tauf-Taui<<" fm/c"<<endl;                                   

    for(int i=0;i<Maxx;i++)      
    {
      for(int j=0;j<Maxy;j++)
      of  << setprecision(10) << setw(22) << DataTable->GetSd(iRap, i, j);  //need revise
      of  << endl;
    }

  cout<<"Entropy density table has been printed out!"<<endl;
  of.close();
}


void LdMatching::OutputTable_ed(const char *filename, const int iRap)
{
  ofstream of;
  of.open(filename, std::ios_base::out);

  of<<"% Energy Density profile for x=(" << Xmin <<", "<<Xmax
                                         << Ymin <<", "<<Ymax
                                         <<"Rapidity "<<rapMin+ iRap* nRap
                                         <<"Matching time "<<Tauf-Taui<<" fm/c"<<endl;                                   

    for(int i=0;i<Maxx;i++)      
    {
      for(int j=0;j<Maxy;j++)
      of  << setprecision(12) << setw(22) << DataTable->GetEd(iRap, i, j);  //need revise
      of  << endl;
    }

  cout<<"Energy density table has been printed out!"<<endl;
  of.close();
}

void LdMatching::OutputTable_BulkPi(const char *filename, const int iRap)
{
  ofstream of;
  of.open(filename, std::ios_base::out);

  of<<"% Bulk viscosity profile for x=(" << Xmin <<", "<<Xmax
                                         << Ymin <<", "<<Ymax
                                         <<"Rapidity "<<rapMin+ iRap* nRap
                                         <<"Matching time "<<Tauf-Taui<<" fm/c"<<endl;                                   

    for(int i=0;i<Maxx;i++)      
    {
      for(int j=0;j<Maxy;j++)
      of  << setprecision(10) << setw(22) << DataTable->GetBulk_Pi(iRap, i, j);  //need revise
      of  << endl;
    }

  cout<<"Bulk Viscosity table has been printed out!"<<endl;
  of.close();
}


//debug
void LdMatching::Output_picontract_comp(const char *filename, const int iRap)
{
  ofstream of;
  of.open(filename, std::ios_base::out);

  of<<"% Viscosity ratio profile for x=(" << Xmin <<", "<<Xmax
                                         << Ymin <<", "<<Ymax
                                         <<"Rapidity "<<rapMin+ iRap* nRap
                                         <<"Matching time "<<Tauf-Taui<<" fm/c"<<endl;                                   

    for(int i=0;i<Maxx;i++)      
    {
      for(int j=0;j<Maxy;j++)
      of  << setprecision(12) << setw(22) << DataTable->GetPiSq(iRap, i, j);  //need revise
      of  << endl;
    }

  cout<<"Vis/Ideal Ratio table has been printed out for comparison!"<<endl;
  of.close();
}

void LdMatching::OutputTmnTable(const char *filename,const int iRap)
{
  ofstream of;
  of.open(filename, std::ios_base::out);

//  Output Tmn table
    for(int i=0;i<Maxx;i++)
      {
         for(int j=0;j<Maxy;j++)
            {
              of <<setprecision(12) << setw(22)<<DataTable->GetTmn(iRap, i, j, 0, 0);
//                  <<setprecision(12) << setw(22)<<T01[iRap][i][j]
//                  <<setprecision(12) << setw(22)<<T02[iRap][i][j]
//                  <<setprecision(12) << setw(22)<<T11[iRap][i][j]
//                  <<setprecision(12) << setw(22)<<T12[iRap][i][j]
//                  <<setprecision(12) << setw(22)<<T22[iRap][i][j]<<endl;
            }
        of << endl;

      }
  cout<<"T^00 has been printed out!"<<endl;
  of.close();

} 


void LdMatching::Lower_matrix_single(gsl_matrix *dest, gsl_matrix *src)
{
// This routine can lower one index of a tensor: T^{\mu\nu} ---->T^{\mu\alpha}*g_{\alpha\nu}----->T^{\mu}_{\nu}
  double data[]={1, -1, -1, -1,
                 1, -1, -1, -1,
                 1, -1, -1, -1,
                 1, -1, -1, -1};
  gsl_matrix_view gmn_converter 
       = gsl_matrix_view_array (data, 4, 4);
  gsl_matrix_memcpy(dest, src);
  gsl_matrix_mul_elements(dest, &gmn_converter.matrix);  //Equivalent to Tmn*gmn
}

void LdMatching::Lower_matrix_double(gsl_matrix *dest, gsl_matrix *src)
{
// This routine can lower two indices of a tensor: T^{\mu\nu} ---->g_{\mu\alpha}* T^{\alpha\beta}*g_{\beta\nu} = T^{\mu}_{\nu}
  double data[]={1, -1, -1, -1,
                -1,  1,  1,  1,
                -1,  1,  1,  1,
                -1,  1,  1,  1};
  gsl_matrix_view gmn_converter2 
       = gsl_matrix_view_array (data, 4, 4);
  gsl_matrix_memcpy(dest, src);
  gsl_matrix_mul_elements(dest, &gmn_converter2.matrix);  //Equivalent to Tmn*gmn
}

double LdMatching::Contract_matrix(gsl_matrix *upper)
{
// Input the upper matrix T^mn, this routine returns the contraction T^mn * T_mn, the matrix should be 4 dimensional
  double tr=0;
  double a[] = { 0.00, 0.00, 0.00, 0.00,
                 0.00, 0.00, 0.00, 0.00,
                 0.00, 0.00, 0.00, 0.00,
                 0.00, 0.00, 0.00, 0.00 }; 

  gsl_matrix_view result = gsl_matrix_view_array(a, 4, 4);

  gsl_matrix * lower = gsl_matrix_alloc (4, 4);
  Lower_matrix_double(lower, upper);

  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                 1.0, upper, lower,
                 0.0, &result.matrix);

  for(int k=0;k<4;k++)
    tr+=gsl_matrix_get(&result.matrix,k,k);

  gsl_matrix_free(lower);
  return tr;

}

void LdMatching::Diagnostic(int iRap, int i, int j)
{
//output the various quantities stored in the current point
  cout<<"Initial diagonosis at the current point: "
                 <<setw(22)<<setprecision(10)<<"i= "<<i<<" j="<<j<<endl;
  cout<<"At this point: "<<Xmin+i*dx<<" "<<Ymin+j*dy<<endl;
  cout<<"Tmn at this point: "<<endl;
  gsl_matrix *Tmn_current = gsl_matrix_alloc(4,4);
  DataTable->Tmn_matrix_get(Tmn_current, iRap, i, j);
  for(int ip=0;ip<4;ip++)
  {
    for(int jp=0;jp<4;jp++)
      cout<<scientific<<setw(18)<<setprecision(10)<<gsl_matrix_get(Tmn_current, ip, jp);
    cout<<endl;
  }

  cout<<"U_mu at this point: "<<endl;
  for(int ip=0;ip<4;ip++)
    cout<<scientific<<setw(18)<<setprecision(10)<<DataTable->GetUm(iRap, i, j, ip);
  cout<<endl;

  cout<<"Energy density at this point: "<<endl
      <<scientific<<setw(18)<<setprecision(10)<<DataTable->GetEd(iRap, i, j)<<endl;
  cout<<"Bulk viscosity at this point: "<<endl
      <<scientific<<setw(18)<<setprecision(10)<<DataTable->GetBulk_Pi(iRap, i, j)<<endl; 

  cout<<"Pi_mn at this point:"<<endl;    
  for(int ir=0;ir<4;ir++)
  {
    for(int jr=0;jr<4;jr++)
      cout<<scientific<<setw(18)<<setprecision(10)<<DataTable->GetPi_mn(iRap, i, j, ir, jr);
    cout<<endl;
  }
  
  gsl_matrix_free(Tmn_current);
  exit(0);
  
};
