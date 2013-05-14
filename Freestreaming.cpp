#include <cmath>
#include <iomanip>
#include <vector>
#include <sstream>
#include "Freestreaming.h"
#include "gauss_quadrature.h"
#include "mistools.h"
#include "arsenal.h"

using namespace std;

FreeStrm::FreeStrm(double xmax, double ymax, double ptmin, double ptmax, double dx0,double dy0,
		    double dpt0, int ny, double rapmin, double rapmax, double taui, double tauf)
{
    Xmax=xmax;
    Ymax=ymax;
    Xmin=-xmax;
    Ymin=-ymax;
    
    PTmax=ptmax;
    PTmin=ptmin;
    dpt=dpt0;
    MaxPT=(int)((PTmax-PTmin)/dpt+0.1)+1;

    dx=dx0;
    dy=dy0;
    Maxx=(int)((Xmax-Xmin)/dx+0.1)+1;
    Maxy=(int)((Ymax-Ymin)/dy+0.1)+1;
    nRap=ny;
    rapMin=rapmin;
    rapMax=rapmax;

    Taui=taui;
    Tauf=tauf;
    
    Xcm = 0.;
    Ycm = 0.;
//density with pt dependence will be read from 
    densityTable  = new double*** [nRap];
    densityTableInterp  = new double*** [nRap];
    for(int iy=0;iy<nRap;iy++) {
    densityTable[iy] =  new double** [Maxx];
    densityTableInterp[iy] =  new double** [Maxx];
    for(int i=0;i<Maxx;i++) {
        densityTable[iy][i] = new double* [Maxy];
        densityTableInterp[iy][i] = new double* [Maxy];
        for(int j=0;j<Maxy;j++) {
            densityTable[iy][i][j]=new double[MaxPT];
                      for(int ipt=0;ipt<MaxPT; ipt++)
                            densityTable[iy][i][j][ipt]=0.0;
        }
      }
    } 

// Free-streamed density profile: in a certain direction
    shiftedTable  = new double*** [nRap];
    for(int iy=0;iy<nRap;iy++) {
    shiftedTable[iy] =  new double** [Maxx];
    for(int i=0;i<Maxx;i++) {
        shiftedTable[iy][i] = new double* [Maxy];
        for(int j=0;j<Maxy;j++) {
            shiftedTable[iy][i][j]=new double[MaxPT];
                      for(int ipt=0;ipt<MaxPT; ipt++)
                            shiftedTable[iy][i][j][ipt]=0.0;
        }
      }
    }

 //vector Pt integrated gluon density profile
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

 //vector Pt integrated gluon density profile
    dEd2rdyTable  = new double** [nRap];
    for(int iy=0;iy<nRap;iy++) {
    dEd2rdyTable[iy] =  new double* [Maxx];
    for(int i=0;i<Maxx;i++) {
        dEd2rdyTable[iy][i] = new double [Maxy];
        for(int j=0;j<Maxy;j++) {
            dEd2rdyTable[iy][i][j]=0.0;
        }
      }
    }
    cout << "Free-streaming procedure Initialized!" << endl;

}


FreeStrm::~FreeStrm()
{
//clean densityTable
    for(int iy=0;iy<nRap;iy++) {
      for(int i=0;i<Maxx;i++) {
        for(int j=0;j<Maxy;j++) {
           delete [] densityTable[iy][i][j];
           delete [] densityTableInterp[iy][i][j];
        }
        delete [] densityTable[iy][i];
        delete [] densityTableInterp[iy][i];
        }
    delete [] densityTable[iy];
    delete [] densityTableInterp[iy];
    }
    delete [] densityTable;
    delete [] densityTableInterp;

//Clean shifted table
    for(int iy=0;iy<nRap;iy++) {
      for(int i=0;i<Maxx;i++) {
        for(int j=0;j<Maxy;j++) delete [] shiftedTable[iy][i][j];
        delete [] shiftedTable[iy][i];
        }
    delete [] shiftedTable[iy];
    }
    delete [] shiftedTable;   

//Clean pt integrated table
    for(int iy=0;iy<nRap;iy++)  {
        {
          for(int i=0;i<Maxx;i++) 
            delete [] dNd2rdyTable[iy][i];
        }
      delete [] dNd2rdyTable[iy];
      }
    delete [] dNd2rdyTable;
//Clean pt integrated table
    for(int iy=0;iy<nRap;iy++)  {
        {
          for(int i=0;i<Maxx;i++) 
            delete [] dEd2rdyTable[iy][i];
        }
      delete [] dEd2rdyTable[iy];
      }
    delete [] dEd2rdyTable;
  

}


void FreeStrm::ReadTable(string filename)
{
  double a,b,c,d;    //temp, ditch them
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
              for(int ipt=0;ipt<MaxPT;ipt++)
                  {
                    DataFile>>densityTable[iy][i][j][ipt];
                  }
    DataFile >> a >> b;  //safty procedure, reach the end of the file
    cout << "Read in complete!" << endl;
  }
  DataFile.close();
  cout << "Data table has been read into memory!" << endl;

}


void FreeStrm::GetDensity(int iRap, int i, int j, double phip, double* results)
{
    double x,y;     //unshifted coordinate
    double shfedx, shfedy;    //shifed coordinate
    double DTau, Phip;
    double it, jt;             //indices corresponding to shifted coordinates

    Phip=phip;

    x=Xmin+i*dx;
    y=Ymin+j*dy;
    DTau=Tauf-Taui;

    shfedx=x-DTau*cos(Phip);
    shfedy=y-DTau*sin(Phip);

    double *A0 = new double [MaxPT];
    double *A1 = new double [MaxPT];
    double *A2 = new double [MaxPT];
    double *A3 = new double [MaxPT];

    if( (shfedx>Xmax||shfedx<Xmin) ||(shfedy>Ymax||shfedy<Ymin) )  //coordinates locate outside of the grid
    {
        for(int ipt = 0; ipt < MaxPT; ipt++)
           results[i] = 0.0;
    }
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

        interpCubic4PointsArray(densityTable[iRap][iti][jti], densityTable[iRap][iti][jti+1], 
          densityTable[iRap][iti][jti+2], densityTable[iRap][iti][jti+3], MaxPT, 1, yfraction, A0);

        interpCubic4PointsArray(densityTable[iRap][iti+1][jti], densityTable[iRap][iti+1][jti+1], 
          densityTable[iRap][iti+1][jti+2], densityTable[iRap][iti+1][jti+3], MaxPT, 1, yfraction, A1);

        interpCubic4PointsArray(densityTable[iRap][iti+2][jti], densityTable[iRap][iti+2][jti+1], 
          densityTable[iRap][iti+2][jti+2], densityTable[iRap][iti+2][jti+3], MaxPT, 1, yfraction, A2);

        interpCubic4PointsArray(densityTable[iRap][iti+3][jti], densityTable[iRap][iti+3][jti+1], 
          densityTable[iRap][iti+3][jti+2], densityTable[iRap][iti+3][jti+3], MaxPT, 1, yfraction, A3);

       // row interpolation + extrapolation first
//        if(A0<0||A1<0||A2<0||A3<0) 
//        {
//          cout<<"hhhaaahaha"<<endl;
//          cout<<"interpCubic4Points messed up, result="<<endl
//              <<A0<<" "<<A1<<" "<<A2<<" "<<A3<<endl; 	 
//         exit(0);
//        }  //debug
       interpCubic4PointsArray(A0,A1,A2,A3, MaxPT, 1, xfraction, results);
    }

    delete [] A0;
    delete [] A1;
    delete [] A2;
    delete [] A3;
    return;


    //revised in Mar.18, 2013
    //switch to linear interpolation since:
    //1. cubic interpolation assumes smooth function, which is unknown for fKLN output
    //2. negative result is given by cubic method
    //interp on boundary
//         if (jti<0) jti=0;
//         if (jti>=Maxy-2) 
//           {
//             jti=Maxy-2; // need 2 points
//             // cout<<"out of y boundary"<<endl;
//           }
// 
//         if (iti<0) iti=0;       
//         if (iti>=Maxx-2) 
//           { 
//             iti=Maxx-2; // need 2 points
//            // cout<<"out of x boundary"<<endl;
//           }
//  	return Bilinear2dInterp(iti, jti, 1, 1, densityTable[iRap][iti][jti][ipt], 
//                    densityTable[iRap][iti][jti+1][ipt], densityTable[iRap][iti+1][jti+1][ipt], 
//                    densityTable[iRap][iti+1][jti][ipt]);

}

void FreeStrm::GetDensityInterp(int iRap, int i, int j, double phip, double* results)
{
    double x,y;     //unshifted coordinate
    double shfedx, shfedy;    //shifed coordinate
    double DTau, Phip;
    double it, jt;             //indices corresponding to shifted coordinates

    Phip=phip;

    x=Xmin+i*dx;
    y=Ymin+j*dy;
    DTau=Tauf-Taui;

    shfedx=x-DTau*cos(Phip);
    shfedy=y-DTau*sin(Phip);

    double *A0 = new double [MaxPT];
    double *A1 = new double [MaxPT];
    double *A2 = new double [MaxPT];
    double *A3 = new double [MaxPT];

    if( (shfedx>Xmax||shfedx<Xmin) ||(shfedy>Ymax||shfedy<Ymin) )  //coordinates locate outside of the grid
    {
        for(int ipt = 0; ipt < MaxPT; ipt++)
           results[i] = 0.0;
    }
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

        interpCubic4PointsArray(densityTableInterp[iRap][iti][jti], densityTableInterp[iRap][iti][jti+1], 
          densityTableInterp[iRap][iti][jti+2], densityTableInterp[iRap][iti][jti+3], NpTinterp, 1, yfraction, A0);

        interpCubic4PointsArray(densityTableInterp[iRap][iti+1][jti], densityTableInterp[iRap][iti+1][jti+1], 
          densityTableInterp[iRap][iti+1][jti+2], densityTableInterp[iRap][iti+1][jti+3], NpTinterp, 1, yfraction, A1);

        interpCubic4PointsArray(densityTableInterp[iRap][iti+2][jti], densityTableInterp[iRap][iti+2][jti+1], 
          densityTableInterp[iRap][iti+2][jti+2], densityTableInterp[iRap][iti+2][jti+3], NpTinterp, 1, yfraction, A2);

        interpCubic4PointsArray(densityTableInterp[iRap][iti+3][jti], densityTableInterp[iRap][iti+3][jti+1], 
          densityTableInterp[iRap][iti+3][jti+2], densityTableInterp[iRap][iti+3][jti+3], NpTinterp, 1, yfraction, A3);

       interpCubic4PointsArray(A0,A1,A2,A3, MaxPT, 1, xfraction, results);
    }

    delete [] A0;
    delete [] A1;
    delete [] A2;
    delete [] A3;
    return;
}

void FreeStrm::ShiftDensityInterp(const int iRap, double phip, double ****shiftedTableInterp)
{

  double Phip=phip;
  double* local_results = new double [NpTinterp];

  for(int i=0;i<Maxx;i++)
  for(int j=0;j<Maxy;j++)
  {
    GetDensityInterp(iRap, i , j, Phip, local_results);
    for(int ipt = 0; ipt < NpTinterp; ipt++)
       shiftedTableInterp[iRap][i][j][ipt]= local_results[ipt];
  }
  delete [] local_results;
//  cout<<"Profile has been shifted to: "<<Phip<<endl;
}

void FreeStrm::ShiftDensity(const int iRap, double phip)
{

  double Phip=phip;
  double* local_results = new double [MaxPT];

  for(int i=0;i<Maxx;i++)
  for(int j=0;j<Maxy;j++)
  {
    GetDensity(iRap, i , j, Phip, local_results);
    for(int ipt = 0; ipt < MaxPT; ipt++)
       shiftedTable[iRap][i][j][ipt]= local_results[ipt];
  }
  delete [] local_results;
//  cout<<"Profile has been shifted to: "<<Phip<<endl;
}


void FreeStrm::InterpDensity(const int iRap, double* pT, int npT)
{
   NpTinterp = npT;
   vector<double>* dens1=new vector<double>(MaxPT,0.0);
   vector<double>* pt0=new vector<double>(MaxPT,0.0);    
   for(int ipt0 = 0; ipt0 < MaxPT; ipt0++)
      (*pt0)[ipt0]= (PTmin + dpt*ipt0);

   for(int iy=0;iy<nRap;iy++)
      for(int i=0;i<Maxx;i++)
         for(int j=0;j<Maxy;j++)
         {
            densityTableInterp[iy][i][j]=new double[npT];
            for(int ipt = 0; ipt < npT; ipt++)
               densityTableInterp[iy][i][j][ipt]=0.0;
         }

   for(int i = 0; i < Maxx; i++)
      for(int j = 0; j < Maxy; j++)
      {
         for(int ipt0 = 0; ipt0 < MaxPT; ipt0++)
            (*dens1)[ipt0] = densityTable[iRap][i][j][ipt0];
         for(int ipt = 0; ipt < npT; ipt++)
            densityTableInterp[iRap][i][j][ipt]=interpCubicDirect(pt0, dens1, pT[ipt]);
      }
   delete dens1;
   delete pt0;
}

void FreeStrm::OutputTable(const char *filename, const int iRap)
{
  ofstream of;
  of.open(filename, std::ios_base::out);

  double rap = rapMin+(rapMax-rapMin)/nRap*iRap;


//  Output dNd2rdy Table
    for(int i=0;i<Maxx;i++)
      for(int j=0;j<Maxy;j++)
        for(int ipt=0;ipt<MaxPT;ipt++)
    {
      of <<  setprecision(3) << setw(10) <<  rap 
          << setprecision(3) << setw(10) <<  Xmin+i*dx 
          << setprecision(3) << setw(10) <<  Ymin+j*dy
          << setprecision(3) << setw(10) <<  PTmin+ipt*dpt
          << setprecision(12) << setw(22) << dNd2rdyTable[iRap][i][j]
          << endl;
    }

  cout<<"Free Streamed Profile complete!"<<endl;
  of.close();
}




void FreeStrm::dumpBlockTable(const char *filename, double ***data, const int iRap)
{
  
  ofstream of;
  of.open(filename, std::ios_base::out);

    for(int i=0;i<Maxx;i++)
    {
       for(int j=0;j<Maxy;j++)
           {
               of <<scientific << setw(22) <<data[iRap][i][j];              
           }
    of << endl;
    }

  cout<<"Block Table complete!"<<endl;
  of.close();

}

//Following functions are used for debugging and testing if free streaming
//generate the desirable profile.
void FreeStrm::CreateDataTable(const char *filename, const int iRap)
{
  ofstream of;
  of.open(filename, std::ios_base::out);
  cout << "Work in test mode: using Gaussian Profile!" << endl;

  double rap = rapMin+(rapMax-rapMin)/nRap*iRap;


    for(int i=0;i<Maxx;i++)
    for(int j=0;j<Maxy;j++)
    for(int ipt=0;ipt<MaxPT;ipt++)
    {
      double x=Xmin+i*dx;
      double y=Ymin+j*dy;
      double ptstep = PTmin + ipt*dpt;
      of <<  setprecision(3) << setw(10) <<  rap 
          << setprecision(3) << setw(10) <<  x 
          << setprecision(3) << setw(10) <<  y 
          << setprecision(3) << setw(10) <<  ptstep 
          << setprecision(12) << setw(22) << GaussProfile(iRap,i,j,ipt)   //debugging
          << endl;
    }
  cout<<"Gaussian initial condition is created!"<<endl;
  of.close();
}


double FreeStrm::GaussProfile(const int iRap, int i, int j, int ipt)
{
    double x=Xmin+i*dx;
    double y=Ymin+j*dy;
    double ptstep = PTmin + ipt*dpt;

    double fxy;
    // prefactor=1/(12.087747372898045e0), Gaussian normalization factor

    fxy= exp(-(x - 0.1)*(x - 0.1)/3.0 - (y - 0.1)*(y - 0.1)/3.0 - ptstep*ptstep)/12.087747372898045;


    return fxy;

}

double FreeStrm::BoxProfile(const int iRap, int i, int j, int ipt)
{   
    double x=Xmin+i*dx;
    double y=Ymin+j*dy;
    double ptstep = PTmin + ipt*dpt;
    double fxy;
    double prefactor=1/18.0;

    fxy= prefactor*(stepfunc(x+3)-stepfunc(x-3))*(stepfunc(y+3)-stepfunc(y-3))*exp(-ptstep*ptstep);

    return fxy;
}


void FreeStrm::generateEdTable(const int iRap)
//Multiply free-streamed gluon density by pt 
//then integrate over Pt and Phi_pt to get dE/d^2rdy spectra
{
//block for gaussian integration
  cout << "Freestreaming: Start to generate dE/dy Table---"<<endl;
  ostringstream ed_filename;
  ed_filename.str("");  //clean the string stream
  int kind=1;
  const int order=100;  //debug
  double alpha=0.0, beta=0.0;
  double xphip[order],wphi[order];
  double xpt[order],wpt[order];

  double delta_tau = Tauf-Taui;
  if(Tauf == Taui)
  {
  delta_tau = 1;
  cout<<"Notice: Energy density is calculated at tauf=taui="<<Tauf<<endl;
  }
  double phipmin=0.0, phipmax=2.0*M_PI;

//initiate intermediate produced gluon density table for a certain phi_p
  double ****dEd2rdyPhipTable;
    dEd2rdyPhipTable  = new double*** [nRap];
    for(int iy=0;iy<nRap;iy++) {
    dEd2rdyPhipTable[iy] =  new double** [Maxx];
    for(int i=0;i<Maxx;i++) {
        dEd2rdyPhipTable[iy][i] = new double* [Maxy];
        for(int j=0;j<Maxy;j++) {
            dEd2rdyPhipTable[iy][i][j]=new double[order];
                      for(int iphip=0;iphip<order; iphip++)
                            dEd2rdyPhipTable[iy][i][j][iphip]=0.0;
        }
      }
    } 

  gauss_quadrature(order, kind, alpha, beta, phipmin, phipmax,xphip, wphi); 
  gauss_quadrature(order, kind, alpha, beta, PTmin, PTmax, xpt, wpt); 

  for(int iphi=0;iphi<order;iphi++)
    {
      ShiftDensity(0, xphip[iphi]);

//Initialize vectors for cubic interpolation
      vector<double>* dens1=new vector<double>(MaxPT,0.0);
      vector<double>* pt0=new vector<double>(MaxPT,0.0);    
      for(int ipt0=0;ipt0<MaxPT;ipt0++)
        (*pt0)[ipt0]=(PTmin+dpt*ipt0);

//inner integration: integrate over pt for a certain phi_p.
      for(int i=0;i<Maxx;i++)  {
        for(int j=0;j<Maxy;j++) {
          for(int ipt=0;ipt<order;ipt++)
          {

            for(int k=0;k<MaxPT;k++)
            {
              (*dens1)[k]=shiftedTable[iRap][i][j][k];
            } 

//Use 1D cubic interpolation interpCubicDirect(vector<double>* x, vector<double>* y, double x0)

          double elem=interpCubicDirect(pt0, dens1, xpt[ipt]);

//use the following line to get dE/dy table at tau=tauf  /(Tauf-Taui)
          dEd2rdyPhipTable[iRap][i][j][iphi]+=elem*wpt[ipt]*xpt[ipt]*xpt[ipt];
         }

       }

      } //for int i=0;i<Maxx;i++



    for(int iy=0;iy<nRap;iy++)
      for(int i=0;i<Maxx;i++)
        for(int j=0;j<Maxy;j++)
          for(int iphi=0;iphi<order;iphi++) {
            dEd2rdyTable[iRap][i][j]+=dEd2rdyPhipTable[iRap][i][j][iphi]*wphi[iphi]/delta_tau;
            dEd2rdyPhipTable[iRap][i][j][iphi]=0.0;    //reset table for next loop
          }  
  //clean up vectors
    delete  dens1;
    delete  pt0;
  } //for int iphi=0;iphi<order;iphi++


  //Clean intermediate table
    for(int iy=0;iy<nRap;iy++) {
      for(int i=0;i<Maxx;i++) {
        for(int j=0;j<Maxy;j++) delete [] dEd2rdyPhipTable[iy][i][j];
        delete [] dEd2rdyPhipTable[iy][i];
        }
    delete [] dEd2rdyPhipTable[iy];
    }
    delete [] dEd2rdyPhipTable; 
  //dump energy density table
  ed_filename<<"data/ed_fs"<<"_tauf_"<<Tauf<<".dat";
  dumpBlockTable(ed_filename.str().c_str(), dEd2rdyTable, 0);

  cout << "Energy density profile has been dumped to: "<< ed_filename.str() 
       <<endl;

  findCM(dEd2rdyTable);

}

void FreeStrm::findCM(double ***source, const int iRap)
{
/*Find the principal axis of the elliptic profile and align it with the x-axis
// {x} = \int{dx * dy * e(x,y) * x} / \int{dx * dy * e(x, y)}
// {y} = \int{dx * dy * e(x,y) * y} / \int{dx * dy * e(x, y)}
*/
  double x_ave = 0.;
  double y_ave = 0.;
  double ed_total = 0.;

  for(int i = 0; i < Maxx; i++)
  {
    for(int j = 0; j< Maxy; j++)
    {
      double x_current = Xmin + dx * i;
      double y_current = Ymin + dy * j;

      x_ave += source[iRap][i][j] * x_current * dx * dy;
      y_ave += source[iRap][i][j] * y_current * dx * dy;
      ed_total += source[iRap][i][j] * dx * dy;
    }
  } //<-> for i=0:Maxx

  x_ave = x_ave/(ed_total + 1e-18);
  y_ave = y_ave/(ed_total + 1e-18);

  Xcm = x_ave;
  Ycm = y_ave;

  cout << "Weighted center of profile is: " << endl
       << "("
       << Xcm << ", " << Ycm 
       << ")" <<endl;  //debug
}


double FreeStrm::getShiftedProfile(double ***data, int i, int j, double x0, double y0, 
  double phi, bool limit, const int iRap)
{
    double x,y;     //unshifted coordinate
    double shfedx, shfedy;    //shifed coordinate
    double Phi;
    double it, jt;             //indices corresponding to shifted coordinates

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
}

// double FreeStrm::getEpx(const int iRap)
// //use gaussian quadrature integration routine to get spatial eccentricity
// {
//   cout<<"Calculating spatial eccentricity Epx-------"<<endl;
//   double Epx=0.0;
//   double Epx_angle = 0.0;
//   double Epx_nu = 0.0, Epx_dn = 0.0;    //Epx = Epx_nu/Epx_dn
//   const int method= 2;  // 1 for gaussian quadrature integration
//                         // 2 for direct summing over x and y
// //block for gaussian integration
//   if(method ==1) 
//   {
//     int kind=1;
//     const int order=100;
//     double alpha=0.0, beta=0.0;
//     double x_pts[order],x_weight[order];
//     double y_pts[order],y_weight[order];

// //Assign gaussian points
//     gauss_quadrature(order, kind, alpha, beta, Xmin, Xmax, x_pts, x_weight); 
//     gauss_quadrature(order, kind, alpha, beta, Ymin, Ymax, y_pts, y_weight);   

//   //initiate intermediate produced gluon density table for a certain x
//     vector<double>* dEpxdxTable=new vector<double>(Maxx,0.0);
//     vector<double>* dEpxdxTable_dn=new vector<double>(Maxx,0.0);  //for denominator
  
//   //do the integration on y to generate a dEpx/dy lookup table
//     vector<double>* y_grids=new vector<double>(Maxy,0.0);
//     vector<double>* dEd2rdy_grids=new vector<double>(Maxy,0.0); 
//     for(int i=0;i<Maxx;i++)
//     {
//   //Initialize vectors for cubic interpolation of energy density on x-axis
//       for(int iy0=0;iy0<Maxy;iy0++)
//       {
//         (*y_grids)[iy0]=(Ymin+dy*iy0);
//         (*dEd2rdy_grids)[iy0]=dEd2rdyTable[iRap][i][iy0];  //debug, arrive here
//       }
//   //Use 1D cubic interpolation interpCubicDirect(vector<double>* x, vector<double>* y, double x0)
//       for(int iy_order=0;iy_order<order;iy_order++)
//       {
//         double elem=interpCubicDirect(y_grids, dEd2rdy_grids, y_pts[iy_order]);
//         double x_current = Xmin + dx*i;
//         double y_current = y_pts[iy_order];
//         double epx_numerator = y_current*y_current-x_current*x_current;
//         double epx_denominator = y_current*y_current+x_current*x_current;
//         (*dEpxdxTable)[i] += elem*y_weight[iy_order]*epx_numerator;   // Integ[(y^2-x^2)*e(x,y),{y, ymin, ymax}]
//         (*dEpxdxTable_dn)[i] += elem*y_weight[iy_order]*epx_denominator;   // Integ[(x^2+y^2)*e(x,y),{y, ymin, ymax}]
//       }
//     }//<->for i=0:Maxx
  
  
//   //Integrate over x
//     vector<double>* x_grids=new vector<double>(Maxx,0.0);
//     for(int ix0=0;ix0<Maxx;ix0++)
//       (*x_grids)[ix0]=(Ymin+dx*ix0);
  
//     for(int ix_order=0;ix_order<order;ix_order++)
//     {
//       double interp_dEpx = interpCubicDirect(x_grids, dEpxdxTable, x_pts[ix_order]);
//       double interp_dEpx_dn = interpCubicDirect(x_grids, dEpxdxTable_dn, x_pts[ix_order]);
//       Epx_nu += interp_dEpx * x_weight[ix_order];
//       Epx_dn += interp_dEpx_dn * x_weight[ix_order];
//     } 
  
//     Epx = Epx_nu/(Epx_dn+1e-18);
//   }  

//   if(method ==2)  //direct sum over x and y
//   {
//     double epx_nu_real = 0.0;  //numerator of epx: real part
//     double epx_nu_img = 0.0; //numerator of epx: imaginary part
//     double epx_dn = 0.0;  //denominator of epx

//     // //prepare file name
//     // ostringstream filename_stream;
//     // filename_stream.str("");  //clean up before use it
//     // filename_stream << "data/ed_shifted_tauf_" << Tauf << ".dat";

//     // //open file to output shifted energy density for the current proper time
//     // ofstream of;
//     // of.open(filename_stream.str().c_str(), std::ios_base::out);

//     //start to calculate epx
//     for(int i=0; i<Maxx; i++)
//     {
//       for(int j=0; j<Maxy; j++)
//       {
//        double x = Xmin + dx*i;
//        double y = Ymin + dy*j;

//        double ed_shifted   = getShiftedProfile(dEd2rdyTable, i, j, Xcm, Ycm, 0, true, 0);
//    // of.close();
//        epx_nu_real+= ed_shifted * ( y * y - x * x) * dx * dy;  // real part
                     
//        epx_nu_img += ed_shifted * ( y * x) * dx * dy;  //imaginary part

//        epx_dn += ed_shifted * ( y * y
//                      + x * x) * dx * dy;  // denominator

//        // of << setw(14) << setprecision(8) << ed_shifted ;  //output for debugging
//       } 
//       // of << endl;     
//     }//<-> for i=0:Maxx


//   Epx = sqrt( epx_nu_real * epx_nu_real + 4.0 * epx_nu_img * epx_nu_img )
//        /(epx_dn + 1e-18);
//   Epx_angle = atan(2.0*epx_nu_img/(epx_nu_real + 1e-18));
//   cout << " phi_2 = " << Epx_angle << endl;
//   }//<-> if method ==2

//   else 
//   {
//     cout<<"Method = "<<method<<", no such method exists!"<<endl;
//     exit(0);
//   }
// //  cout<<"Numerator get!"<<endl<<endl;
//   return Epx;

// } 
  
double FreeStrm::getEpx(int nth_order, const int iRap)
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
       double r = sqrt(x*x + y*y);
       double phi = atan2(y,x);

       double ed_temp = dEd2rdyTable[iRap][i][j];

       epx_nu_real += pow(r, nth_order) 
                  * cos(nth_order * phi) * ed_temp * dx*dy;
       epx_nu_img += pow(r, nth_order) 
                  * sin(nth_order * phi) * ed_temp * dx*dy;
       epx_dn += ed_temp * pow(r, nth_order)* dx*dy;
     }       
  }//<-> for i=0:Maxx
  
  Epx = sqrt( epx_nu_real * epx_nu_real + epx_nu_img * epx_nu_img )
       /(epx_dn + 1e-18);
  Epx_angle = -atan2( epx_nu_img, (epx_nu_real + 1e-18))/nth_order;
/*
  cout<<"Spatial Eccentricity complete!"<<endl;
  cout << "Epx_angle =" <<Epx_angle << endl;*/

  return Epx;
}



void FreeStrm::GetInteDensity(const int iRap)
//integrate free-streamed gluon density profile over Pt and Phi_pt to get dNd^2rdy spectra
//can be further used in routines calculating T^{\mu \nu}
{

//block for gaussian integration
  int kind=1;
  const int order=100;
  double alpha=0.0, beta=0.0;
  double xphip[order],wphi[order];
  double xpt[order],wpt[order];
//  double elem;

  double phipmin=0.0, phipmax=2.0*M_PI;

//initiate intermediate produced gluon density table for a certain phi_p
  double ****dNd2rdyPhipTable;
    dNd2rdyPhipTable  = new double*** [nRap];
    for(int iy=0;iy<nRap;iy++) {
    dNd2rdyPhipTable[iy] =  new double** [Maxx];
    for(int i=0;i<Maxx;i++) {
        dNd2rdyPhipTable[iy][i] = new double* [Maxy];
        for(int j=0;j<Maxy;j++) {
            dNd2rdyPhipTable[iy][i][j]=new double[order];
                      for(int iphip=0;iphip<order; iphip++)
                            dNd2rdyPhipTable[iy][i][j][iphip]=0.0;
        }
      }
    } 

  gauss_quadrature(order, kind, alpha, beta, phipmin, phipmax,xphip, wphi); 
  gauss_quadrature(order, kind, alpha, beta, PTmin, PTmax, xpt, wpt); 

  for(int iphi=0;iphi<order;iphi++)
    {
      ShiftDensity(0, xphip[iphi]);

//Initialize vectors for cubic interpolation
      vector<double>* dens1=new vector<double>(MaxPT,0.0);
      vector<double>* pt0=new vector<double>(MaxPT,0.0);    
      for(int ipt0=0;ipt0<MaxPT;ipt0++)
        (*pt0)[ipt0]=(PTmin+dpt*ipt0);

//inner integration: integrate over pt for a certain phi_p.
      for(int i=0;i<Maxx;i++)  {
        for(int j=0;j<Maxy;j++) {
          for(int ipt=0;ipt<order;ipt++)
          {

            for(int k=0;k<MaxPT;k++)
            {
              (*dens1)[k]=shiftedTable[iRap][i][j][k];
            } 

//Use 1D cubic interpolation interpCubicDirect(vector<double>* x, vector<double>* y, double x0)

          double elem=interpCubicDirect(pt0, dens1, xpt[ipt]);
//           if(elem<0) {cout<<"elem = "<<elem<<endl; }
//This line gives dN/dy at tau0 table******************************************************
           dNd2rdyPhipTable[iRap][i][j][iphi]+=elem*wpt[ipt]*xpt[ipt];
//******************************************************************************************

//use the following line to get dE/dy table at tau=tauf  /(Tauf-Taui)
//            dNd2rdyPhipTable[iRap][i][j][iphi]+=elem*wpt[ipt]*xpt[ipt]*xpt[ipt];
//******************************************************************************************
         }

       }

      } //for int i=0;i<Maxx;i++



    for(int iy=0;iy<nRap;iy++)
      for(int i=0;i<Maxx;i++)
        for(int j=0;j<Maxy;j++)
          for(int iphi=0;iphi<order;iphi++) {
            dNd2rdyTable[iRap][i][j]+=dNd2rdyPhipTable[iRap][i][j][iphi]*wphi[iphi];
            dNd2rdyPhipTable[iRap][i][j][iphi]=0.0;    //reset table for next loop
          }    
  } //for int iphi=0;iphi<order;iphi++


  //Clean intermediate table
    for(int iy=0;iy<nRap;iy++) {
      for(int i=0;i<Maxx;i++) {
        for(int j=0;j<Maxy;j++) delete [] dNd2rdyPhipTable[iy][i][j];
        delete [] dNd2rdyPhipTable[iy][i];
        }
    delete [] dNd2rdyPhipTable[iy];
    }
    delete [] dNd2rdyPhipTable; 

  cout<<"Pt integrated table complete!"<<endl;

}
