#ifndef Freestreaming_h
#define Freestreaming_h

#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include "stdlib.h"
/*
revise history:
Apr.02, 2013
Add the function to calculate spatial eccentricity: firstly find the 

Mar.22, 2013
1. To maintain higher precision, use cubic interpolation again, but set negative 
   number to zero by hand;
2. Add dEdy table, function generateEdTable() generateEpx() to calculate eccentricity in free-streaming
stage

Mar.18, 2013, 
In function switch to linear interpolation to avoid negative value returned by cubic method
the output file of both methods have been compaired:
relative error 2.76%, difference of total entropy densities -2.7345e-04;
*/

using namespace std;


class FreeStrm
{
protected:
	double ****densityTable, ****shiftedTable;
      double ****densityTableInterp;
	double ***dNd2rdyTable;
	double ***dEd2rdyTable;
	double Xmax,Ymax,Xmin,Ymin,dx,dy, PTmax, PTmin, dpt;
	int    nRap;
	double rapMin, rapMax;
	int    Maxx, Maxy, MaxPT;
      int    NpTinterp;
	double Taui, Tauf, Phip;
	double Xcm, Ycm;    //center of the energy density profile

	void findCM(double ***source, const int iRap = 0);  //find the center of the profile

	//A general function to shift data table to a new center (x0, y0) and/or rotate phi clockwisely
	double getShiftedProfile(double ***data, int i, int j, double x0, double y0, 
		double phi = 0, bool limit=false, const int iRap=0);  


public:
	FreeStrm(double xmax, double ymax, double ptmin, double ptmax, double dx0,double dy0,
			    double dpt, int nrap, double rmin, double rmax, double taui, double tauf);
	~FreeStrm();

	void ReadTable(string filename);
	void GetDensity(int iy, int i, int j, double phip, double* results);  //debugging
      void GetDensityInterp(int iRap, int i, int j, double phip, double* results);
      void ShiftDensityInterp(const int iRap, double phip, double ****shiftedTableInterp);
	void ShiftDensity(const int iy, double phip);   //free stream density to another coordinate with a specific angle
      void InterpDensity(const int iRap, double* pT, int npT);
	double GetShiftdeDensity(int iy, int i, int j, int ipt) {return shiftedTable[iy][i][j][ipt];};

	void OutputTable(const char* filename, const int iy);  //output free steamed and pt integrated gluon density
	void dumpBlockTable(const char *filename, double ***data, const int iy);

	void CreateDataTable(const char *filename, const int iy=0); //testing
	double GaussProfile(const int iy, int i, int j, int ipt);  		//generating test profile
	double BoxProfile(const int iRap, int i, int j, int ipt);
	void GetInteDensity(const int iy);   //test integration

	void generateEdTable(const int iRap=0);
	double getEpx(int n, const int iRap=0);    //calculate eccentricity in the free-streaming stage


};


#endif
