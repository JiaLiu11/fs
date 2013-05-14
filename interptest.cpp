#include <iomanip>
#include <string>
#include <cstring>
#include "mistools.h"
#include <iostream>

using namespace std;

int main()
{
	double v1,v2,v3,v4,result;
	double dx=0.5,dy=0.5;
	double x1,y1;

	v1=12.0, v2=19.0, v3=23.5, v4=16.5;
	x1=0.5, y1=0.5;

	result=Bilinear2dInterp(x1,y1,dx,dy,v1,v2,v3,v4);
	cout<<"result="<<result<<endl;

	double result2=Linear1dInterp(1.6,1,2,3);
	cout<<"result for 1d is:"<<result2<<endl;

	return 0;

}