#include <statistics.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>
using namespace std;
using namespace alglib;
const double R = 0.95;
const double Y = 0.9;
const vector<double> dop_x = {1.5,1.6,1.7};
void mat_interval(const real_1d_array& x, fstream& fout) 
{
	int n = x.length();
	double inv_norm_disp = invnormaldistribution((Y +1)/2);
	double mat_o = samplemean(x);
	double disp = samplevariance(x);
	fout << "Расчет оценки " << endl;
	fout << mat_o << endl;
	fout << disp << endl;
	vector <double> k;
	for(auto s : dop_x)
	{
		double d = (s - mat_o) / pow(disp,0.5);
		k.push_back(d);
	}
	fout << "Вычисление доверительного интервала " <<endl;
	for(auto s : k)
	{
		fout <<"K: "<< s << endl;
		double sl = pow(s,2)/(double)(2*(n - 1)) + 1/n;
		fout  << "Ur" << ": ";
		fout <<setw(10)<< s -inv_norm_disp * pow(sl,0.5)<< " ";
		fout <<setw(10)<< s + inv_norm_disp * pow(sl,0.5)<<endl;
		fout  << "R" << ": ";
		fout <<setw(10)<<normaldistribution( s - inv_norm_disp * pow(sl,0.5))<<" ";
		fout <<setw(10)<<normaldistribution( s + inv_norm_disp * pow(sl,0.5))<<endl;
	}
	fout <<"Построить области принятия гипотез "<<endl;
	for(auto s : k)
	{
		fout << "K: "<< s << endl;
		double sl = pow(s,2)/(double)(2*(n - 1)) + 1/n;
		fout  << "Ur" << ": ";
		fout <<setw(10)<<normaldistribution(R)  -inv_norm_disp * pow(sl,0.5)<< " ";
		fout <<setw(10)<<normaldistribution(R) + inv_norm_disp * pow(sl,0.5)<<endl;
		fout  << "R" << ": ";
		fout <<setw(10)<<normaldistribution( normaldistribution(R) - inv_norm_disp * pow(sl,0.5))<<" ";
		fout <<setw(10)<<normaldistribution( normaldistribution(R) + inv_norm_disp * pow(sl,0.5))<<endl;
	}
	fout << " R > R3" << endl;
	for(auto s : k)
	{
		fout << "K: "<< s << endl;
		double sl = pow(s,2)/(double)(2*(n - 1)) + 1/n;
		fout  << "R"  << ": ";
		double U_dop = normaldistribution( normaldistribution(R) + inv_norm_disp * pow(sl,0.5));
		fout <<setw(10)<<normaldistribution( normaldistribution(R) + inv_norm_disp * pow(sl,0.5))<<" ";
		fout <<setw(10)<<normaldistribution( U_dop + inv_norm_disp * pow(sl,0.5))<<endl;
		fout  << "X" << ": ";
		fout <<setw(10)<<normaldistribution(R) + inv_norm_disp * pow(sl,0.5)<<" ";
		fout <<setw(10)<<U_dop + inv_norm_disp * pow(sl,0.5)<<endl;
	}
	fout << " R < R3" << endl;
	for(auto s : k)
	{
		fout << "K: "<< s << endl;
		double sl = pow(s,2)/(double)(2*(n - 1)) + 1/n;
		fout  << "R"  << ": ";
		fout <<setw(10)<<normaldistribution( 0)<<" ";
		fout <<setw(10)<<normaldistribution( normaldistribution(R) - inv_norm_disp * pow(sl,0.5))<<endl;
		fout  << "X" << ": ";
		fout <<setw(10)<<0<<" ";
		fout <<setw(10)<<normaldistribution(R) + inv_norm_disp * pow(sl,0.5)<<endl;		
	}
	fout << " R = R3" << endl;
	for(auto s : k)
	{
		fout << "K: "<< s << endl;
		double sl = pow(s,2)/(double)(2*(n - 1)) + 1/n;
		fout  << "R"  << ": ";
		double U_dop = normaldistribution( normaldistribution(R) + inv_norm_disp * pow(sl,0.5));
		fout <<setw(10)<<normaldistribution( U_dop - inv_norm_disp * pow(sl,0.5))<<" ";
		fout <<setw(10)<<U_dop<<endl;
		fout  << "X" << ": ";
		fout <<setw(10)<<U_dop - inv_norm_disp * pow(sl,0.5)<<" ";
		fout <<setw(10)<< normaldistribution(R) + inv_norm_disp * pow(sl,0.5)<<endl;	
	}
	fout << " R != R3" << endl;
	for(auto s : k)
	{
		fout << "K: "<< s << endl;
		double sl = pow(s,2)/(double)(2*(n - 1)) + (1/n);
		fout  << "R"  << ": ";
		double U_dop = normaldistribution( normaldistribution(R) + inv_norm_disp * pow(sl,0.5));
		fout <<setw(10)<<normaldistribution( normaldistribution(R) - inv_norm_disp * pow(sl,0.5))<<" ";
		fout <<setw(10)<<normaldistribution( U_dop - inv_norm_disp * pow(sl,0.5));
		fout  << " U ";
		fout <<setw(10)<<normaldistribution( normaldistribution(R) + inv_norm_disp * pow(sl,0.5))<<" ";
		fout <<setw(10)<<normaldistribution( U_dop + inv_norm_disp * pow(sl,0.5))<<endl;
		
		fout  << "X" << ": ";
		fout <<setw(10)<< normaldistribution(R) - inv_norm_disp * pow(sl,0.5)<<" ";
		fout <<setw(10)<< U_dop - inv_norm_disp * pow(sl,0.5);
		fout  << " U ";
		fout << setw(10)<< normaldistribution(R) + inv_norm_disp * pow(sl,0.5)<<" ";
		fout << setw(10)<< U_dop + inv_norm_disp * pow(sl,0.5)<<endl;
	}
	fout << "Расчет точности статистических решений" <<endl;
	for(auto s : k)
	{
		fout << "K: "<<s<<endl;
		for(double y = 0.9; y < 0.96; y +=0.05)
		{
			fout << "Y: "<<y << endl;
			fout << "n:  ";
			for(int i = 10; i <=100; i += 10)
			{
				fout << setw(10) << i;
		
			}
			fout <<endl;
			fout << "U:  ";
			for(int i = 10; i <=100; i += 10)
			{
				double sl = (pow(s,2) /(double)(2 * (i -1))) + (1/i);
				double sl1 = normaldistribution(y) + normaldistribution(1 - y);
				fout << setw(10) << sl1 * pow(sl,0.5);		
			}
			fout << endl;			
		}
	}

}