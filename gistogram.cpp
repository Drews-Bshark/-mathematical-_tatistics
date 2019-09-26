#include <statistics.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <iomanip>
#include <algorithm>
#define MAX_SERIA 3
using namespace std;
using namespace alglib;
void gistogram (const map <int,real_1d_array>& x,const real_1d_array& y, fstream& fout)
{
	
	map<int,vector<double>> result;
	
	for(int i = 1; i <= MAX_SERIA; i++)
	{
		for(int j = 0; j < x.at(i).length(); j++)
		{
			result[i].push_back(x.at(i)[j]);
		}
	}
	for(int i = 1; i <= MAX_SERIA; i++)
	{
		fout << "Для серии "<< i <<endl;
		map<double, int> count;
		for(auto c : result[i])
			++count[c];
		fout <<endl;
		for(auto s : count)
			fout<<setw(3) << s.second <<" ";
		fout <<endl;
		count.clear();
	}
	map<double, int> count;
	fout << "Для общей выборки "<<endl;
	for(int i = 0; i < y.length(); i++)
	{
		++count[y[i]];
	}
	fout <<endl;
	for(auto s : count)
		fout<<setw(3) << s.second <<" ";
	fout <<endl;
	fout << " Значения для общей выборки"<< endl;
	double mat;
	double disp;
	double assmetri;
	double eksess;
	samplemoments(y,mat,disp,assmetri,eksess);
	fout <<"Математическое ожидание: M = "<< mat<< endl;
	fout <<"Не смешеное оценка дисперсии: Si = "<< disp << endl;
	fout <<"Коэффициент ассиметрии : В1 = "<< assmetri<< endl;
	fout <<"Коэффициент эксцесса: В2 = "<< eksess + 3<< endl;
	/*eksess+= 3;
	double l;
	double q; 
	l = 0.313 * pow(2* eksess/ (3 - eksess), 0.5);
	q = (5*eksess)/ (double)(2*(3 - eksess));
	fout << "L: "<<l <<endl;
	fout << "q:"<<q <<endl;
	double n;
	double p;
	n =  1 - (pow(0.0039, 2) /  pow(l, 2));
	n = pow(n,q);
	n*=y.length();
	fout << "n: "<<n<<endl 
	double down_p = 1 - (pow(,00, 2)/pow(l, 2));
	down_p = pow(down_p,q);
	double up = 0;*/
	
}