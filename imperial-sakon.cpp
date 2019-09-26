#include <statistics.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>
#define ABS(x)  ( (x < 0) ? - (x) : x )
using namespace std;
using namespace alglib;
void imperial_sakon(const real_1d_array& x, fstream& fout,const map<double, double>& lamda)
{	
	vector<double> result;
	for(int i = 0; i < x.length(); i++)
		result.push_back(x[i]);
	sort(result.begin(), result.end());
	map<double, int> count;
	for(auto s : result)
		++count[s];
	for(auto s : count)
		fout<<setw(10)<< s.first ;
	fout <<endl;
	for(auto s : count)
		fout<<setw(10)<< s.second;
	fout << endl;
	int sum = 0;
	for(auto s : count)
	{
		sum += s.second;
		fout <<  setw(10)<< sum / (double)(result.size() + 1);
	}
	sum = 0;
	fout<<endl << "Значения нормального рспределения"<< endl;
	for(auto s : count)
		fout << setw(10)<< normaldistribution(s.first);
	fout << endl;
	fout << "Проверка согласованности эмперического и нормального распредеения(Критерий Колмагорова): "<< endl<<endl;
	
	double d = 0; 
	for(auto s : count)
	{
		sum += s.second;
		if(d < ABS(sum / (double)(result.size() + 1) - normaldistribution(s.first)))
			d =  ABS(sum / (double)(result.size() + 1) - normaldistribution(s.first));
	}
	double lam = d * pow(result.size(),0.5);
	fout << "Лямда равна: " << lam << endl;
	if(lam <= 2)
	{
		for(auto s : lamda)
			if(s.first > lam)
			{
				fout << "P(l): "<<"Выыод: имперический и теоретический законы равны c долей вероятности: " << s.second << endl;
				break;
			}
	}
	else 
		fout << "P(l): -> 0"<<endl<< "Выыод: имперический и теоретический законы не равны"<<endl;
}