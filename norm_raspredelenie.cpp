/*
 * mat.cpp
 *
 *  Created on: 15 мар. 2019 г.
 *      Author: andrey
 */
#include <statistics.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>
using namespace std;
using namespace alglib;
vector<double> shapir_uilk_txt(int number, string name)
{
	fstream date;
	date.open("date_shapir_uilk.txt",ios::out);
	if(!date)
			cout << "NO FILE OPEN"<<endl;
		fstream fopen;
		fopen.open(name);
		if(!fopen)
			cout << "NO FILE OPEN"<<endl;
		vector<string> str; // замена запятых на точки
		while(!fopen.eof())
		{
			string s;
			getline(fopen,s);
			int i = -1;
			while(++i < s.size())
				if(s[i] == ',')
					s[i] = '.';
			str.push_back(s);
		}
		fopen.close();
		date << str[number] <<endl;
		date.close();
		fstream open_date;
		open_date.open("date_shapir_uilk.txt");
		double d = 0;
		vector<double> result;
		while(open_date >> d)
			result.push_back(d);
		open_date.close();
		return result;
}

void norm_raspredelenie(const real_1d_array& x, fstream& fout)
{
	vector<double> result;
	real_1d_array x_all;
	for(int i = 0; i < x.length(); i++)
		result.push_back(x[i]);
	sort(result.begin(), result.end());
	double mat_o = samplemean(x);
	double disp = samplevariance(x);
	int i = 0;
	fout << "Общая выборка после сортировки по возрастанию:"<< endl;
	for(auto s : result)
	{
		
		fout << s << " ";
		if( i % 10 == 0)
			fout << endl;
		i ++;
	}
	fout << endl;
	fout <<	"Дисперсия: " << disp<<endl;
	fout << "Мат. ожидание: " << mat_o<<endl;
	vector<double> koff_a = shapir_uilk_txt((result.size() % 2 == 0 ? result.size() / 2 : (result.size() - 1) / 2),  "shapir-uilk-koff.txt" );
	vector<double> wer = shapir_uilk_txt((result.size() % 2 == 0 ? result.size() / 2 : (result.size() - 1) / 2), "shapir-uilk-wer.txt" );
	double w = 0;
	if(result.size() <= 50)
	{
		fout << "Критерий Шапира-уилка: "<< endl;
		
		
		for(int i  = 0; i <= koff_a.size(); i++)
		{
			double p1 = koff_a[koff_a.size() - i - 1] / 10000;
			double p2 = result[i] - result[result.size() - 1 - i];
			w +=  p1 * p2 ;
		}
		//cout << w << " "; 
		w = pow(w,2);
		//cout << disp << endl;
		w /= disp;
		fout <<"статистика W : "<<  w <<endl;
		vector<double> a = {0.9,0.95,0.99};
		int i = 0;
		for(auto s : wer)
		{
			fout << " W(a): "<< s<<endl;
			if (w > s)
				fout << "Выборка починяется нормальному закону распределения с долей вероятности равной " << a[i] << endl;
			else
				fout << "Выборка НЕ починяется нормальному закону распределения с долей вероятности равной " << a[i] << endl;
			i++;
		}
	}	
	else 
	{
		int n = result.size();
		fout<<"Так как значений в выборке больше чем 50, применяем критерий Шапира-Френчиа" << endl;
		double down_c;
		for(int i = 0 ; i < n / 2; i++)
		{
			double p = ((i + 1)-(3/8)) / (n + (1/4));
			down_c += pow(4.91 * (pow(p, 0.14)  - pow(1 - p, 0.14)),2);
		}
		
		down_c = pow(down_c, 0.5);
		double up_c;
		for(int i  = 1; i <= n / 2; i++)
		{	
			double p = ((i + 1)-(3/8)) / (n + (1/4));	
			up_c = 4.91 * (pow(p, 0.14)  - pow(1 - p, 0.14));
			double up = (up_c / down_c);
			double down = result[i - 1] - result[result.size() - i];
			w +=up * down;
		}
		w /= disp;
		fout <<"статистика W : "<<  w <<endl;
		double w1 = (-0.0084*pow(result.size(), 4) + 1.2513 * pow(result.size(),3) - 70.724 * pow(result.size(), 2) + 1890 * result.size() + 73840)/ pow(10,5);
		double w2 = (-0.0113*pow(result.size(), 4) + 1.656 * pow(result.size(),3) - 91.88 * pow(result.size(), 2) + 2408.6 * result.size() + 67608)/pow(10,5);
		double w3 = (-0.0148*pow(result.size(), 4) + 2.1875 * pow(result.size(),3) - 122.61 * pow(result.size(), 2) +  3257.73* result.size() + 55585)/pow(10,5);
		if (w > w1)
				fout<<"W'0.1: "<<w1 <<endl  << "Выборка починяется нормальному закону распределения с долей вероятности равной " << 0.1 << endl;
		else 
			fout<<"W'0.1: "<<w1 <<endl  << "Выборка НЕ починяется нормальному закону распределения с долей вероятности равной "<< 0.1 << endl;
		if(w > w2 )
				fout <<"W'0.05: " << w2 <<endl <<"Выборка починяется нормальному закону распределения с долей вероятности равной " << 0.05 << endl;
		else 
			fout<<"W'0.05: " << w2 <<endl << "Выборка НЕ починяется нормальному закону распределения с долей вероятности равной "<< 0.05 << endl;
		if(w >w3 )
			fout <<"W'0.01: "<<w3 <<endl<<"Выборка починяется нормальному закону распределения с долей вероятности равной " << 0.01 << endl;
		else 
			fout<<"W'0.01: "<<w3 <<endl << "Выборка НЕ починяется нормальному закону распределения с долей вероятности равной "<< 0.01 << endl;	
	}

}