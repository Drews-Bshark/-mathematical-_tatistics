/*
 * file.cpp
 *
 *  Created on: 15 мар. 2019 г.
 *      Author: andrey
 */

#include <statistics.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
# define ABS(x)  ( (x < 0) ? - (x) : x )
using namespace std;
using namespace alglib;
real_1d_array delet (const real_1d_array& x,int i)
{
	int j = 0;
	real_1d_array y;
	y.setlength(x.length() - 1);
	while(j < x.length() - 1)
	{
		if(j >= i)
			y[j] = x[j + 1];
		else
			y[j] = x[j];
		j++;
	}
	return y;

}
real_1d_array static_funk(const vector<double>& number,const vector <double>& a,fstream& fout)
{
	
	real_1d_array x;
	x.setcontent(number.size(),&(number[0]));
	fout<< "Процедура отбраковки аномальных измерений :"<< endl;
	double mat;
	double disp;
	double assmetri;
	double eksess;
	int n = x.length();

	samplemoments(x,mat,disp,assmetri,eksess);
	eksess += 3;
	double disp_sm = disp *(x.length() - 1)/(x.length());
	double m3 = 0;
	double m4 = 0;
	for (int i = 0; i < x.length(); i++)
		m4 += pow(x[i] - mat,4);
	m4 /= x.length();
	for(int i  = 0 ; i < x.length(); i++)
		m3 +=pow(x[i] - mat,3);
	m3 /= x.length();
	fout <<"Количество Элементов: n = "<< x.length()<< endl;
	fout <<"Математическое ожидание: M = "<< mat<< endl;
	fout <<"Смешеное оценка дисперсии: Si = "<< disp_sm<< endl;
	fout <<"Не смешеное оценка дисперсии: Si = "<< disp << endl;
	fout <<"Третий момент распределение: m3 = "<< m3<< endl;
	fout <<"Четвертый момент распределение: m4 = "<< m4<< endl;
	fout <<"Коэффициент ассиметрии : В1 = "<< assmetri<< endl;
	fout <<"Коэффициент эксцесса: В2 = "<< eksess<< endl;
	bool flag = false;
	//cout << invstudenttdistribution(x.length() - 2, 0.95)<<endl;
	for(int i = 0; i < x.length(); i++)
	{
		//cout<<"test: "<< x[i]<< " " << ABS((x[i] - mat))/pow(disp,0.5)<< " "<<pow(disp_sm,0.5)<<endl;
		if(invstudenttdistribution(n - 2, 0.95)  <  ABS((x[i] - mat))/pow(disp,0.5))
			{
			fout<< "Аномальное измерение: " << x[i] << endl;		
			x = delet(x,i);
			flag = true;
			}
	}
	fout << endl;
	if(flag)
	{
		samplemoments(x,mat,disp,assmetri,eksess);
		m4 = 0;
		m3 = 0;
		for (int i = 0; i < x.length(); i++)
			m4 += pow(x[i] - mat,4);
		m4 /= x.length();
		for(int i  = 0 ; i < x.length(); i++)
			m3 +=pow(x[i] - mat,3);
		m3 /= x.length();
		eksess += 3;
		disp_sm = disp *(x.length() - 1)/(x.length());	
		fout <<"Количество Элементов: n = "<< x.length()<< endl;
		fout <<"Математическое ожидание: M = "<< mat<< endl;
		fout <<"Смешеное оценка дисперсии: Si = "<< disp_sm<< endl;
		fout <<"Не смешеное оценка дисперсии: Si = "<< disp << endl;
		fout <<"Третий момент распределение: m3 = "<<  m3<< endl;
		fout <<"Четвертый момент распределение: m4 = "<< m4<< endl;
		fout <<"Коэффициент ассиметрии : В1 = "<< assmetri<< endl;
		fout <<"Коэффициент эксцесса: В2 = "<< eksess<< endl;
	}
	fout << " Проверка на однородность :"<< endl;
	for(int i = 0; i < a.size(); i++)
	{
		if(invnormaldistribution(1- ((1 - a[i])/ 2)) >= (ABS(sampleskewness(x))/pow((6/x.length()),1/2)))
				fout << "результат серии не противоречит нормальному закоу распределения с уровнем доверия :" << a[i] << endl;
	}
	return(x);
}
