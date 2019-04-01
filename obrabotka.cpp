/*
 * file.cpp
 *
 *  Created on: 15 мар. 2019 г.
 *      Author: andrey
 */

#include </home/andrey/Документы/mat.stat/Mat.stat/includes/statistics.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
# define ABS(x)  ( (x < 0) ? -(x) : x )
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
	bool flag = true;
	fout<< "Процедура отбраковки аномальных измерений :"<< endl;
	double mat;
	double disp;
	double assmetri;
	double eksess;
	while(flag)
	{
		flag = false;
		samplemoments(x,mat,disp,assmetri,eksess);
		double disp_sm = disp *(x.length() - 1)/(x.length());
		fout <<"Количество Элементов: n = "<< x.length()<< endl;
		fout <<"Математическое ожидание: M = "<< mat<< endl;
		fout <<"Смешеное оценка дисперсии: Si = "<< disp_sm<< endl;
		fout <<"Не смешеное оценка дисперсии: Si = "<< disp << endl;
		fout <<"Третий момент распределение: m3 = "<< assmetri * pow(disp,3/2)<< endl;
		fout <<"Четвертый момент распределение: m4 = "<< eksess * pow(disp,2)<< endl;
		fout <<"Коэффициент ассиметрии : В1 = "<< assmetri<< endl;
		fout <<"Коэффициент эксцесса: В2 = "<< eksess + 3<< endl;

		for(int i = 0; i < x.length(); i++)
		{
			if(invstudenttdistribution(x.length() - 2, 0.95)  <  ABS((x[i] - mat))/pow(disp,0.5))
				{
					fout<< "Аномальное измерение: " << x[i] << endl;
					flag = true;
					x = delet(x,i);
				}
		}
		fout << endl;
	}
	fout << " Проверка на однородность :"<< endl;
	for(int i = 0; i < a.size(); i++)
	{
		if(invnormaldistribution(1- ((1 - a[i])/ 2)) >= (ABS(sampleskewness(x))/pow((6/x.length()),1/2)))
				fout << "результат серии не противоречит нормальному закоу распределения с уровнем доверия :" << a[i] << endl;
	}
	return(x);
}
