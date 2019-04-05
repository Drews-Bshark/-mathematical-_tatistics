/*
 * stat_touchnost.cpp
 *
 *  Created on: 3 апр. 2019 г.
 *      Author: andrey
 */

#include </home/andrey/Документы/mat.stat/Mat.stat/-mathematical-_tatistics/includes/statistics.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
using namespace alglib;
using namespace std;
void stat_touchnost(const double& a,const double& k, fstream& fout)
{
	fout << "Ошибки первого и второго рода равна: "<<  1 - a <<endl;
	for(int i = 1; i <=10; i++)
	{
		
		fout << setw(10);
		fout << i * 10 ;
		if(i*10 < k && k < ((i + 1) * 10 ))
		{
				fout << setw(10);
				fout << k ;
		}
	}
	fout << endl;
	//cout << "Ошибка"<< endl;
	//cout<< 	"мат.ожидания";
	for(int i = 1; i <=10; i++)
	{
		if(i * 10 == k)
			fout << "M: ";
		fout << setw(10);
		fout <<  (invnormaldistribution(a) - invnormaldistribution(1 - a))/pow((i*10)/2,0.5);
		if(i*10 < k && k <((i + 1) * 10) )
		{
			fout << setw(10);
			fout<<"M: " << (invnormaldistribution(a) - invnormaldistribution(1 - a))/pow((k)/2,0.5);
		}
	}
	fout << endl;
	//cout << "ошибка"<<endl;
	//cout<<"дисперсии";
	for(int i = 1; i <=10; i++)
	{
	
		if(i * 10 == k)
			fout << "K: ";
		fout << setw(10);
		fout << invfdistribution(i *10 - 1,i * 10 - 1, 1 - a) / invfdistribution(i *10 - 1,i * 10 - 1,a);

		if(i*10 < k && k <((i * 10 + 10)))
		{
			fout << setw(10);
			fout<<"K: " <<invfdistribution(k - 1,k - 1, 1 - a) / invfdistribution(k - 1,k - 1,a);
		}
		}
	fout << endl;

}

