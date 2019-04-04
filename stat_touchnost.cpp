/*
 * stat_touchnost.cpp
 *
 *  Created on: 3 апр. 2019 г.
 *      Author: andrey
 */

#include </home/andrey/Документы/mat.stat/Mat.stat/-mathematical-_tatistics/includes/statistics.h>
#include <iostream>
#include <vector>
#include <iomanip>
using namespace alglib;
using namespace std;
int main()
{
	double a = 0.9;
	cout << "РАССЧИТЫВАЕМ ТОЧНОСТЬ СТАТИСТИЧЕСКИХ РЕШЕНИЙ ДЛЯ ОШТБОК 1 И 2 РОДА И КОЛИЧЕСТВО ИЗМЕРЕНИЙ: "<< endl;
	cout << " Постоить график зависимости точности статисческих решений от объема исследованной выборки" << endl;
	//cout << "N";
	for(int i = 1; i <=10; i++)
	{
		cout << setw(10);
		cout << i * 10 ;
	}
	cout << endl;
	//cout << "Ошибка"<< endl;
	//cout<< 	"мат.ожидания";
	for(int i = 1; i <=10; i++)
	{
		cout << setw(10);
		cout <<  (invnormaldistribution(a) - invnormaldistribution(1 - a))/pow((i*10)/2,0.5);
	}
	cout << endl;
	//cout << "ошибка"<<endl;
	//cout<<"дисперсии";
	for(int i = 1; i <=10; i++)
	{
		cout << setw(10);
		cout << invfdistribution(a,1,2);
	}
	cout << endl;
	return 0;

}

