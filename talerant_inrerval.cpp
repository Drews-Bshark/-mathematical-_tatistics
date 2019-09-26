#include <statistics.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>
using namespace std;
using namespace alglib;
// талерантный интервал
 void talerant_interval(double a, int n, double disp, double mat, fstream& fout)
 {
	 fout << "Мат ожидание " << mat << endl;
	 fout << "Дисперсия "<< disp << endl;
	 
	double InfersNormalDistruution_k_z = invnormaldistribution(a);// квантиль распеределение заданого уровня доверия 
	double InfersNormalDistruution_y = invnormaldistribution(a);
	double up = 1  - pow( InfersNormalDistruution_y, 2)/ (2 * (n - 1));
	double up1 =  pow(InfersNormalDistruution_k_z, 2) / (2 * (n - 1));
	double down = 1 - pow(InfersNormalDistruution_y,2)/(2 * (n - 1));
	double k = (InfersNormalDistruution_k_z + InfersNormalDistruution_y - pow((up/n) + up1,0.5))/down;
	
	fout <<"Для нормального распределения K: " << k << endl;
	double x_min = mat - k * disp;
	if(x_min < 0)
		x_min  = 0;
	double x_max = mat  + k * disp;
	fout << "Максимальный X:"<< x_max << endl;
	fout <<"Минимальный Х:"<<x_min << endl<<endl<<endl<<endl;
	vector<double> b2 = {3, 1.8, 4.2, 3.869, 2.1};
	for(int i = 0; i < b2.size(); i++)
	{
		//cout << endl;
		double down_n = (b2[i] - 3) / 2; 
		//cout <<down_n << endl;
		double down_n_2 = (n - 1) /(double) n;
		//cout << down_n_2 << endl;
		double n_ekv = (n - 1) /  (1 + down_n * down_n_2);
		//cout << n_ekv<<endl;
		//cout << invchisquaredistribution(n - 1,a)<<endl;
		up = 1  - pow( InfersNormalDistruution_y, 2)/ (2* (n_ekv));
		//cout << "UP  = " << up << endl;
		up1 =  pow(invchisquaredistribution(1,1 - a), 2) / (2 * (n_ekv ));
		//cout << "UP1  = " << up1 << endl;
		down = 1 - pow(InfersNormalDistruution_y ,2)/(2 * (n_ekv ));
		//cout << "down  = " << down << endl;
		k = (invchisquaredistribution(1,1 - a) + InfersNormalDistruution_y - pow((up/n) + up1,0.5))/down;
		//cout << "K  = " << k << endl;
		fout << "Для распределений пирсона ";
		if(i == 0)
			fout << "Нoрмальное "<<endl<< "K: "<< k<<endl;
		if(i == 1)
			fout << "Равномерное "<<endl<<"K: "<< k<<endl;
		if(i == 2)
			fout << "Логическое  "<<endl<< "K: "<< k<<endl;
		if(i == 3)
			fout << "Полунормальное  "<<endl<< "K: "<< k<<endl;
		if(i == 4)
			fout << "Рaпределение Пирсона (эксперементальный) "<<endl<< "K: "<< k<<endl;	
		fout << "Максимальный X:"<< mat + k * disp << endl;
		if(mat - k * disp > 0 )
			fout <<"Минимальный Х:"<<mat - k * disp << endl; 
		else
			fout <<"Минимальный Х:"<< 0 << endl; 
	}
 }