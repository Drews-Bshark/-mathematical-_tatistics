#include </home/andrey/Документы/mat.stat/Mat.stat/-mathematical-_tatistics/includes/statistics.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
using namespace std;
using namespace alglib;
// талерантный интервал
 int main()
 {
	double a = 0.95;
	double InfersNormalDistruution_k_z = invnormaldistribution(a);// квантиль распеределение заданого уровня доверия 
	double InfersNormalDistruution_y = invnormaldistribution(a);
	double disp  = 0.249;
	double mat = 1.12;
	int n = 14;
	cout << InfersNormalDistruution_k_z<<endl;
	double up = 1  - pow( InfersNormalDistruution_y, 2)/ (2* (n - 1));
	double up1 =  pow(InfersNormalDistruution_k_z, 2) / (2 * (n - 1));
	double down = 1 - pow(InfersNormalDistruution_y,2)/(2 * (n - 1));
	double k = (InfersNormalDistruution_k_z + InfersNormalDistruution_y - pow((up/n) + up1,0.5))/down;
	//cout << pow(up/n + up1,0.5)<<endl;
	//cout << down<< endl;
	cout << k<<endl;
	double x_min = mat - k * disp;
	if(x_min < 0)
		x_min  = 0;
	double x_max = mat  + k * disp;
	cout << x_max<<endl;
	cout << x_min<< endl;
	//x = xs +- k*S // выбор спреднего 
	return 0;
 }