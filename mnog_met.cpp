#include </home/andrey/Документы/mat.stat/Mat.stat/includes/statistics.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

using namespace std;
using namespace alglib;
void mnog_met(map <int,real_1d_array>& x, double a, fstream& fout,fstream& problam)
{
	int seri = x.size();
	fout << " ПРОВЕРКА СТАТИСТИЧЕСКОЙ ОДНОРОДНОСТИ РЕЗУЛЬТАТОВ ИСПЫТАНИЯ МНОГОМЕРНЫМ МЕТОДОМ "<< endl<<endl;
	fout << "ПОВЕРКА ГИПОТЕЗЫ ПО КРИТЕРИЮ БАРТЛЕТТА"<<endl;
	double s = 0;
	double f = 0;
	for(int i = 1; i <= seri; i++)
		f += (x[i].length() - 1);
	fout<< "		f:" << f<< endl;
	for(int i = 1; i <= seri; i++)
		s += (x[i].length() - 1)*samplevariance(x[i]);
	s = s/f;
	fout << "		S:" <<s<<endl;
	double c = 0;
	double sum_fi = 0;
	double v = 0;
	for(int i = 1; i <= seri; i++)
		sum_fi += (1.0 / (x[i].length())); 
	double ln_s = 0;
	for(int i = 1; i <= seri; i++)
		ln_s  += (x[i].length() - 1) * log(samplevariance(x[i]));
	fout <<"LN_S " <<ln_s<<endl;
	fout << "f*ln(Si) "<< f*log(s)<<endl;
	v = f*log(s)- ln_s;
	c = 1 + ((sum_fi - 1.0/f)/ (3*(seri - 1)));;
	fout << " 		C: "<<c << endl;
	fout << " 		V: "<< v<<endl;
	fout << " 		T: "<< v/c<< endl;
	fout << "Квантильнормльного распределения равен "<< invchisquaredistribution(seri - 1,1 - a)<<endl;
	if(v/c <= invchisquaredistribution(seri - 1,1 - a))
	{
		fout << "Так как статистический критерий Бартлетта меньше "<<" Квантиля хи-квадратного распределения со степенями свободы ";
		fout<<seri - 1<<" То можно принять гипотезу о равенстве всех "<<seri<< " дисперсий"<< endl;
		fout <<endl;
		fout<<"Проверка равенства мат ожиданий многомерном методом"<< endl;
		double sk = 0;
		double xi = 0;
		double sum_n = 0;
		for(int i = 1; i <=seri; i++)
			sum_n += x[i].length() - 1;
		for(int i = 1; i <= seri; i++)
			xi += x[i].length() * samplemean(x[i]); 
		xi /=  (f + seri);
		for(int  i = 1; i <= seri; i++)
		{		
			sk += x[i].length() * pow(samplemean(x[i]) - xi,2);
		}
		sk  /= (seri - 1);
		fout << " Sk: "<<sk << endl;
		fout <<"Квантиль распределение ФИШЕРА при степеням свободы: "<<endl<<"		"<<seri - 1 <<endl<<"		"<< f<<endl<<"Равен:		"<<invfdistribution(seri - 2,f,1 - a)<<endl;
		if(sk/s <=invfdistribution(seri - 2,f,1 - a))
			fout << "На уровне ожидания "<<a<< " по критерию Фишера принять гипотезу о равенстве мат ожиданий МОЖНО!!"<< endl;
		else
			problam <<  " На уровне ожидания "<<a<< " по критерию Фишера принять гипотезу о равенстве мат ожиданий НЕЛЬЗЯ!!"<< endl;
			fout <<  " На уровне ожидания "<<a<< " по критерию Фишера принять гипотезу о равенстве мат ожиданий НЕЛЬЗЯ!!"<< endl;

	}
	else
	{
		fout <<"Так как статистический критерий Бартлетта больше ";
		fout<<" Квантиля хи-квадратного распределения со степенями свободы "<< seri - 1;
		problam<<" То НЕЛЬЗЯ принять гипотезу о равенстве всех "<<seri<< " дисперсий (Критерий Батлера)"<< endl;
		fout<<" То НЕЛЬЗЯ принять гипотезу о равенстве всех "<<seri<< " дисперсий (Критерий Батлера)"<< endl;
	}
	fout << endl << endl;
	
}
