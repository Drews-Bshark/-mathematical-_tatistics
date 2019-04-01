/*
 * mat.cpp
 *
 *  Created on: 15 мар. 2019 г.
 *      Author: andrey
 */
#include </home/andrey/Документы/mat.stat/Mat.stat/includes/statistics.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#define MAX_STUDENT 30
#define MAX_SERIA 3
using namespace std;
using namespace alglib;
real_1d_array static_funk(const vector<double>& number, const vector <double>& a,fstream& fout);
string data_txt(int student_number);
void mnog_met(map <int,real_1d_array>& x, double a, fstream& fout, fstream& problam);
const int N[] = {16,27,28};
const vector<double> a = {0.9,0.95,0.99};
int main()
{
	fstream problam;
	string n;
	problam.open("problam.txt", ios::out);
	fstream fout;
	fout.open("reshenie.txt", ios ::out);
	for(int k = 1; k < 30; k++)
	{
		n = data_txt(k);
		problam <<endl<< n <<endl;
		map <int,real_1d_array> x;
		fstream file;
		file.open("data.txt");
		if(!file)
			cout << "No file open"<< endl;
		map<int,vector<double>> seria;
		int j = 0;
		while(++j <= MAX_SERIA)
		{
			vector <double> d_v;
			int i = -1;
			while(++i < N[j - 1])
			{
				double d;
				file >> d;
				d_v.push_back(d);
			}
			seria[j] = d_v;
		}
		file.close();
		fout<<endl<<endl;
		fout <<"Ф.И.О "<< n <<endl<<endl<<endl;
		fout<<"1. ПРОЦЕДУРА ОТБРАКОВКИ АНОМАЛЬНЫХ ИЗМЕРЕНИЙ, ОПРЕДЕЛЕНИЕ ОДНОРОДНОСТИ РЕЗУЛЬТАТОВ ИЗМЕРЕНИЙ :"<< endl;
		for(int i = 1; i <= MAX_SERIA; i++)
			x[i] = (static_funk(seria[i],a,fout));
		double disp_min =samplevariance(x[1]);
		ae_int_t count_min = 1;
		double disp_max = samplevariance(x[1]);
		ae_int_t count_max = 1;
		int i = 1;
		while(i <= x.size())
		{
			if(samplevariance(x[i]) < disp_max)
			{
				disp_max = samplevariance(x[i]);
				count_max = i;
			}
			if(samplevariance(x[i]) > disp_min)
			{
				disp_min = samplevariance(x[i]);
				count_min = i;
			}
			i++;
		}

		fout<< endl<< endl;
		fout<<"2 ПРОВЕРКА ОДНОРОДНОСТИ РЕЗУЛЬТАТОВ ЭКСПЕРЕМЕНТА : "<< endl;
		fout<<"Проверка гипотезы о равенстве наибольше отличающихсядисперсий(критерий Фишера): "<< endl;
		fout<<"ДИСПЕРСИЯ МАКСИМАЛЬНАЯ: " << disp_min<< " Количесвтво элементов: "<< x[count_min].length() <<endl; // вывод маскимальной и минимальной дисперсии
		fout<<"ДИСПЕРСИЯ МИНИМАЛЬНАЯ: "<< disp_max<<" Количество элементов: "<< x[count_max].length() << endl;
		for(int i = 0; i < a.size(); i++)
		{
			fout << "При уровне доверия: "<< a[i]<< " квантиль распределение Фишера равен: ";
			fout << invfdistribution(x[count_max].length(),x[count_min].length(), 1 -  a[i])<< endl;

			fout<<"МАТ. ОЖИДАНИЯ МИНИМАЛЬНОЕ: "<<samplemean(x[count_min])<< endl;
			fout<<"МАТ. ОЖИДАНИЕ МАКСИМАЛЬНОЕ: "<<samplemean(x[count_max])<< endl;
			fout << "Проверка гипотезы о равенстве для наиболее отличающихся оценок математических ожиданий (критерий Стьюндента)"<<endl;
			int v;
			double up = pow((disp_max/x[count_max].length()) + (disp_min/x[count_min].length()),2);
			double down_1= pow(disp_max/x[count_max].length(),2) / (x[count_max].length()-1);
			double down_2 = pow(disp_min/x[count_min].length(),2) / (x[count_min].length()-1);

			v = ceil(up/(down_1 + down_2));
			fout<< "		V :" << v<<endl;
			fout <<"Квантиль стьюндента: " << invstudenttdistribution(v,a[i])<< endl;
			if(disp_max/disp_min <= invfdistribution(x[count_max].length(),x[count_min].length(), 1 -  a[i]))
			{
				fout << "Критерий ФИШЕРА ВЕРНО :"<< endl;
				double s;
				s = pow((disp_max*(x[count_max].length()-1) + disp_min*(x[count_min].length()-1))/(x[count_max].length() + x[count_min].length() - 2),0.5);
				fout << "		S:"<< s << endl;
				if((samplemean(x[count_max]) - samplemean(x[count_min]))/s*(pow((x[count_max].length() + x[count_min].length())/x[count_max].length()*x[count_min].length(),0.5))<=invstudenttdistribution(v,a[i]))
					fout << "Гипотеза о равенстве мат.ожиданий серий 1-3 на уровне доверия "<< a[i]<< " принять можно"<<endl;
				else
				{
					fout/*<<endl<<endl*/<<"Критерий Стьюнднта НЕ ВЕРНО"<<endl<<"Надо применять доп критерии"<<" При уровне доверии:"<< a[i]<<endl;
					mnog_met( x,a[i], fout,problam);
					}
			}
			else
			{
				fout << "Критерий ФИШЕРА НЕ ВЕРНО"<<endl;
				if((samplemean(x[count_max]) - samplemean(x[count_min])) / pow(up,1/4) <= invstudenttdistribution(v,a[i]))
					fout << "Критерий Стьюнднта  ВЕРНО"<<endl<< "Гипотеза о равенстве мат.ожиданий серий 1-3 на уровне доверия "<< a[i]<< " принять можно"<<endl;
				else
				{
					fout/*<<endl<<endl*/<<"Критерий Стьюнднта НЕ ВЕРНО"<<endl<< "Надо применять доп критерии"<<" При уровне доверии: "<<a[i]<<endl;
					mnog_met( x,a[i], fout,problam);
				}
			}
		}
		fout << endl <<endl<< "3. ПРОВЕРКА НОРМАЛЬНОСТИ РАСПЕРДЕЛЕНИЕЯ ВЕРОЯТНОСТИ РЕЗУЛЬТАТОВ ИСПЫТАНИЯ"<<endl << endl;
		int n_sum;
		for(auto s : x)
			n_sum += x.length;
		if(n_sum <= 50)
		{
			
		}
	
	}

	return 0;
}

