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
#define ABS(x)  ( (x < 0) ? - (x) : x )
#define MAX_STUDENT 30
#define MAX_SERIA 3
using namespace std;
using namespace alglib;

void gistogram (const map <int,real_1d_array>& x,const real_1d_array& y, fstream& fout);
void mat_interval(const real_1d_array& x, fstream& fout);
real_1d_array static_funk(const vector<double>& number, const vector <double>& a,fstream& fout);
string data_txt(int student_number);
void mnog_met(map <int,real_1d_array>& x, double a, fstream& fout, fstream& problam);
void stat_touchnost(const double& a,const double& k, fstream& fout);
void norm_raspredelenie(const real_1d_array& x, fstream& fout);
void imperial_sakon(const real_1d_array& x, fstream& fout,const map<double, double>& l);
const map<double, double > lamda ={{0, 1},{0.1, 1},{0.2, 1},{0.3, 1},{0.4, 0.997},{0.5, 0.964},{0.6, 0.864},{0.7, 0.711},{0.8, 0.544},{1, 0.27},{1.1, 0.178},{1.2, 0.112},{1.3, 0.068},{1.4, 0.04},{1.5, 0.022},{1.6, 0.012},{1.7, 0.006},{1.8, 0.003},{1.9, 0.002},{2, 0.001}};
void talerant_interval(double a, int n, double disp, double mat, fstream& fout);
const int N[] = {16,27,28};
const vector<double> a = {0.9,0.95,0.99};


int main()
{
	fstream problam;
	string n;
	problam.open("problam.txt", ios::out);
	fstream fout;
	fout.open("reshenie.txt", ios ::out);
	for(int k = 2; k <=2; k++)
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
		double mat_min =samplemean (x[1]);
		ae_int_t count_min = 1;
		double mat_max = samplemean (x[1]);
		ae_int_t count_max = 1;
		int i = 1;
		while(i <= x.size())
		{
			if(samplemean(x[i]) > mat_max)
			{
				mat_max = samplemean (x[i]);
				count_max = i;
			}
			if(samplemean (x[i]) < mat_min)
			{
				mat_min = samplemean (x[i]);
				count_min = i;
			}
			i++;
		}

		fout<< endl<< endl;
		fout<<"2 ПРОВЕРКА ОДНОРОДНОСТИ РЕЗУЛЬТАТОВ ЭКСПЕРЕМЕНТА : "<< endl;
		fout<<"Проверка гипотезы о равенстве наибольше отличающихсядисперсий(критерий Фишера): "<< endl;
		fout<<"МАТ. ОЖИДАНИЕ МАКСИМАЛЬНАЯ: " << mat_max<< " Колич	есвтво элементов: "<< x[count_min].length() <<endl; // вывод маскимальной и минимальной дисперсии
		fout<<"МАТ. ОЖИДАНИЕ МИНИМАЛЬНАЯ: "<< mat_min<<" Количество элементов: "<< x[count_max].length() << endl;
		for(int i = 0; i < a.size(); i++)
		{
			int v;
			double disp_min = samplevariance(x[count_min]);
			double disp_max = samplevariance(x[count_max]);
			double up = pow((disp_max/x[count_max].length()) + (disp_min/x[count_min].length()),2);
			double down_1= pow(disp_max/x[count_max].length(),2) / (x[count_max].length()-1);
			double down_2 = pow(disp_min/x[count_min].length(),2) / (x[count_min].length()-1);
			v = ceil(up/(down_1 + down_2));
			double up_stud = ABS(mat_max - mat_min);
			double down_stud_true1 = x[count_max].length() + x[count_min].length();
			double down_stud_true2 = x[count_max].length() * x[count_min].length();
			fout<<"ДИСПЕРСИЯ МИНИМАЛЬНОЕ: "<<samplevariance(x[count_min])<< endl;
			fout<<"ДИСПЕРСИЯ МАКСИМАЛЬНОЕ: "<<samplevariance(x[count_max])<< endl;
			fout << "При уровне доверия: "<< a[i]<< " квантиль распределение Фишера равен: ";
			fout << invfdistribution(x[count_max].length(),x[count_min].length(), 1 -  a[i])<< endl;			
			fout << "Отношение дисперсий: "<< disp_max/disp_min<<endl;
			bool flag;
			if(disp_max/disp_min <= invfdistribution(x[count_max].length() - 1,x[count_min].length(), 1 -  a[i]))
			{
				flag = true;
				fout << "Критерий ФИШЕРА ВЕРНО :"<< endl;
			}
			else 
			{
				flag = false;
				fout << "Критерий ФИШЕРА НЕ ВЕРНО"<<endl;
			}
			fout << "Проверка гипотезы о равенстве для наиболее отличающихся оценок математических ожиданий (критерий Стьюндента)"<<endl;			
			
			if(flag)
			{
				
				double s;
				fout<< "		V :" << v<<endl;
				fout <<"Квантиль стьюндента: " << invstudenttdistribution(v,a[i])<< endl;
				s = disp_max*(x[count_max].length()-1) + disp_min*(x[count_min].length()-1);
				s /= x[count_max].length() + x[count_min].length() - 2;
				s = pow(s,0.5);
				fout << "		S:"<< s << endl;
				//cout << 
				fout<< "Левая часть результат вычислений" << up_stud / s / pow(down_stud_true1/ down_stud_true2, 0.5)<<endl;		
				if( up_stud / s / pow(down_stud_true1/ down_stud_true2, 0.5)<=invstudenttdistribution(v,a[i]))
					fout << "Гипотеза о равенстве мат.ожиданий серий 1-3 на уровне доверия "<< a[i]<< " принять можно"<<endl;
				else
				{
					fout/*<<endl<<endl*/<<"Критерий Стьюнднта НЕ ВЕРНО"<<endl<<"Надо применять доп критерии"<<" При уровне доверии:"<< a[i]<<endl;
					mnog_met( x,a[i], fout,problam);
					}
			}
			else
			{
				
				fout<< "		V :" << v<<endl;
				fout <<"Квантиль стьюндента: " << invstudenttdistribution(v,a[i])<< endl;
				fout<< "Левая часть результат вычислений" << (samplemean(x[count_max]) - samplemean(x[count_min])) / pow(up,1/4)<<endl;
				if((samplemean(x[count_max]) - samplemean(x[count_min])) / pow(up,1/4) <= invstudenttdistribution(v,a[i]))
					fout << "Критерий Стьюнднта  ВЕРНО"<<endl<< "Гипотеза о равенстве мат.ожиданий серий 1-3 на уровне доверия "<< a[i]<< " принять можно"<<endl;
				else
				{
					fout/*<<endl<<endl*/<<"Критерий Стьюнднта НЕ ВЕРНО"<<endl<< "Надо применять доп критерии"<<" При уровне доверии: "<<a[i]<<endl;
					mnog_met( x,a[i], fout,problam);
				}
			}
		}
		vector<double> vec;
		for(int i = 1 ; i <= MAX_SERIA; i++ )
			for(int j = 0; j < x[i].length(); j++)
				vec.push_back(x[i][j]);
		real_1d_array y;
		y.setcontent(vec.size(),&(vec[0]));
		fout << "РАССЧИТЫВАЕМ ТОЧНОСТЬ СТАТИСТИЧЕСКИХ РЕШЕНИЙ ДЛЯ ОШТБОК 1 И 2 РОДА И КОЛИЧЕСТВО ИЗМЕРЕНИЙ: "<< endl;
		fout << " Постоить график зависимости точности статисческих решений от объема исследованной выборки" << endl;
		int sum_seria = 0;
		for(int i = 1; i <=MAX_SERIA;i++)
			sum_seria += x[i].length();
		fout <<endl<<  "Для общей выборки"<<endl;
		for(auto s : a)
			stat_touchnost(s,sum_seria,fout);
		fout<<endl<<endl<<endl<< "КРИТЕРИЙ ШАПИРА-УИЛКА"<<endl;
		for(int i = 1; i <= MAX_SERIA; i++)
		{
			fout<<endl<< "Расчет для серии "<< i << endl<<endl;
			norm_raspredelenie(x[i],fout);
		}
		fout <<endl<<  "Для общей выборки"<<endl;
		norm_raspredelenie(y,fout);
		fout << "Постоение Имперического закона распееделения эксперементальных данных"<<endl << endl << endl;
		for(int i = 1; i <= MAX_SERIA; i++)
		{
			fout<<endl<< "Данные для построение графика серии "<< i << endl<<endl;
			imperial_sakon(x[i],fout,lamda);
		}
		fout << "Для общей выборки" <<endl<<endl; 
		imperial_sakon(y,fout,lamda);
		fout <<endl<< "Построене талератного интервалов "<< endl << endl;
		for(int i = 1; i <= MAX_SERIA; i++)
		{
			fout << "Для серии " << i <<endl;
			talerant_interval(0.95, x[i].length(),samplevariance(x[i]) , samplemean(x[i]), fout);	
			fout << endl<<endl<<endl;
		}			
		fout << "Для общей выборки" <<endl<<endl; 
		talerant_interval(0.95, y.length(),samplevariance(y) , samplemean(y), fout);	
		mat_interval(y,fout);
		gistogram (x,y,fout);
		fout.close();
		fout.open("reshenie.txt");
		fstream fin;
		fin.open("reshenie_new.txt",ios ::out);
		while(!fout.eof())
		{
			string str;
			getline(fout,str);
			str+='\n';
			for(auto& s : str)
			{
				if(s == '.')
					s = ',';
			}
			fin <<str;
						
		}
		fin.close();
		fout.close();
	}

	return 0;
}
