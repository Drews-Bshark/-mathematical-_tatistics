/*
 * Odnorodnist.cpp
 *
 *  Created on: 17 мар. 2019 г.
 *      Author: andrey
 */

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
# define ABS(x)  ( (x < 0) ? -(x) : x )
using namespace std;
const int N[] = {16,27,28};


string data_txt(int student_number)
{
	string name;
	fstream NAME;
	NAME.open("name.txt");
	int i = 0;
	while(!NAME.eof())
	{
		getline(NAME,name);
		if(i == student_number - 1)
			break;
		i++;
	}
	NAME.close();
	int count = 0;
	vector<string> file = {"seria1.txt","seria2.txt","seria3.txt"};
	map <string,vector<double>> number;
	fstream data;
	data.open("data.txt",ios::out);
	while(count < file.size())
	{
		fstream fopen;
		fopen.open(file[count]);
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
		data<< str[student_number - 1]<<endl;
		count++;
	}
	return (name);

}


