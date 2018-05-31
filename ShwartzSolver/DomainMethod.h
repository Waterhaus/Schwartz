#pragma once
#include <vector>
#include <iostream>
#include <string>
using namespace std;



double GetStep(int GridSize, double a_border, double b_border)
{

	if (a_border > b_border)
	{
		double temp = a_border;
		a_border = b_border;
		b_border = temp;
	}

	return (b_border - a_border) / (GridSize - 1);
}

double B(double x, double xi, double h)
{
	if (abs(x - xi) > h) return 0;
	double temp = (x - xi) / h;

	if (x - xi <= 0) return temp + 1;
	else return -temp + 1;
	
}



struct IndexPair {

	int index_i;
	int index_j;

};


class BoxDomain 
{
public:
	vector<double> u;
	vector<double> b;
	vector<IndexPair> border_in;
	vector<IndexPair> border_out;
	vector<double> out_x;
	vector<double> out_y;


	double x_start;
	double y_start;

	double x_end;
	double y_end;

	int SIZE;
	double h;
	string name;
	

	void SetIndex(double R)
	{
		//i == 0; -> y = y_start;
		double x_temp = 0;
		double y_temp = 0;
		double mod;
		IndexPair index;
		for (size_t j = 0; j < SIZE; j++)
		{
			x_temp = x_start + j*h;
			mod = sqrt( y_start*y_start + x_temp*x_temp );
			index.index_i = 0;
			index.index_j = j;
			if (mod <= R) border_in.push_back(index);
			else border_out.push_back(index);
		}

		//i == SIZE - 1 -> y = y_end
		for (size_t j = 0; j < SIZE; j++)
		{
			x_temp = x_start + j*h;
			mod = sqrt(y_end*y_end + x_temp*x_temp);
			index.index_i = SIZE - 1;
			index.index_j = j;
			if (mod <= R) border_in.push_back(index);
			else border_out.push_back(index);
		}

		//j == 0 - > x = x_start
		for (size_t i = 0; i < SIZE; i++)
		{
			y_temp = y_start + i*h;
			mod = sqrt(y_temp*y_temp + x_start*x_start);
			index.index_i = i;
			index.index_j = 0;
			if (mod <= R) border_in.push_back(index);
			else border_out.push_back(index);
		
		}

		//j == SIZE - 1 -> x = x_end
		for (size_t i = 0; i < SIZE; i++)
		{
			y_temp = y_start + i*h;
			mod = sqrt(y_temp*y_temp + x_end*x_end);
			index.index_i = i;
			index.index_j = SIZE - 1;
			if (mod <= R) border_in.push_back(index);
			else border_out.push_back(index);

		}
	}

	void FillBorderDots()
	{
		out_x.resize(border_out.size());
		out_y.resize(border_out.size());
		for (size_t i = 0; i < out_x.size(); i++)
		{
			out_x[i] = x_start + h*border_out[i].index_j;
			out_y[i] = y_start + h*border_out[i].index_i;
		}
	}

	void UpdateBorder(double* border,int size)
	{
		int s = 0;
		int i = 0;
		int j = 0;
		for (size_t k = 0; k < size; k++)
		{
			i = border_out[k].index_i;
			j = border_out[k].index_j;
			s = i*SIZE + j;
			b[s] = border[k];
		}
		cout << name + " - border updated" << endl;
	}

	int FindIndex(double val, double val_start)
	{
		return (int)floor((val - val_start) / h);
	}

	void GetInterpolatedPoints(double* x, double* y, double* g, int size)
	{
		int i = 0;
		int j = 0;
		int s1, s2, s3, s4;
		double S = 0;
		for (size_t k = 0; k < size; k++)
		{
			i = FindIndex(y[k], y_start);
			j = FindIndex(x[k], x_start);

			s1 = i*SIZE + j;
			s2 = (i + 1)*SIZE + j;
			s3 = i*SIZE + j + 1;
			s4 = (i + 1)*SIZE + j + 1;

			g[k] = u[s1] * B(x[k], x_start + j*h, h)*B(y[k], y_start + i*h, h) +
				u[s2] * B(x[k], x_start + j*h, h)*B(y[k], y_start + (i + 1)*h, h) +
				u[s3] * B(x[k], x_start + (j + 1)*h, h)*B(y[k], y_start + i*h, h) +
				u[s4] * B(x[k], x_start + (j + 1)*h, h)*B(y[k], y_start + (i + 1)*h, h);

		}
		cout << name + " - interpolatin done";
	}

};
