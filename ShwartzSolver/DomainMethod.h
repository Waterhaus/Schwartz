#pragma once
#include <vector>
#include <iostream>
#include <string>
#include "LaplasSolver.h"
#include "Relax.h"
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

double F(double x, double y)
{
	return 0;
}

struct IndexPair {

	int index_i;
	int index_j;

};


class BoxDomain 
{
public:
	Vector<double> u;
	Vector<double> b;
	vector<IndexPair> border_in;
	vector<IndexPair> border_out;
	vector<double> out_x;
	vector<double> out_y;

	Matrix<double> L;

	double x_start;
	double y_start;

	double x_end;
	double y_end;

	int SIZE;
	double h;
	string name;
	

	BoxDomain(	string b_name,int N,double X_start, double Y_start,
				double X_end, double Y_end, double Radius)
	{
		name = b_name;
		SIZE = N;
		h = GetStep(SIZE, X_start, X_end);
		double hy = GetStep(SIZE, Y_start, Y_end);
		if (abs(h - hy) > 0.001) cout << "fail of box borders " + name << endl;

		FillRightPart();
		
		SetIndex(Radius);
		FillBorderCoords();
	}

	void FillRightPart()
	{
		b.Resize(SIZE*SIZE);
		int s = 0;
		for (size_t i = 0; i < SIZE; i++)
		{
			for (size_t j = 0; j < SIZE; j++)
			{
				s = i*SIZE + j;
				
				if (i == 0 || j == 0 || j == SIZE - 1 || i == SIZE - 1)
					b[s] = 0; 
				else
					b[s] = F(x_start + j*h,y_start + i*h);
				
			}
		}
	
	}

	//заполняем два массива, отвечающих за границу области.
	// один за индексы внутри круга, другой снаружи 
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

	//заполняем массив х и у - точек на границе
	void FillBorderCoords()
	{
		out_x.resize(border_out.size());
		out_y.resize(border_out.size());
		for (size_t i = 0; i < out_x.size(); i++)
		{
			out_x[i] = x_start + h*border_out[i].index_j;
			out_y[i] = y_start + h*border_out[i].index_i;
		}
	}

	//получаем массив значений, которые устанавливаем на границе(заносим в правую часть)
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

	//определяем индекс i*h <= val <= (i+1)*h
	int FindIndex(double val, double val_start)
	{
		return (int)floor((val - val_start) / h);
	}

	//интерполируем значения в любой точке внутри сетки
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
		cout << name + " - interpolating done";
	}

	void CorrectLaplasMatrix()
	{
		int s = 0;

		for (size_t k = 0; k < border_in.size(); k++)
		{
			s = border_in[k].index_i*SIZE + border_in[k].index_j;

			for (size_t i = 0; i < SIZE*SIZE; i++)
			{
				L[s][i] = 0;
			}
			L[s][s] = 1;
		}

		for (size_t k = 0; k < border_out.size(); k++)
		{
			s = border_out[k].index_i*SIZE + border_out[k].index_j;

			for (size_t i = 0; i < SIZE*SIZE; i++)
			{
				L[s][i] = 0;
			}
			L[s][s] = 1;
		}

	}

	void SolveLaplas(double w, double EPS)
	{
		if (L.Size() != SIZE*SIZE) 
		{
			L.Resize(SIZE);
			Laplas(&L, SIZE, SIZE, h, h);
			u.Resize(SIZE);
			CorrectLaplasMatrix();
		}
		Vector<double> x0 = u;
		RelaxFast(&L, &b, &x0, &u, w, EPS,SIZE);
		cout << "new solve complete for " + name << endl;
	}
};
