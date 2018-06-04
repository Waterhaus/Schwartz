﻿#pragma once
#include <iostream>
#include "Matrix.h"


using namespace std;
Vector<double> Relax(Matrix<double> A, Vector<double> b, Vector<double> x0, double w, double EPS)
{

	Vector<double> x = A*x0;
	
	double S = 0;
	int iter = 0;
	while(iter == 0 || (x - x0).norm2() > EPS )
	{
		
		for (size_t i = 0; i < x.Size(); i++)
		{
			S = 0;
			for (size_t j = 0; j < x.Size(); j++)
			{
				if (i != j) S = S + A[i][j] * x[j];
			}
			//cout << "S = " << S << endl;
			x0[i] = x[i];
			x[i] = (1.0 - w)*x[i] + (w / A[i][i])*(b[i] - S);

		}
		iter++;
		cout << "||x - x_prev|| = " << (x - x0).norm2() << '\r';
		//cout << "x[" << iter << "] = " << x << endl;
	}
	cout << endl;
	cout << "solve with w = " << w << "; x = " << x << "; " << iter <<" iterations." << endl;
	return x;
}

void RelaxMini(Matrix<double> *A, Vector<double> *b, Vector<double> *x0, double w, double EPS, int N)
{

	Vector<double> x = (*A)*(*x0);

	

	cout << (*A) << endl;
	cout << "Ax0 = " << x << endl;
	char aaa;
	cin >> aaa;
	double S = 0;
	int iter = 0;
	int k = 0;
	while (iter == 0 || (x - (*x0)).norm2() > EPS)
	{

		for (size_t i = 0; i < x.Size(); i++)
		{
			S = 0;
			if (i > N + 1 && i < x.Size() - N - 2) k = N + 1;
			else k = 0;

			for (size_t j = i - k; j < i + k; j++)
			{
				if (i != j) S = S + (*A)[i][j] * x[j];
			}
			//cout << "S = " << S << endl;
			(*x0)[i] = x[i];
			x[i] = (1.0 - w)*x[i] + (w / (*A)[i][i])*((*b)[i] - S);

		}
		iter++;
		cout << (x - (*x0)).norm2() << endl;
		cout << "x[" << iter << "] = " << x << endl;
	}

	cout << "solve with w = " << w << "; x = " << x << "; " << iter << " iterations." << endl;
}



void RelaxFast(Matrix<double> *A, 
	Vector<double> *b,
	Vector<double> *x0,
	Vector<double> *solution, double w, double EPS, int N)
{

	Vector<double> x = (*A)*(*x0);
	cout  << "b = "<< (*b) << endl;
	//cout << "Ax0 = " << x << endl;

	
	double S = 0;
	int iter = 0;
	int k = 0;
	while (iter == 0 || (x - (*x0)).norm2() > EPS)
	{

		for (size_t i = 0; i < x.Size(); i++)
		{
			S = 0;	
			//if (i < x.Size() - N - 1 && i > N)
			//cout << i <<")  " << (*A)[i][i - N] << "  " << (*A)[i][i - 1] << "  " << (*A)[i][i] << "  " << (*A)[i][i + 1] << "  " << (*A)[i][i + N] << endl;

			if (i > 0) S = S + (*A)[i][i - 1] * x[i - 1]; 
			if (i < x.Size() - 1) S = S + (*A)[i][i + 1] * x[i + 1];
			if (i >= N) S = S + (*A)[i][i - N] * x[i - N];
			if (i < x.Size() - N ) S = S + (*A)[i][i + N] * x[i + N];

			//cout << "S = " << S << endl;
			(*x0)[i] = x[i];

			

			x[i] = (1.0 - w)*x[i] + (w / (*A)[i][i])*((*b)[i] - S);
			
		}
		iter++;
		
		cout << "||x - x_prev|| = " << (x - (*x0)).norm2() << '\r';
		//cout << "x[" << iter << "] = " << x << endl;
	}
	for (size_t i = 0; i < x.Size(); i++)
	{
		(*solution)[i] = x[i];
	}
	cout << endl;
	cout << "solve with w = " << w << "; x = " << (*solution) << "; " << iter << " iterations." << endl;
	
}



void RelaxFastBad(Matrix<double> *A,
	Vector<double> *b,
	Vector<double> *x0,
	Vector<double> *solution, double w, double EPS, int N)
{

	Vector<double> x = (*A)*(*x0);
	cout << "b = " << (*b) << endl;
	cout << "Ax0 = " << x << endl;
	x = Vector<double>::Identity(x.Size());

	double S = 0;
	int iter = 0;
	int k = 0;
	while (iter == 0 || (x - (*x0)).norm2() > EPS)
	{
		int i = 0;
		int j = 0;
		for (size_t s = 0; s < x.Size(); s++)
		{

			if (j == N) j = 0;
			i = (s - j + 0.0) / (double)N;
			S = 0;
			S = S + (*A)[s][s - 1] * x[s - 1];
			S = S + (*A)[s][s + 1] * x[s + 1];
			S = S + (*A)[s][s - N] * x[s - N];
			S = S + (*A)[s][s + N] * x[s + N];
			S = S + (*A)[s][N*(N - 2) + j] * x[N*(N - 2) + j];
			S = S + (*A)[s][N*(N - 1) + j] * x[N*(N - 1) + j];
			S = S + (*A)[s][j]*x[]
			
			//cout << "S = " << S << endl;
			(*x0)[s] = x[s];


		
			x[s] = (1.0 - w)*x[s] + (w / (*A)[s][s])*((*b)[s] - S);
			j++;
		}
		iter++;
		cout << "||x - x_prev|| = " << (x - (*x0)).norm2() << '\r';
		//cout << "x[" << iter << "] = " << x << endl;
	}
	for (size_t i = 0; i < x.Size(); i++)
	{
		(*solution)[i] = x[i];
	}
	cout << "RelaxBad:  solve with w = " << w << "; x = " << (*solution) << "; " << iter << " iterations." << endl;

}