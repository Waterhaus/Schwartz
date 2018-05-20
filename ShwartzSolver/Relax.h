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
		cout << (x - x0).norm2() << endl;
		//cout << "x[" << iter << "] = " << x << endl;
	}

	cout << "solve with w = " << w << "; x = " << x << "; " << iter <<" iterations." << endl;
	return x;
}

Vector<double> RelaxMini(Matrix<double> A, Vector<double> b, Vector<double> x0, double w, double EPS,int N)
{

	Vector<double> x = A*x0;

	double S = 0;
	int iter = 0;
	int k = 0;
	while (iter == 0 || (x - x0).norm2() > EPS)
	{

		for (size_t i = 0; i < x.Size(); i++)
		{
			S = 0;
			if (i > N + 1 && i < x.Size() - N - 2) k = N + 1;
			else k = 0;
			for (size_t j = i - k; j < i + k; j++)
			{
				if (i != j) S = S + A[i][j] * x[j];
			}
			//cout << "S = " << S << endl;
			x0[i] = x[i];
			x[i] = (1.0 - w)*x[i] + (w / A[i][i])*(b[i] - S);

		}
		iter++;
		cout << (x - x0).norm2() << endl;
		//cout << "x[" << iter << "] = " << x << endl;
	}

	cout << "solve with w = " << w << "; x = " << x << "; " << iter << " iterations." << endl;
	return x;
}
