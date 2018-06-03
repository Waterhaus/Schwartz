#pragma once
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

		for (size_t i = 0; i < x.Size(); i++)
		{


			S = 0;
			//if (i < x.Size() - N - 1 && i > N)
			//cout << i <<")  " << (*A)[i][i - N] << "  " << (*A)[i][i - 1] << "  " << (*A)[i][i] << "  " << (*A)[i][i + 1] << "  " << (*A)[i][i + N] << endl;
			
			if(i < 2*N|| i > x.Size() - 2*N)
			{
				if( i < 2*N )
				{
					int index = N*(N - 1) + i - N;
					if (i > 0)S = S + (*A)[i][i - 1] * x[i - 1];
					if (i < x.Size() - 1) S = S + (*A)[i][i + 1] * x[i + 1];
					if (i >= N) S = S + (*A)[i][index] * x[index];
					//cout << "i<N  " << index << "  " <<(*A)[i][index] << endl;
					if (i < x.Size() - N) S = S + (*A)[i][i + N] * x[i + N];

				}
				
				if(i > x.Size() - 2*N)
				{
					int index = i - N*(N - 2);
					if(i > x.Size() - N) index = i - N*(N - 2);
					if (i > 0)S = S + (*A)[i][i - 1] * x[i - 1];
					if (i < x.Size() - 1) S = S + (*A)[i][i + 1] * x[i + 1];
					if (i >= N) S = S + (*A)[i][i - N] * x[i - N];
					if (i < x.Size() - N) S = S + (*A)[i][index] * x[index];
					//cout << "i in end  " << (*A)[i][index] << endl;

				}
			}
			else
			{
				if (i > 0) S = S + (*A)[i][i - 1] * x[i - 1];
				if (i < x.Size() - 1) S = S + (*A)[i][i + 1] * x[i + 1];
				if (i >= N) S = S + (*A)[i][i - N] * x[i - N];
				if (i < x.Size() - N) S = S + (*A)[i][i + N] * x[i + N];
			
			}
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
	cout << "RelaxBad:  solve with w = " << w << "; x = " << (*solution) << "; " << iter << " iterations." << endl;

}