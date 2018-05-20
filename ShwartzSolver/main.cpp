// ShwartzSolver.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "Relax.h"
#include "LaplasSolver.h"
using namespace std;

int main()
{
	Matrix<double> A(3);

	for (size_t i = 0; i < 3; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			A[i][j] = i + j*j + 1;
		}
		A[i][i] = (i+1) * 100;
	}
	
	double mas[3] = { 1, 1, 1 };
	Vector<double> v(mas, 3);
	Vector<double> x0(mas, 3);
	x0[0] = 1;
	Vector<double> Av = A*v;
	//cout << Relax(A,Av,x0,1.2,0.00001);
	Matrix<int> L = CreateLaplasMatrixForCircle(23);
	cout << L;

	char a;
	cin >> a;
    return 0;
}

