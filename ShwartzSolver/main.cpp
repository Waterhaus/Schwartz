// ShwartzSolver.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "Relax.h"
#include "LaplasSolver.h"
#include "SaveFile.h"
using namespace std;

int main()
{
	int N = 23;
	vector<int> index;
	Matrix<double> L = CreateLaplasMatrixForCircle(N,&index);
	Vector<double> b = CreateRightPart(&index, N, 0.1);
	Vector<double> x(b.Size());
	Vector<double> anser = Relax(L, b, x, 0.8, 0.0001);
	
	SaveToFile("file1", &anser, N);
	char a;
	cin >> a;
    return 0;
}

