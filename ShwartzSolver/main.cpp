// ShwartzSolver.cpp: ���������� ����� ����� ��� ����������� ����������.
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
	int N = 79;
	vector<int> index;
	Matrix<double> L = CreateLaplasMatrixForCircle(N,&index);

	print_index(index, N);
	Vector<double> b = CreateRightPart(&index, N, 0.1);
	Vector<double> x(b.Size());
	Vector<double> anser = RelaxMini(L, b, x, 1.4, 0.0001,N);
	
	SaveToFile("file1", &anser, N);
	char a;
	cin >> a;
    return 0;
}

