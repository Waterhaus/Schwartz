// ShwartzSolver.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "Relax.h"
#include "LaplasSolver.h"
#include "SaveFile.h"
#include "ShwartzSolver.h"
using namespace std;

int main()
{
	
	/*vector<int> index;
	Matrix<double> L = CreateLaplasMatrixForCircle(N,&index);

	print_index(index, N);
	Vector<double> b = CreateRightPart(&index, N, 0.1);
	Vector<double> x(b.Size());
	Vector<double> anser = RelaxMini(L, b, x, 1.4, 0.0001,N);
	
	SaveToFile("file1", &anser, N);
	
	*/
	int N = 80;
	double r = 1.5;
	ShwartzSolver solver(2,2,N,N,0.5,0.5,r);
	solver.SaveSolution();
	/*int N = 7;
	int M = 7;
	Matrix<double> L = Laplas(N, M, 0.1, 0.2);
	
	Vector<double> x0 = Vector<double>::Identity(L.Size());
	x0[0] = 50;
	Vector<double> Lx0 = L*x0;
	Vector<double> temp(Lx0.Size());
	Vector<double> solut = RelaxFast(&L, &Lx0, &temp, 1, 0.000001, M);
	cout << "true err = " <<(solut - x0).norm2();*/
	char a;
	cin >> a;
    return 0;
}

