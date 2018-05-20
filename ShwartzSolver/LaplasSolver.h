#pragma once
#include "Matrix.h"
#include <vector>

 Matrix<double> CreateLaplasMatrix(int N)
{
	int dim = (int)pow((N - 1), 2); //Размерность матрицы
	Matrix<double> A(dim);
	int n = N - 1;

	int k = 0;
	for (int i = 0; i < dim; i++)
	{
		A[i][i] = 4;

		if (k != n - 1 && i + 1 != dim) A[i + 1][ i] = -1;
		if (k != n && i - 1 >= 0) A[i - 1][ i] = -1;


		if ((i + n - 1) < dim) A[i][ i + n - 1] = -1;
		if ((i + n - 1) < dim) A[i + n - 1][ i] = -1;


		if (k >= n) k = 0;
		k++;


	}

	return A;
}

 Matrix<int> CreateLaplasMatrixForCircle(int N)
 {
	 Matrix<int> L(N); //= CreateLaplasMatrix(N); 
	 
	 Matrix<double> M = CreateLaplasMatrix(N);
	 vector<int> s_index(M.Size());
	 int n = N / 2;
	 double R = (double)N / 2.0;
	 double I = 0;
	 int index = 0;
	 for (size_t i = 0; i < N; i++)
	 {
		 for (size_t j = 0; j < N; j++)
		 {
			 I = sqrt((double)((i - n)*(i - n)) + (double)((j - n)*(j - n)));
			 if (I <= R) 
			 { 
				 L[i][j] = 1; 
			 }
			 else L[i][j] = 0;

			 if ((i == 0 || i == N - 1) && L[i][j] == 1) 
				 L[i][j] = 2;
			 if ((j == 0 || j == N - 1) && L[i][j] == 1) 
				 L[i][j] = 2;

			 if( j > 0  )
			 {
				 if (L[i][j] == 1 && L[i][j - 1] == 0) L[i][j] = 2;
				 if (L[i][j] == 0 && L[i][j - 1] == 1) L[i][j - 1] = 2;
			
			 }

			 s_index[index] = L[i][j];
			 index++;
			
			
		 }
	 }

	 return L;
 }

 /*
 Matrix<double> CreateMapForCircle(int N)
 {
	 vector<vector<int>> map;
	 vector<int> data;

	 int n = N / 2;
	 double R = (double)N / 2.0;
	 double I = 0;
	 double I_old = 0;
	 double I_new = 0;
	 int isGamma;
	 int s = 0;
	 for (int i = 0; i < N; i++)
	 {
		 for (int j = 0; j < N; j++)
		 {
			 s++;
			 I = sqrt((double)((i - n)*(i - n)) + (double)((j - n)*(j - n)));
			 if (I >= R)
			 {
				 data = {s,,i,j};
			 }
		 }
	 }

	 return L;
 }
*/