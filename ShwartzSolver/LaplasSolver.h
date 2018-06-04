#pragma once
#include "Matrix.h"
#include <vector>

double f(double x, double y)
{
	return 0;
}
double g(double x, double y)
{
	return 10*sin(x) + 20;
}

Matrix<double> Laplas(int N1,  double hx, double hy)
{

	int dim = N1 * N1; //–азмерность матрицы
	Matrix<double> A(dim);



	int i = 0;
	for (int s = 0; s < dim; s++)
	{

		//проверка дл€ j 
		if (i == N1) i = 0;
		// найдем теперь i
		int j = (int)(((double)(s - i)) / ((double)N1));

		A[s][ s] = 2 * (1.0f / (hx*hx) + 1.0f / (hy*hy));

		if (i < N1 - 1) A[s + 1][s] = -1.0f / (hx*hx);
		if (i > 0) A[s - 1][s] = -1.0f / (hx*hx);


		if (j < N1 - 1) A[s + N1][s] = -1.0f / (hy*hy);
		if (j > 0) A[s - N1][s] = -1.0f / (hy*hy);


		i++;
	}

	return A;

}



Matrix<double> Laplas(int N, int M, double hx, double hy)
{

	int dim = N * M; //–азмерность матрицы
	Matrix<double> A(dim);


	

	int i = 0;
	for (int s = 0; s < dim; s++)
	{
		if (s == dim - 1)
		{
			int c = 0;
		}
		//проверка дл€ j 
		if (i == N) i = 0;
		// найдем теперь i
		int j = (int)(((double)(s - i)) / ((double)N));

		A[s][ s] = 2 * (1.0f / (hx*hx) + 1.0f / (hy*hy));

		if (i < N - 1) 
			A[s][s + 1] = -1.0f / (hx*hx);
		if (i > 0)
			A[s ][s - 1] = -1.0f / (hx*hx);


		if (j < M - 1) 
			A[s][s + M] = -1.0f / (hy*hy);
		if (j > 0) 
			A[s][s - M] = -1.0f / (hy*hy);

		cout << A[s][s] << endl;
		i++;
	}

	return A;

}


void Laplas(Matrix<double> *A, int N, int M, double hx, double hy)
{

	int dim = N * M; //–азмерность матрицы

	int i = 0;
	for (int s = 0; s < dim; s++)
	{
		
		//проверка дл€ j 
		if (i == N) i = 0;
		// найдем теперь i
		int j = (int)(((double)(s - i)) / ((double)N));

		(*A)[s][s] = 2 * (1.0f / (hx*hx) + 1.0f / (hy*hy));

		if (i < N - 1)
			(*A)[s][s + 1] = -1.0f / (hx*hx);
		if (i > 0)
			(*A)[s][s - 1] = -1.0f / (hx*hx);


		if (j < M - 1)
			(*A)[s][s + N] = -1.0f / (hy*hy);
		if (j > 0)
			(*A)[s][s - N] = -1.0f / (hy*hy);

		
		i++;
	}

	

}



void LaplasPolar(Matrix<double> *A, int N, int M, double hr, double hphi)
{

	int dim = N * M; //–азмерность матрицы

	int j = 0;
	for (int s = M; s < dim - M; s++)
	{

		//проверка дл€ j 
		if (j == M) j = 0;
		// найдем теперь i
		int i = (int)((double)(s - j) / (double)M);

		(*A)[s][s] = -2.0*( 1.0/(hr*hr) + 1.0/(i*i*hr*hr*hphi*hphi));

		if (i < N - 1)
			(*A)[s][s + 1] = (1.0 / (hr*hr) - 1.0 /(2.0*i*hr*hr) );
		if (i > 0)
			(*A)[s][s - 1] = (1.0 / (hr*hr) + 1.0 / (2.0*i*hr*hr));


		if (j < M - 1)
			(*A)[s][s + N] = 1.0/(i*i*hr*hr*hphi*hphi);
		if (j > 0)
			(*A)[s][s - N] = 1.0 / (i*i*hr*hr*hphi*hphi);

		j++;
	}



}


void LaplasPolar2(Matrix<double> *A, int N, int M, double hr, double hphi)
{

	int dim = N * M; //–азмерность матрицы
	int s = 0;
	
	for (size_t i = 2; i < N - 2; i++)
	{
		for (size_t j = 1; j < M - 1; j++)
		{
			s = i*M + j;

			(*A)[s][s] = -2.0*(1.0 / (hr*hr) + 1.0 / (j*j*hr*hr*hphi*hphi));
			(*A)[s][s + 1] = (1.0 / (hr*hr) + 1.0 / (2.0*j*hr*hr));
			(*A)[s][s - 1] = (1.0 / (hr*hr) - 1.0 / (2.0*j*hr*hr));
			(*A)[s][s + N] = 1.0 / (j*j*hr*hr*hphi*hphi);
			(*A)[s][s - N] = 1.0 / (j*j*hr*hr*hphi*hphi);

			

		}
	}


}
/*

A[s][s] = -2.0 * (1.0 / ( h_r * h_r ) + 1.0 / ( ri * ri * h_phi * h_phi ) );

A[s][s + 1] = 1.0 / ( h_r * h_r ) - 1.0 / ( 2.0 * ri * h_r * h_r );
A[s][s - 1] = 1.0 / ( h_r * h_r ) + 1.0 / ( 2.0 * ri * h_r * h_r );

A[s][GRID_DIM_N + j] = 1.0 / ( ri * ri * h_phi * h_phi );
A[s][s - GRID_DIM_N] = 1.0 / ( ri * ri * h_phi * h_phi );
*/
void LaplasPolar3(Matrix<double> *A, int N, int M, double hr, double hphi)
{

	int dim = N * M; //–азмерность матрицы
	int s = 0;
	double ri = 0;
	double h_phi = hphi;
	double h_r = hr;
	for (size_t i = 2; i < N - 2; i++)
	{
		for (size_t j = 1; j < M - 1; j++)
		{
			s = i*M + j;
			ri = j*hr;
			(*A)[s][s] = -2.0 * (1.0 / (h_r * h_r) + 1.0 / (ri * ri * h_phi * h_phi));
			(*A)[s][s + 1] = 1.0 / (h_r * h_r) - 1.0 / (2.0 * ri * h_r * h_r);
			(*A)[s][s - 1] = 1.0 / (h_r * h_r) + 1.0 / (2.0 * ri * h_r * h_r);
			(*A)[s][s + N] = 1.0 / (ri * ri * h_phi * h_phi);
			(*A)[s][s - N] = 1.0 / (ri * ri * h_phi * h_phi);



		}
	}


}



 Matrix<double> CreateLaplasMatrix(int N)
{
	int dim = (int)pow((N - 1), 2); //–азмерность матрицы
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


 void LaplasMatrixForCircle(Matrix<double>* L,int N_x,int M_y, vector<int>* s_index, double hx, double hy)
 {
	 Matrix<int> C(N_x,M_y); //= CreateLaplasMatrix(N); 

	 Laplas(L,N_x,M_y, hx, hy);
	 (*s_index).resize(L->Size());
	 int n = N_x / 2;
	 double R = (double)N_x / 2.0;
	 double I = 0;
	 int index = 0;
	 for (size_t i = 0; i < N_x; i++)
	 {
		 for (size_t j = 0; j < M_y; j++)
		 {
			 I = sqrt((double)((i - n)*(i - n)) + (double)((j - n)*(j - n)));
			 if (I <= R)
			 {
				 C[i][j] = 1;
				 (*s_index)[index] = 1;
			 }
			 else {
				 C[i][j] = 0;
				 (*s_index)[index] = 0;
			 }

			 if ((i == 0 || i == N_x - 1) && C[i][j] == 1)
			 {
				 C[i][j] = 2;
				 (*s_index)[index] = 2;
			 }
			 if ((j == 0 || j == N_x - 1) && C[i][j] == 1)
			 {
				 C[i][j] = 2;
				 (*s_index)[index] = 2;
			 }

			 if (j > 0)
			 {
				 if (C[i][j] == 1 && C[i][j - 1] == 0)
				 {
					 C[i][j] = 2;
					 (*s_index)[index] = 2;
				 }
				 if (C[i][j] == 0 && C[i][j - 1] == 1) {
					 C[i][j - 1] = 2;
					 (*s_index)[index - 1] = 2;
				 }

			 }


			 index++;


		 }
	 }

	 for (size_t i = 0; i < (*s_index).size(); i++)
	 {
		 if ((*s_index)[i] == 0 || (*s_index)[i] == 2)
		 {
			 for (size_t j = 0; j < (*s_index).size(); j++)
			 {
				 if (i != j) (*L)[i][j] = 0;
				 else (*L)[i][j] = 1;
			 }
		 }
	 }
	 // cout << C;
	 
 }


 Matrix<double> CreateLaplasMatrixForCircle(int N, vector<int>* s_index)
 {
	 Matrix<int> C(N); //= CreateLaplasMatrix(N); 
	 
	 Matrix<double> L = Laplas(N,1,1);
	 (*s_index).resize(L.Size());
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
				 C[i][j] = 1; 
				 (*s_index)[index] = 1;
			 }
			 else {
				 C[i][j] = 0; 
				 (*s_index)[index] = 0;
			 }

			 if ((i == 0 || i == N - 1) && C[i][j] == 1)
			 {
				 C[i][j] = 2;
				 (*s_index)[index] = 2;
			 }
			 if ((j == 0 || j == N - 1) && C[i][j] == 1)
			 {
				 C[i][j] = 2;
				 (*s_index)[index] = 2;
			 }

			 if( j > 0  )
			 {
				 if (C[i][j] == 1 && C[i][j - 1] == 0) 
				 { 
					 C[i][j] = 2; 
					 (*s_index)[index] = 2;
				 }
				 if (C[i][j] == 0 && C[i][j - 1] == 1) {
					 C[i][j - 1] = 2; 
					 (*s_index)[index - 1] = 2;
				 }
			
			 }

			 
			 index++;
			
			
		 }
	 }

	 for (size_t i = 0; i < (*s_index).size(); i++)
	 {
		 if((*s_index)[i] == 0 || (*s_index)[i] == 2)
		 {
			 for (size_t j = 0; j < (*s_index).size(); j++)
			 {
				 if (i != j) L[i][j] = 0;
				 else L[i][j] = 1;
			 }
		 }
	 }
	// cout << C;
	 return L;
 }

 void print_index(vector<int> index, int N) {

	 for (size_t i = 0; i < N; i++)
	 {
		 for (size_t j = 0; j < N; j++)
		 {
			 cout << index[i*N + j] << " ";
		 }
		 cout << endl;
	 }

 }

void CreateRightPart(Vector<double>* b,vector<int>* index, int N, double hx,double hy, double Rx, double Ry)
{
	

	int i = 0;
	int j = 0;
	for (size_t s = 0; s < b->Size(); s++)
	{
		//проверка дл€ j 
		if (i == N) i = 0;
		// найдем теперь i
		j = (int)(((double)(s - i)) / ((double)N));

		if ((*index)[s] == 2) (*b)[s] = g(i*hx - Rx, j*hy - Ry);
		if ((*index)[s] == 1) (*b)[s] = f(i*hx - Rx, j*hy - Ry);
	}

}