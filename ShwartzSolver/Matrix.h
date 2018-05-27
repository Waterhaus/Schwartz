#pragma once
#include "Vector.h"
template <class T>
class Matrix
{
	size_t size_i, size_j;
	
	vector<vector<T>> Mat;

public:
	Matrix()
	{
		size_i = size_j = 4;
		Mat.resize(size_i);

		for (int i = 0; i < size_i; ++i)
			Mat[i].resize(size_j);

		for (int i = 0; i < size_i; ++i)
		{
			for (int j = 0; j < size_j; ++j)
			{
				Mat[i][j] = 0;
			}
		}
	}

	Matrix(size_t n)
	{
		size_i = size_j = n;
		Mat.resize(size_i);

		for (int i = 0; i < size_i; ++i)
			Mat[i].resize(size_j);

		for (int i = 0; i < size_i; ++i)
		{
			for (int j = 0; j < size_i; ++j)
			{
				Mat[i][j] = 0;
			}
		}
	}

	Matrix(size_t N, size_t M)
	{
		size_i = N;
		size_j = M;
		Mat.resize(N);
		for (int i = 0; i < N; ++i)
			Mat[i].resize(M);

		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < M; ++j)
			{
				Mat[i][j] = 0;
			}
		}
	}
	Matrix(const Matrix &rV)
	{
		this->size_i = rV.size_i;
		this->size_j = rV.size_j;
		Mat.resize(size_i);
		for (int i = 0; i < size_i; ++i)
			Mat[i].resize(size_j);


		for (int i = 0; i < size_i; ++i)
		{
			for (int j = 0; j < size_j; ++j)
			{
				Mat[i][j] = rV.Mat[i][j];
			}
		}
	}


	Matrix operator=(const Matrix &rV)
	{
		this->size_i = rV.size_i;
		this->size_j = rV.size_j;

		Mat.resize(size_i);
		for (int i = 0; i < size_i; ++i)
			Mat[i].resize(size_j);

		for (int i = 0; i < size_i; ++i)
		{
			for (int j = 0; j < size_i; ++j)
			{
				Mat[i][j] = rV[i][j];
			}
		}
		return *this;
	}

	Matrix operator + (Matrix &a)
	{
		Matrix c(size_i);

		for (size_t i = 0; i < size_i; i++)
		{
			//c[i] = v[i] + a[i];
			for (size_t j = 0; j < size_i; j++)
			{
				c[i][j] = v[i][j] + a[i][j];
			}
		}

		return c;
	}

	Matrix operator - (Matrix &a)
	{
		Matrix c(size_i);

		for (size_t i = 0; i < size_i; i++)
		{
			//c[i] = v[i] + a[i];
			for (size_t j = 0; j < size_i; j++)
			{
				c[i][j] = v[i][j] + a[i][j];
			}
		}

		return c;
	}

	Matrix operator * (Matrix &B)
	{
		Matrix C = new Matrix(B.size_i);

		for (int i = 0; i < size_i; i++)
			for (int j = 0; j < size_i; j++)
				for (int r = 0; r < size_i; r++)
					C[i][j] += Mat[i][r] * B[r][j];

		return C;
	}

	Vector<T> operator * (Vector<T> &x)
	{
		Vector<T> b(x.Size());

		
		//A*x = b
		for (int i = 0; i < size_i; i++)
		{
			for (int j = 0; j < size_i; j++)
			{
				b[i] += Mat[i][j] * x[j];
			}
		}
		return b;
	}

	friend istream &operator >> (istream &stream, Matrix &a)
	{
		for (int i = 0; i < a.size_i; ++i)
		{
			for (int j = 0; j < a.size_j; ++j)
				stream >> a.Mat[i][j];
		}
		return stream;
	}
	friend ostream &operator <<(ostream &stream, Matrix &a)
	{
		for (int i = 0; i < a.size_i; ++i)
		{
			for (int j = 0; j < a.size_j; ++j)
				stream << a.Mat[i][j] << " ";
			stream << endl;
		}
		return stream;
	}
	vector<T> &operator [](int k)
	{
		if (k < size_i)
			return Mat[k];
			
	}
	~Matrix()
	{
		 
	}

	size_t Size()
	{
		return size_i;
	}

	size_t SizeN()
	{
		return size_i;
	}
	size_t SizeM()
	{
		return size_j;
	}
	
};