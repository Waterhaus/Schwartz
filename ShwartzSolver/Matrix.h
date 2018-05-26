#pragma once
#include "Vector.h"
template <class T>
class Matrix
{
	size_t size_i, size_j;
	T **m;
public:
	Matrix()
	{
		size_i = size_j = 4;
		m = new T*[size_i];
		for (int i = 0; i < size_i; ++i)
			m[i] = new T[size_j];
		for (int i = 0; i < size_i; ++i)
		{
			for (int j = 0; j < size_j; ++j)
			{
				m[i][j] = 0;
			}
		}
	}

	Matrix(size_t n)
	{
		size_i = size_j = n;
		m = new T*[size_i];
		for (int i = 0; i < size_i; ++i)
			m[i] = new T[size_i];
		for (int i = 0; i < size_i; ++i)
		{
			for (int j = 0; j < size_i; ++j)
			{
				m[i][j] = 0;
			}
		}
	}

	Matrix(size_t N, size_t M)
	{
		size_i = N;
		size_j = M;
		m = new T*[N];
		for (int i = 0; i < N; ++i)
			m[i] = new T[M];
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < M; ++j)
			{
				m[i][j] = 0;
			}
		}
	}
	Matrix(const Matrix &rV)
	{
		this->size_i = rV.size_i;
		this->size_j = rV.size_j;
		m = new T*[size_i];
		for (int i = 0; i < size_i; ++i)
			m[i] = new T[size_j];
		for (int i = 0; i < size_i; ++i)
		{
			for (int j = 0; j < size_j; ++j)
			{
				m[i][j] = rV.m[i][j];
			}
		}
	}
	Matrix operator=(const Matrix &rV)
	{
		this->size_i = rV.size_i;
		for (int i = 0; i < size_i; ++i)
		{
			for (int j = 0; j < size_i; ++j)
			{
				m[i][j] = rV[i][j];
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
					C[i][j] += m[i][r] * B[r][j];

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
				b[i] += m[i][j] * x[j];
			}
		}
		return b;
	}

	friend istream &operator >> (istream &stream, Matrix &a)
	{
		for (int i = 0; i < a.size_i; ++i)
		{
			for (int j = 0; j < a.size_j; ++j)
				stream >> a.m[i][j];
		}
		return stream;
	}
	friend ostream &operator <<(ostream &stream, Matrix &a)
	{
		for (int i = 0; i < a.size_i; ++i)
		{
			for (int j = 0; j < a.size_j; ++j)
				stream << a.m[i][j] << " ";
			stream << endl;
		}
		return stream;
	}
	T* operator [](int k)
	{
		if (k < size_i)
			return m[k];
		else
			return NULL;
	}
	~Matrix()
	{
		for (int i = 0; i < size_i; ++i)
		{
			delete[] m[i];
		}
		delete[]m;
	}

	size_t Size()
	{
		return size_i;
	}

	
};