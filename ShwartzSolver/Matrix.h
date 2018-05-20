#pragma once
#include "Vector.h"
template <class T>
class Matrix
{
	size_t size;
	T **m;
public:
	Matrix()
	{
		size = 4;
		m = new T*[size];
		for (int i = 0; i < size; ++i)
			m[i] = new T[size];
		for (int i = 0; i < size; ++i)
		{
			for (int j = 0; j < size; ++j)
			{
				m[i][j] = 0;
			}
		}
	}
	Matrix(size_t n)
	{
		size = n;
		m = new T*[size];
		for (int i = 0; i < size; ++i)
			m[i] = new T[size];
		for (int i = 0; i < size; ++i)
		{
			for (int j = 0; j < size; ++j)
			{
				m[i][j] = 0;
			}
		}
	}
	Matrix(const Matrix &rV)
	{
		this->size = rV.size;
		m = new T*[size];
		for (int i = 0; i < size; ++i)
			m[i] = new T[size];
		for (int i = 0; i < size; ++i)
		{
			for (int j = 0; j < size; ++j)
			{
				m[i][j] = rV.m[i][j];
			}
		}
	}
	Matrix operator=(const Matrix &rV)
	{
		this->size = rV.size;
		for (int i = 0; i < size; ++i)
		{
			for (int j = 0; j < size; ++j)
			{
				m[i][j] = rV[i][j];
			}
		}
		return *this;
	}

	Matrix operator + (Matrix &a)
	{
		Matrix c(size);

		for (size_t i = 0; i < size; i++)
		{
			//c[i] = v[i] + a[i];
			for (size_t j = 0; j < size; j++)
			{
				c[i][j] = v[i][j] + a[i][j];
			}
		}

		return c;
	}

	Matrix operator - (Matrix &a)
	{
		Matrix c(size);

		for (size_t i = 0; i < size; i++)
		{
			//c[i] = v[i] + a[i];
			for (size_t j = 0; j < size; j++)
			{
				c[i][j] = v[i][j] + a[i][j];
			}
		}

		return c;
	}

	Matrix operator * (Matrix &B)
	{
		Matrix C = new Matrix(B.size);

		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				for (int r = 0; r < size; r++)
					C[i][j] += m[i][r] * B[r][j];

		return C;
	}

	Vector<T> operator * (Vector<T> &x)
	{
		Vector<T> b(x.Size());

		
		//A*x = b
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				b[i] += m[i][j] * x[j];
			}
		}
		return b;
	}

	friend istream &operator >> (istream &stream, Matrix &a)
	{
		for (int i = 0; i < a.size; ++i)
		{
			for (int j = 0; j < a.size; ++j)
				stream >> a.m[i][j];
		}
		return stream;
	}
	friend ostream &operator <<(ostream &stream, Matrix &a)
	{
		for (int i = 0; i < a.size; ++i)
		{
			for (int j = 0; j < a.size; ++j)
				stream << a.m[i][j] << " ";
			stream << endl;
		}
		return stream;
	}
	T* operator [](int k)
	{
		if (k < size)
			return m[k];
		else
			return NULL;
	}
	~Matrix()
	{
		for (int i = 0; i < size; ++i)
		{
			delete[] m[i];
		}
		delete[]m;
	}

	size_t Size()
	{
		return size;
	}

	
};