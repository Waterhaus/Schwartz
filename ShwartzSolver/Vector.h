#pragma once

#include <vector>

using namespace std;


template <class T>
class Vector
{
	size_t size;
	vector<T> v;
public:
	Vector()
	{
		size = 4;
		v.resize(size);
		for (int i = 0; i < size; ++i)
			v[i] = 0;
	}

	Vector(size_t n)
	{
		size = n;
		v.resize(size);
		for (int i = 0; i < size; ++i)
			v[i] = 0;
	}

	Vector(const  Vector &rV)
	{
		this->size = rV.size;
		v.resize(size);
		for (int i = 0; i < this->size; ++i)
			this->v[i] = rV.v[i];
	}

	Vector(T *mas, size_t n)
	{
		size = n;
		v.resize(size);
		for (int i = 0; i < size; ++i)
			v[i] = mas[i];
	}


	Vector operator = (const  Vector &rV)
	{
		this->size = rV.size;
		for (int i = 0; i < this->size; ++i)
			this->v[i] = rV.v[i];
		return *this;
	}

	size_t Size()
	{
		return size;
	}

	T &operator[] (int i)
	{
		return v[i];
	}

	~Vector()
	{
		
	}


	friend ostream &operator << (ostream &stream, Vector &a)
	{
		for (int i = 0; i < a.size; ++i) stream << a.v[i] << " ";
		stream << endl;
		return stream;
	}
	friend istream &operator >> (istream &stream, Vector &a)
	{
		for (int i = 0; i < a.size; ++i) stream >> a.v[i];
		return stream;
	}


	T norm2()
	{
		T S = 0.0;
		for (size_t i = 0; i < size; i++)
		{
			S = S + v[i] * v[i];
		}
		return sqrt(S);
	}

	T normInf()
	{
		T max = abs(v[0]);
		T temp = 0;
		for (size_t i = 1; i < size; i++)
		{
			temp = abs(v[i]);
			if (temp > max) max = temp;
		}

		return max;
	}




	static T ScalarProduct(Vector a, Vector b)
	{
		size_t size = a.Size();
		T S = 0.0;
		for (size_t i = 0; i < size; i++)
		{
			S = S + a[i] * b[i];
		}
		return S;

	}

	static Vector<T> Identity(int N)
	{
		size_t size = N;
		Vector<T> vec(size);
		T S = 0.0;
		for (size_t i = 0; i < size; i++)
		{
			vec[i] = 1;
		}
		return vec;

	}

	Vector operator + (Vector &a)
	{
		Vector c(size);

		for (size_t i = 0; i < size; i++)
		{
			c[i] = v[i] + a[i];
		}

		return c;
	}

	Vector operator - (Vector &a)
	{
		Vector c(size);

		for (size_t i = 0; i < size; i++)
		{
			c[i] = v[i] - a[i];
		}

		return c;
	}



	Vector operator * (T alpha)
	{
		Vector b(size);

		for (size_t i = 0; i < size; i++)
		{
			b[i] = v[i] * alpha;
		}

		return b;
	}

	void Resize(int SIZE)
	{
		v.resize(SIZE);
		size = SIZE;
	}

};
