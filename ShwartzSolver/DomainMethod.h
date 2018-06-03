#pragma once
#include <vector>
#include <iostream>
#include <string>

#include "LaplasSolver.h"
#include "Relax.h"
using namespace std;
# define M_PI           3.14159265358979323846  /* pi */


double GetStep(int GridSize, double a_border, double b_border)
{

	if (a_border > b_border)
	{
		double temp = a_border;
		a_border = b_border;
		b_border = temp;
	}

	return (b_border - a_border) / (GridSize - 1);
}

double B(double x, double xi, double h)
{
	if (abs(x - xi) > h) return 0;
	double temp = (x - xi) / h;

	if (x - xi <= 0) return temp + 1;
	else return -temp + 1;
	
}

double F(double x, double y)
{
	return 0;
}

double G(double x, double y)
{
	return 1.0;
}

double Fpolar(double r, double phi)
{
	return F(r*cos(phi),r*sin(phi));
}


double Gpolar(double r, double phi)
{
	return G(r*cos(phi), r*sin(phi));
}

struct IndexPair {

	int index_i;
	int index_j;

};


class BoxDomain 
{
public:
	Vector<double> u;
	Vector<double> b;
	vector<IndexPair> border_in;
	vector<IndexPair> border_out;
	vector<double> out_x;
	vector<double> out_y;

	Matrix<double> L;

	double x_start;
	double y_start;

	double x_end;
	double y_end;

	int SIZE;
	double h;
	string name;
	

	BoxDomain(	string b_name,int N,double X_start, double Y_start,
				double X_end, double Y_end, double Radius)
	{
		name = b_name;
		SIZE = N;
		x_start = X_start;
		y_start = Y_start;
		x_end = X_end;
		y_end = Y_end;
		h = GetStep(SIZE, X_start, X_end);
		double hy = GetStep(SIZE, Y_start, Y_end);
		if (abs(h - hy) > 0.001) cout << "fail of box borders " + name << endl;



		FillRightPart();
		
		SetIndex(Radius);
		FillBorderCoords();
		vector<double> bord(border_in.size());
		UpdateBorder(bord.data(), bord.size());
		cout << name + " is generated" << endl;
	}

	void FillRightPart()
	{
		b.Resize(SIZE*SIZE);
		int s = 0;
		for (size_t i = 0; i < SIZE; i++)
		{
			for (size_t j = 0; j < SIZE; j++)
			{
				s = i*SIZE + j;
				
				if (i == 0 || j == 0 || j == SIZE - 1 || i == SIZE - 1)
					b[s] = G(x_start + j*h, y_start + i*h);
				else
					b[s] = F(x_start + j*h,y_start + i*h);
				
			}
		}
	
	}

	//заполняем два массива, отвечающих за границу области.
	// один за индексы внутри круга, другой снаружи 
	void SetIndex(double R)
	{
		//i == 0; -> y = y_start;
		double x_temp = 0;
		double y_temp = 0;
		double mod;
		IndexPair index;
		for (size_t j = 0; j < SIZE; j++)
		{
			x_temp = x_start + j*h;
			mod = sqrt( y_start*y_start + x_temp*x_temp );
			index.index_i = 0;
			index.index_j = j;
			if (mod <= R) border_in.push_back(index);
			else border_out.push_back(index);
		}

		//i == SIZE - 1 -> y = y_end
		for (size_t j = 0; j < SIZE; j++)
		{
			x_temp = x_start + j*h;
			mod = sqrt(y_end*y_end + x_temp*x_temp);
			index.index_i = SIZE - 1;
			index.index_j = j;
			if (mod <= R) border_in.push_back(index);
			else border_out.push_back(index);
		}

		//j == 0 - > x = x_start
		for (size_t i = 0; i < SIZE; i++)
		{
			y_temp = y_start + i*h;
			mod = sqrt(y_temp*y_temp + x_start*x_start);
			index.index_i = i;
			index.index_j = 0;
			if (mod <= R) border_in.push_back(index);
			else border_out.push_back(index);
		
		}

		//j == SIZE - 1 -> x = x_end
		for (size_t i = 0; i < SIZE; i++)
		{
			y_temp = y_start + i*h;
			mod = sqrt(y_temp*y_temp + x_end*x_end);
			index.index_i = i;
			index.index_j = SIZE - 1;
			if (mod <= R) border_in.push_back(index);
			else border_out.push_back(index);

		}
	}

	//заполняем массив х и у - точек на границе
	void FillBorderCoords()
	{
		out_x.resize(border_out.size());
		out_y.resize(border_out.size());
		for (size_t i = 0; i < out_x.size(); i++)
		{
			out_x[i] = x_start + h*border_out[i].index_j;
			out_y[i] = y_start + h*border_out[i].index_i;
		}
	}

	//получаем массив значений, которые устанавливаем на границе(заносим в правую часть)
	void UpdateBorder(double* border,int size)
	{
		int s = 0;
		int i = 0;
		int j = 0;
		for (size_t k = 0; k < size; k++)
		{
			i = border_in[k].index_i;
			j = border_in[k].index_j;
			s = i*SIZE + j;
			b[s] = border[k];
		}
		cout << name + " - border updated" << endl;
	}

	//определяем индекс i*h <= val <= (i+1)*h
	int FindIndex(double val, double val_start)
	{
		return (int)floor((val - val_start) / h);
	}

	//интерполируем значения в любой точке внутри сетки
	void GetInterpolatedPoints(double* x, double* y, double* g, int size)
	{
		int i = 0;
		int j = 0;
		int s1, s2, s3, s4;
		double S = 0;
		for (size_t k = 0; k < size; k++)
		{
			i = FindIndex(y[k], y_start);
			j = FindIndex(x[k], x_start);

			s1 = i*SIZE + j;
			s2 = (i + 1)*SIZE + j;
			s3 = i*SIZE + j + 1;
			s4 = (i + 1)*SIZE + j + 1;

			g[k] = u[s1] * B(x[k], x_start + j*h, h)*B(y[k], y_start + i*h, h) +
				u[s2] * B(x[k], x_start + j*h, h)*B(y[k], y_start + (i + 1)*h, h) +
				u[s3] * B(x[k], x_start + (j + 1)*h, h)*B(y[k], y_start + i*h, h) +
				u[s4] * B(x[k], x_start + (j + 1)*h, h)*B(y[k], y_start + (i + 1)*h, h);

		}
		cout << name + " - interpolating done";
	}

	void CorrectLaplasMatrix()
	{
		int s = 0;

		for (size_t k = 0; k < border_in.size(); k++)
		{
			s = border_in[k].index_i*SIZE + border_in[k].index_j;

			for (size_t i = 0; i < SIZE*SIZE; i++)
			{
				L[s][i] = 0;
			}
			L[s][s] = 1;
		}

		for (size_t k = 0; k < border_out.size(); k++)
		{
			s = border_out[k].index_i*SIZE + border_out[k].index_j;

			for (size_t i = 0; i < SIZE*SIZE; i++)
			{
				L[s][i] = 0;
			}
			L[s][s] = 1;
		}

	}

	void SolveLaplas(double w, double EPS)
	{
		if (L.Size() != SIZE*SIZE) 
		{
			L.Resize(SIZE*SIZE);
			Laplas(&L, SIZE, SIZE, h, h);
			u.Resize(SIZE*SIZE);
			CorrectLaplasMatrix();
		}
		Vector<double> x0 = u;
		RelaxFast(&L, &b, &x0, &u, w, EPS,SIZE);
		cout << "new solve complete for " + name << endl;
	}


	void Save()
	{

		ofstream myfile;
		myfile.open(name + ".txt");
		myfile << u.Size() << "\n";

		double x = 0;
		double y = 0;
		int s = 0;
		for (size_t i = 0; i < SIZE; i++)
		{
			for (size_t j = 0; j < SIZE; j++)
			{
				s = i*SIZE + j;
				x = x_start + j*h;
				y = y_start + i*h;

				myfile << x << " " << y << " " << u[s] << "\n";
			}
		}

		myfile.close();
		cout << "save to " + name + ".txt" << endl;


	}
};



struct Sector 
{
	int begin;
	int end;

	vector<double> x;
	vector<double> y;
};

class CircleDomain
{
public:
	Vector<double> u;
	Vector<double> b;
	
	vector<double> out_x;
	vector<double> out_y;
	vector<int> index;
	Sector leftSec;
	Sector topSec;
	Sector rigthSec;
	Sector bottomSec;

	Matrix<double> L;

	double RAD;
	int SIZE;
	double hr;
	double hphi;
	string name;


	CircleDomain(string b_name, int N,vector<double> x, vector<double> y, double Radius, double A0)
	{
		name = b_name;
		SIZE = N;
		RAD = Radius;
		hr = GetStep(SIZE, 0.0, RAD);
		hphi = GetStep(SIZE, 0.0, 2*M_PI);
		
		rigthSec.end = 0;
		if(x.size() == 8)
		{
			double phi = acos(x[0] / RAD);
			cout << 2.0*M_PI - phi << endl;
			cout << phi << endl;
			rigthSec.end = FindIndex(phi, hphi);
			rigthSec.begin = FindIndex(2.0*M_PI - phi, hphi);


			phi = acos(x[2] / RAD);
			cout << phi << endl;
			topSec.begin = FindIndex(phi, hphi);
			
			cout << "pi/2 = " << M_PI / 2 << endl;
			
			phi =  acos(x[3] / RAD);
			cout << phi << endl;
			topSec.end = FindIndex(phi, hphi);



			phi = acos(x[4] / RAD);
			cout << phi << endl;
			leftSec.begin = FindIndex(phi, hphi);
			
			cout << "pi = " << M_PI << endl;

			phi = 2*M_PI - acos(x[5] / RAD);
			cout << phi << endl;
			leftSec.end = FindIndex(phi, hphi);


			phi = 2 * M_PI - acos(x[6] / RAD);
			cout << phi << endl;
			bottomSec.begin = FindIndex(phi, hphi);
			
			cout << "3pi/2 = " << 3*M_PI/2 << endl;

			phi = 2 * M_PI - acos(x[7] / RAD);
			cout << phi << endl;
			bottomSec.end = FindIndex(phi, hphi);


			cout << rigthSec.begin << "  " << rigthSec.end << endl;
			cout << topSec.begin << "  " << topSec.end << endl;
			cout << leftSec.begin << "  " << leftSec.end << endl;
			cout << bottomSec.begin << "  " << bottomSec.end << endl;

		}

		FillBorderCoords();
		FillRightPart(A0);

		cout << name + " is generated" << endl;
		
	}

	void FillBorderCoords() 
	{
		int N = topSec.end - topSec.begin;
		topSec.x.resize(N);
		topSec.y.resize(N);
		for (size_t i = 0; i < N; i++)
		{
			topSec.x[i] = RAD*cos((topSec.begin + i)*hphi);
			topSec.y[i] = RAD*sin((topSec.begin + i)*hphi);
		}


		N = leftSec.end - leftSec.begin;
		leftSec.x.resize(N);
		leftSec.y.resize(N);
		for (size_t i = 0; i < N; i++)
		{
			leftSec.x[i] = RAD*cos((leftSec.begin + i)*hphi);
			leftSec.y[i] = RAD*sin((leftSec.begin + i)*hphi);
		}

		N = bottomSec.end - bottomSec.begin;
		bottomSec.x.resize(N);
		bottomSec.y.resize(N);
		for (size_t i = 0; i < N; i++)
		{
			bottomSec.x[i] = RAD*cos((bottomSec.begin + i)*hphi);
			bottomSec.y[i] = RAD*sin((bottomSec.begin + i)*hphi);
		}

		N = rigthSec.end + SIZE - rigthSec.begin;
		rigthSec.x.resize(N);
		rigthSec.y.resize(N);
		size_t i = 0;
		for ( i = 0; rigthSec.begin + i < SIZE; i++)
		{
			rigthSec.x[i] = RAD*cos((rigthSec.begin + i)*hphi);
			rigthSec.y[i] = RAD*sin((rigthSec.begin + i)*hphi);
		}
		int k = i;
		for (;  i < N; i++)
		{
			rigthSec.x[i] = RAD*cos((i - k)*hphi);
			rigthSec.y[i] = RAD*sin((i - k)*hphi);
		}
	}

	bool IsInside(int k)
	{
		if (k >= topSec.begin && k <= topSec.end)
		{
			cout << "1" << endl;
			return true;
		}
		if (k >= leftSec.begin && k <= leftSec.end)
		{
			cout << "2" << endl;
			return true;
		}
		if (k >= bottomSec.begin && k <= bottomSec.end)
		{
			cout << "3" << endl;
			return true;
		}
		if (k >= rigthSec.begin && k < SIZE)
		{
			cout << "4" << endl;
			return true;
		}
		if (k >= 0 && k <= rigthSec.end)
		{
			cout << "5" << endl;
			return true;
		}


			return false;
	}

	void FillRightPart(double A0)
	{
		b.Resize(SIZE*SIZE);
		int s = 0;
		for (size_t i = 0; i < SIZE; i++)
		{
			for (size_t j = 0; j < SIZE; j++)
			{
				s = i*SIZE + j;
				b[s] = Fpolar(j*hr, i*hphi);

					if (j == SIZE - 1) 
					{
						//if (!IsInside(i))
							b[s] = Gpolar(RAD, i*hphi);
							//index.push_back(s);
						//else
						//	b[s] = 0;
					}
					if (j == 0)
					{
						//b[s] = A0;
						//index.push_back(s);
					}
						

					if(i == 0)
					{
						b[s] = 0;
					}
				
				
					

			}
		}
		//cout << b << endl;
	}

	
	//получаем массив значений, которые устанавливаем на границе(заносим в правую часть)
	void UpdateBorder(double* right_border, int right_size,
						double* top_border, int top_size,
						double* left_border, int left_size,
						double* bottom_border, int bottom_size)
	{
		
		int s = 0;
		int j = SIZE - 1;
		for (size_t i = 0; i < top_size; i++)
		{
			s = (topSec.begin + i)*SIZE + j;
			b[s] = top_border[i];
		}



		for (size_t i = 0; i < left_size; i++)
		{
			s = (leftSec.begin + i)*SIZE + j;
			b[s] = left_border[i];
		}


		for (size_t i = 0; i < bottom_size; i++)
		{
			s = (bottomSec.begin + i)*SIZE + j;
			b[s] = bottom_border[i];
		}

		for (size_t i = 0; i < SIZE - rigthSec.begin; i++)
		{
			s = (rigthSec.begin + i)*SIZE + j;
			b[s] = right_border[i];
		}

		int k = SIZE - rigthSec.begin;
		for (size_t i = 0; i < rigthSec.end; i++)
		{
			s = ( i)*SIZE + j;
			b[s] = right_border[i + k];
		}


		cout << name + " - border updated" << endl;
	}

	//определяем индекс i*h <= val <= (i+1)*h
	int FindIndex(double val, double h)
	{
		return (int)floor((val ) / h);
	}

	void GetInterpolatedPointsXY(double* x, double* y, double* g, int size) 
	{
		vector<double> r(size);
		vector<double> phi(size);

		for (size_t i = 0; i < size; i++)
		{
			r[i] = sqrt(x[i] * x[i] + y[i] * y[i]);
			if(y[i] >= 0)
			phi[i] = acos(x[i] / r[i]);
			else
				phi[i] = 2*M_PI - acos(x[i] / r[i]);
		}
		GetInterpolatedPoints(r.data(), phi.data(), g, size);
	}

	//интерполируем значения в любой точке внутри сетки
	void GetInterpolatedPoints(double* r, double* phi, double* g, int size)
	{
		int i = 0;
		int j = 0;
		int s1, s2, s3, s4;
		double S = 0;
		for (size_t k = 0; k < size; k++)
		{
			i = FindIndex(phi[k], hphi);
			j = FindIndex(r[k], hr);

			s1 = i*SIZE + j;
			s2 = (i + 1)*SIZE + j;
			s3 = i*SIZE + j + 1;
			s4 = (i + 1)*SIZE + j + 1;

			g[k] = u[s1] * B(r[k], j*hr, hr)*B(phi[k], i*hphi, hphi) +
				u[s2] * B(r[k], j*hr, hr)*B(phi[k],(i + 1)*hphi, hphi) +
				u[s3] * B(r[k], (j + 1)*hr, hr)*B(phi[k], i*hphi, hphi) +
				u[s4] * B(r[k], (j + 1)*hr, hr)*B(phi[k],(i + 1)*hphi, hphi);

		}
		cout << name + " - interpolating done";
	}

	void CorrectLaplasPolarMatrix()
	{
/*
		A[s][s] = -2.0 * (1.0 / (h_r * h_r) + 1.0 / (ri * ri * h_phi * h_phi));

		A[s][s + 1] = 1.0 / (h_r * h_r) - 1.0 / (2.0 * ri * h_r * h_r);
		A[s][s - 1] = 1.0 / (h_r * h_r) + 1.0 / (2.0 * ri * h_r * h_r);

		A[s][GRID_DIM_N + j] = 1.0 / (ri * ri * h_phi * h_phi);
		A[s][s - GRID_DIM_N] = 1.0 / (ri * ri * h_phi * h_phi);*/
		int s = 0;
		
		//phi == 2pi
		double ri = 0;
		double h_r = hr;
		double h_phi = hphi;

		int i = SIZE - 1;
		for (size_t j = 1; j < SIZE - 1; j++)
		{
			s = i*SIZE + j;
			
			L[s][s] = -2.0*(1.0 / (hr*hr) + 1.0 / (j*j*hr*hr*hphi*hphi));
			L[s][s + 1] = (1.0 / (hr*hr) - 1.0 / (2.0*j*hr*hr));
			L[s][s - 1] = (1.0 / (hr*hr) + 1.0 / (2.0*j*hr*hr));
			L[s][ SIZE + j ] = 1.0 / (j*j*hr*hr*hphi*hphi);
			L[s][s - SIZE] = 1.0 / (j*j*hr*hr*hphi*hphi);
			
		}
		i = SIZE - 2;
		for (size_t j = 1; j < SIZE - 1; j++)
		{
			s = i*SIZE + j;

			L[s][s] = -2.0*(1.0 / (hr*hr) + 1.0 / (j*j*hr*hr*hphi*hphi));
			L[s][s + 1] = (1.0 / (hr*hr) - 1.0 / (2.0*j*hr*hr));
			L[s][s - 1] = (1.0 / (hr*hr) + 1.0 / (2.0*j*hr*hr));
			L[s][ j] = 1.0 / (j*j*hr*hr*hphi*hphi);
			L[s][s - SIZE] = 1.0 / (j*j*hr*hr*hphi*hphi);

		}

		//phi == 0
		i = SIZE - 1;
		for (size_t j = 1; j < SIZE - 1; j++)
		{
			s = i*SIZE + j;

			L[j][j] = -2.0*(1.0 / (hr*hr) + 1.0 / (j*j*hr*hr*hphi*hphi));
			L[j][j + 1] = (1.0 / (hr*hr) - 1.0 / (2.0*j*hr*hr));
			L[j][j - 1] = (1.0 / (hr*hr) + 1.0 / (2.0*j*hr*hr));
			L[j][j + SIZE] = 1.0 / (j*j*hr*hr*hphi*hphi);
			L[j][s - SIZE] = 1.0 / (j*j*hr*hr*hphi*hphi);

		}
		i = 1;
		for (size_t j = 1; j < SIZE - 1; j++)
		{
			s = i*SIZE + j;
			
			L[s][s] = -2.0*(1.0 / (hr*hr) + 1.0 / (j*j*hr*hr*hphi*hphi));
			L[s][s + 1] = (1.0 / (hr*hr) - 1.0 / (2.0*j*hr*hr));
			L[s][s - 1] = (1.0 / (hr*hr) + 1.0 / (2.0*j*hr*hr));
			L[s][s + SIZE] = 1.0 / (j*j*hr*hr*hphi*hphi);
			L[s][ SIZE*(SIZE-1) + j ] = 1.0 / (j*j*hr*hr*hphi*hphi);

		}
		//for (size_t j = 1; j < SIZE - 1; j++)
		//{
		//	s = i*SIZE + j;

		//	L[j][j] = 1;
		//	//L[j][s] = -1;
		//	b[j] = 1;
		//}

		//r == 0
		////r == R

		//for (size_t k = 0; k < index.size(); k++)
		//{
		//	for (size_t l = 0; l < SIZE*SIZE; l++)
		//	{
		//		L[index[k]][l] = 0;
		//	}
		//	L[index[k]][index[k]] = 1;
		//	cout << "L[s][s] = " << L[index[k]][index[k]] << "; b[s] = " << b[index[k]] << endl;
		//	//	
		//}
		for (i = 0; i < SIZE ; i++)
		{
			s = i*SIZE;

			for (size_t SS = 0; SS < SIZE*SIZE; SS++)
			{
				L[s][SS] = 0;
				L[s + SIZE - 1][SS] = 0;
			}

			L[s][s] = 1;
			L[s + SIZE - 1][s + SIZE - 1] = 1;
			b[s] = 1;
			b[s + SIZE - 1] = 1;

		
			cout << "L[s][s] = " << L[s][s] << "; b[s] = " << b[s] << endl;
			cout << "L[s + ][s + ] = " << L[s + SIZE - 1][s + SIZE - 1] << "; b[s +] = " << b[s + SIZE - 1] << endl;
		}

/*
		for (size_t k = 1; k < SIZE*SIZE - 1; k++)
		{
			cout << L[k][k - 1] << " " << L[k][k] << " " << L[k][k + 1] <<endl;
		}
*/
	}

	void SolveLaplas(double w, double EPS)
	{
		if (L.Size() != SIZE*SIZE)
		{
			L.Resize(SIZE*SIZE);
			LaplasPolar2(&L, SIZE, SIZE, hr, hphi);
			u.Resize(SIZE*SIZE);
			CorrectLaplasPolarMatrix();
			if(SIZE == 9)cout << L << endl;
		}
		Vector<double> x0 = u;
		RelaxFastBad(&L, &b, &x0, &u, w, EPS, SIZE);
		cout << "new solve complete for " + name << endl;
	}

	void Save()
	{

		ofstream myfile;
		myfile.open(name + ".txt");
		myfile << u.Size() << "\n";
		
		double x = 0;
		double y = 0;
		int s = 0;
		for (size_t i = 0; i < SIZE; i++)
		{
			for (size_t j = 0; j < SIZE; j++)
			{
				s = i*SIZE + j;
				x = j*hr*cos(i*hphi);
				y = j*hr*sin(i*hphi);

				myfile << x << " " << y << " " << u[s] << "\n";
			}
		}

		myfile.close();
		cout << "save to " + name + ".txt" << endl;


	}
};

