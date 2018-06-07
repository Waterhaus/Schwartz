#pragma once
#include <vector>
#include <iostream>
#include <string>

#include "LaplasSolver.h"
#include "Relax.h"
using namespace std;



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

	//Matrix<double> L;
	SparseMatrix Lap;

	double RAD;
	int SIZE;
	double hr;
	double hphi;
	string name;


	CircleDomain(string b_name, int N, vector<double> x, vector<double> y, double Radius, double A0)
	{
		name = b_name;
		SIZE = N;
		RAD = Radius;
		hr = GetStep(SIZE, 0.0, RAD);
		hphi = GetStep(SIZE, 0.0, 2 * M_PI);

		rigthSec.end = 0;
		if (x.size() == 8)
		{
			double phi = acos(x[0] / RAD);
			rigthSec.end = FindIndex(phi, hphi) + 1;
			rigthSec.begin = FindIndex(2.0*M_PI - phi, hphi);


			phi = acos(x[2] / RAD);
			topSec.begin = FindIndex(phi, hphi);


			phi = acos(x[3] / RAD);
			topSec.end = FindIndex(phi, hphi) ;



			phi = acos(x[4] / RAD);
			leftSec.begin = FindIndex(phi, hphi) ;


			phi = 2 * M_PI - acos(x[5] / RAD);
			leftSec.end = FindIndex(phi, hphi)  + 1;


			phi = 2 * M_PI - acos(x[6] / RAD);
			bottomSec.begin = FindIndex(phi, hphi);


			phi = 2 * M_PI - acos(x[7] / RAD);
			bottomSec.end = FindIndex(phi, hphi) + 1;


			cout << rigthSec.begin << "  " << rigthSec.end << endl;
			cout << topSec.begin << "  " << topSec.end << endl;
			cout << leftSec.begin << "  " << leftSec.end << endl;
			cout << bottomSec.begin << "  " << bottomSec.end << endl;

		}

		FillBorderCoords();
		FillRightPart(A0);
		TestSetBorder();
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
		for (i = 0; rigthSec.begin + i < SIZE; i++)
		{
			rigthSec.x[i] = RAD*cos((rigthSec.begin + i)*hphi);
			rigthSec.y[i] = RAD*sin((rigthSec.begin + i)*hphi);
		}
		int k = i;
		for (; i < N; i++)
		{
			rigthSec.x[i] = RAD*cos((i - k)*hphi);
			rigthSec.y[i] = RAD*sin((i - k)*hphi);
		}
	}

	bool IsInside(int k)
	{
		if (k >= topSec.begin && k < topSec.end)
		{
			
			return true;
		}
		if (k >= leftSec.begin && k < leftSec.end)
		{
			
			return true;
		}
		if (k >= bottomSec.begin && k < bottomSec.end)
		{
			
			return true;
		}
		if (k >= rigthSec.begin && k < SIZE)
		{
			
			return true;
		}
		if (k >= 0 && k < rigthSec.end)
		{
			
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
					if (!IsInside(i))
					b[s] = Gpolar(RAD, i*hphi);
					//index.push_back(s);
					//else
					//	b[s] = 0;
				}
				if (j == 0)
				{
					b[s] = A0;
					//index.push_back(s);
				}


				if (i == 0)
				{
					b[s] = 0;
				}




			}
		}
		//cout << b << endl;
	}


	
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
			s = (i)*SIZE + j;
			b[s] = right_border[i + k];
		}


		cout << name + " - border updated" << endl;
	}


	void TestSetBorder() {

		vector<double> temp1(rigthSec.x.size());
		for (size_t i = 0; i < temp1.size(); i++)
		{
			temp1[i] = 1;
		}
		

		vector<double> temp2(topSec.x.size());
		for (size_t i = 0; i < temp2.size(); i++)
		{
			temp2[i] = 1;
		}

		vector<double> temp3(leftSec.x.size());
		for (size_t i = 0; i < temp3.size(); i++)
		{
			temp3[i] = 1;
		}

		vector<double> temp4(bottomSec.x.size());
		for (size_t i = 0; i < temp4.size(); i++)
		{
			temp4[i] = 1;
		}
		UpdateBorder(temp1.data(), temp1.size(), temp2.data(), temp2.size(), temp3.data(), temp3.size(), temp4.data(), temp4.size());
	}
	//i*h <= val <= (i+1)*h
	int FindIndex(double val, double h)
	{
		return (int)floor((val) / h);
	}

	void GetInterpolatedPointsXY(double* x, double* y, double* g, int size)
	{
		vector<double> r(size);
		vector<double> phi(size);

		for (size_t i = 0; i < size; i++)
		{
			r[i] = sqrt(x[i] * x[i] + y[i] * y[i]);
			if (y[i] >= 0)
				phi[i] = acos(x[i] / r[i]);
			else
				phi[i] = 2 * M_PI - acos(x[i] / r[i]);
		}
		GetInterpolatedPoints(r.data(), phi.data(), g, size);
	}

	//
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
				u[s2] * B(r[k], j*hr, hr)*B(phi[k], (i + 1)*hphi, hphi) +
				u[s3] * B(r[k], (j + 1)*hr, hr)*B(phi[k], i*hphi, hphi) +
				u[s4] * B(r[k], (j + 1)*hr, hr)*B(phi[k], (i + 1)*hphi, hphi);

		}
		cout << name + " - interpolating done";
	}

	void CorrectLaplasPolarMatrix()
	{

		int s = 0;

		//phi == 2pi

		double temp;
		int i = SIZE - 1;
		for (size_t j = 1; j < SIZE - 1; j++)
		{
			s = i*SIZE + j;

			//L[s][s] = -2.0*(1.0 / (hr*hr) + 1.0 / (j*j*hr*hr*hphi*hphi));
			temp = -2.0*(1.0 / (hr*hr) + 1.0 / (j*j*hr*hr*hphi*hphi));
			Lap.Mat[s].push_back(temp);
			Lap.index[s].push_back(s);

			//L[s][s + 1] = (1.0 / (hr*hr) + 1.0 / (2.0*j*hr*hr));
			temp = (1.0 / (hr*hr) + 1.0 / (2.0*j*hr*hr));
			Lap.Mat[s].push_back(temp);
			Lap.index[s].push_back(s + 1);


			//L[s][s - 1] = (1.0 / (hr*hr) - 1.0 / (2.0*j*hr*hr));
			temp = (1.0 / (hr*hr) - 1.0 / (2.0*j*hr*hr));
			Lap.Mat[s].push_back(temp);
			Lap.index[s].push_back(s - 1);

			//L[s][ SIZE + j ] = 1.0 / (j*j*hr*hr*hphi*hphi);
			temp = 1.0 / (j*j*hr*hr*hphi*hphi);
			Lap.Mat[s].push_back(temp);
			Lap.index[s].push_back(SIZE + j);

			//L[s][s - SIZE] = 1.0 / (j*j*hr*hr*hphi*hphi);
			temp = 1.0 / (j*j*hr*hr*hphi*hphi);
			Lap.Mat[s].push_back(temp);
			Lap.index[s].push_back(s - SIZE);

		}
		i = SIZE - 2;
		for (size_t j = 1; j < SIZE - 1; j++)
		{
			s = i*SIZE + j;

			/*	L[s][s] = -2.0*(1.0 / (hr*hr) + 1.0 / (j*j*hr*hr*hphi*hphi));
			L[s][s + 1] = (1.0 / (hr*hr) + 1.0 / (2.0*j*hr*hr));
			L[s][s - 1] = (1.0 / (hr*hr) - 1.0 / (2.0*j*hr*hr));
			L[s][ j] = 1.0 / (j*j*hr*hr*hphi*hphi);
			L[s][s - SIZE] = 1.0 / (j*j*hr*hr*hphi*hphi);*/

			temp = -2.0*(1.0 / (hr*hr) + 1.0 / (j*j*hr*hr*hphi*hphi));
			Lap.Mat[s].push_back(temp);
			Lap.index[s].push_back(s);


			temp = (1.0 / (hr*hr) + 1.0 / (2.0*j*hr*hr));
			Lap.Mat[s].push_back(temp);
			Lap.index[s].push_back(s + 1);



			temp = (1.0 / (hr*hr) - 1.0 / (2.0*j*hr*hr));
			Lap.Mat[s].push_back(temp);
			Lap.index[s].push_back(s - 1);


			temp = 1.0 / (j*j*hr*hr*hphi*hphi);
			Lap.Mat[s].push_back(temp);
			Lap.index[s].push_back(j);


			temp = 1.0 / (j*j*hr*hr*hphi*hphi);
			Lap.Mat[s].push_back(temp);
			Lap.index[s].push_back(s - SIZE);

		}

		//phi == 0
		i = SIZE - 1;
		for (size_t j = 1; j < SIZE - 1; j++)
		{
			s = i*SIZE + j;

			/*L[j][j] = -2.0*(1.0 / (hr*hr) + 1.0 / (j*j*hr*hr*hphi*hphi));
			L[j][j + 1] = (1.0 / (hr*hr) + 1.0 / (2.0*j*hr*hr));
			L[j][j - 1] = (1.0 / (hr*hr) - 1.0 / (2.0*j*hr*hr));
			L[j][j + SIZE] = 1.0 / (j*j*hr*hr*hphi*hphi);
			L[j][s - SIZE] = 1.0 / (j*j*hr*hr*hphi*hphi);*/
			temp = -2.0*(1.0 / (hr*hr) + 1.0 / (j*j*hr*hr*hphi*hphi));
			Lap.Mat[j].push_back(temp);
			Lap.index[j].push_back(j);


			temp = (1.0 / (hr*hr) + 1.0 / (2.0*j*hr*hr));
			Lap.Mat[j].push_back(temp);
			Lap.index[j].push_back(j + 1);



			temp = (1.0 / (hr*hr) - 1.0 / (2.0*j*hr*hr));
			Lap.Mat[j].push_back(temp);
			Lap.index[j].push_back(j - 1);


			temp = 1.0 / (j*j*hr*hr*hphi*hphi);
			Lap.Mat[j].push_back(temp);
			Lap.index[j].push_back(j + SIZE);


			temp = 1.0 / (j*j*hr*hr*hphi*hphi);
			Lap.Mat[j].push_back(temp);
			Lap.index[j].push_back(s - SIZE);

		}
		i = 1;
		for (size_t j = 1; j < SIZE - 1; j++)
		{
			s = i*SIZE + j;

			/*	L[s][s] = -2.0*(1.0 / (hr*hr) + 1.0 / (j*j*hr*hr*hphi*hphi));
			L[s][s + 1] = (1.0 / (hr*hr) + 1.0 / (2.0*j*hr*hr));
			L[s][s - 1] = (1.0 / (hr*hr) - 1.0 / (2.0*j*hr*hr));
			L[s][s + SIZE] = 1.0 / (j*j*hr*hr*hphi*hphi);
			L[s][ SIZE*(SIZE-1) + j ] = 1.0 / (j*j*hr*hr*hphi*hphi);*/
			temp = -2.0*(1.0 / (hr*hr) + 1.0 / (j*j*hr*hr*hphi*hphi));
			Lap.Mat[s].push_back(temp);
			Lap.index[s].push_back(s);


			temp = (1.0 / (hr*hr) + 1.0 / (2.0*j*hr*hr));
			Lap.Mat[s].push_back(temp);
			Lap.index[s].push_back(s + 1);



			temp = (1.0 / (hr*hr) - 1.0 / (2.0*j*hr*hr));
			Lap.Mat[s].push_back(temp);
			Lap.index[s].push_back(s - 1);


			temp = 1.0 / (j*j*hr*hr*hphi*hphi);
			Lap.Mat[s].push_back(temp);
			Lap.index[s].push_back(s + SIZE);


			temp = 1.0 / (j*j*hr*hr*hphi*hphi);
			Lap.Mat[s].push_back(temp);
			Lap.index[s].push_back(SIZE*(SIZE - 1) + j);

		}

		for (i = 0; i < SIZE; i++)
		{
			s = i*SIZE;

			for (size_t SS = 0; SS < Lap.Mat[s].size(); SS++)
			{
				/*L[s][SS] = 0;
				L[s + SIZE - 1][SS] = 0;*/
				Lap.Mat[s][SS] = 0;
				Lap.Mat[s + SIZE - 1][SS] = 0;


			}

			/*L[s][s] = 1;
			L[s + SIZE - 1][s + SIZE - 1] = 1;*/
			Lap.Mat[s].resize(1);
			Lap.index[s].resize(1);
			Lap.Mat[s + SIZE - 1].resize(1);
			Lap.index[s + SIZE - 1].resize(1);

			Lap.Mat[s][0] = 1;
			Lap.Mat[s + SIZE - 1][0] = 1;
			Lap.index[s][0] = s;
			Lap.index[s + SIZE - 1][0] = s + SIZE - 1;

			//b[s] = 1;
			//b[s + SIZE - 1] = 1;


			//cout << "L[s][s] = " << L[s][s] << "; b[s] = " << b[s] << endl;
			//cout << "L[s + ][s + ] = " << L[s + SIZE - 1][s + SIZE - 1] << "; b[s +] = " << b[s + SIZE - 1] << endl;
		}
		cout << "sparse matrix corrected" << endl;
		/*
		for (size_t k = 1; k < SIZE*SIZE - 1; k++)
		{
		cout << L[k][k - 1] << " " << L[k][k] << " " << L[k][k + 1] <<endl;
		}
		*/
	}

	void SolveLaplas(double w, double EPS)
	{
		cout << "open SolveLaplas" << endl;
		if (Lap.Mat.size() != SIZE*SIZE)
		{
			
			LaplasPolar(&Lap, SIZE, SIZE, hr, hphi);
			u.Resize(SIZE*SIZE);
			CorrectLaplasPolarMatrix();
			
		}
		Vector<double> x0 = u;
		RelaxSparse(&Lap, &b, &x0, &u, w, EPS);
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

