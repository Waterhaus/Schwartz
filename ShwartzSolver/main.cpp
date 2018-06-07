// ShwartzSolver.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "Relax.h"
#include "LaplasSolver.h"
#include "SaveFile.h"
#include "BoxDomain.h"
#include "CircleDomain.h"
//using namespace std;

double Intercetion(double R, double Val)
{
	return sqrt(R*R - Val*Val);
}

int main()
{
	int SIZE = 40;
	double Radius = 4.5;

	double HorY = 2;
	double HorX = Intercetion(Radius, HorY);
	double VerX = 2;
	double VerY = Intercetion(Radius, VerX);
	vector<double> x;
	vector<double> y;

	//rightbox
	x.push_back(HorX);
	y.push_back(HorY);

	x.push_back(HorX);
	y.push_back(-HorY);

	//topbox
	x.push_back(VerX);
	y.push_back(VerY);

	x.push_back(-VerX);
	y.push_back(VerY);

	//leftbox
	x.push_back(-HorX);
	y.push_back(HorY);

	x.push_back(-HorX);
	y.push_back(-HorY);

	//bottombox
	x.push_back(-VerX);
	y.push_back(-VerY);

	x.push_back(VerX);
	y.push_back(-VerY);

	CircleDomain circle("circle", SIZE, x, y, Radius, 1);
	
	
	double boxside = 2 * VerX;
	BoxDomain rightbox("rightbox", SIZE, VerX, -HorY, VerX + boxside, -HorY + boxside, Radius);
	
	BoxDomain topbox("topbox", SIZE, -VerX,HorY ,-VerX  + boxside,HorY  + boxside, Radius);
	BoxDomain leftbox("leftbox", SIZE, -VerX - boxside,-HorY ,-VerX,-HorY +boxside, Radius);
	BoxDomain bottombox("bottombox", SIZE, -VerX,-HorY - boxside ,-VerX + boxside, -HorY , Radius);
	
	/*for (size_t i = 0; i < SIZE; i++)
	{
		cout << i << "  - " << circle.IsInside(i) << endl;
	}*/

	double EPS = 0.00001;
	circle.SolveLaplas(1.2, 0.00001);
	rightbox.SolveLaplas(1.4, EPS);
	topbox.SolveLaplas(1.4, EPS);
	leftbox.SolveLaplas(1.4, EPS);
	bottombox.SolveLaplas(1.4, EPS);

	vector<double> temp(rightbox.out_x.size());
	circle.GetInterpolatedPointsXY(rightbox.out_x.data(), rightbox.out_x.data(), temp.data(), temp.size());
	rightbox.UpdateBorder(temp.data(), temp.size());

	temp.resize(topbox.out_x.size());
	circle.GetInterpolatedPointsXY(topbox.out_x.data(), topbox.out_x.data(), temp.data(), temp.size());
	topbox.UpdateBorder(temp.data(), temp.size());

	temp.resize(leftbox.out_x.size());
	circle.GetInterpolatedPointsXY(leftbox.out_x.data(), leftbox.out_x.data(), temp.data(), temp.size());
	leftbox.UpdateBorder(temp.data(), temp.size());

	temp.resize(bottombox.out_x.size());
	circle.GetInterpolatedPointsXY(bottombox.out_x.data(), bottombox.out_x.data(), temp.data(), temp.size());
	bottombox.UpdateBorder(temp.data(), temp.size());


	vector<double> temp_r(circle.rigthSec.x.size());
	vector<double> temp_l(circle.leftSec.x.size());
	vector<double> temp_t(circle.topSec.x.size());
	vector<double> temp_b(circle.bottomSec.x.size());


	rightbox.GetInterpolatedPoints(circle.rigthSec.x.data(), circle.rigthSec.y.data(), temp_r.data(), temp_r.size());
	leftbox.GetInterpolatedPoints(circle.leftSec.x.data(), circle.leftSec.y.data(), temp_l.data(), temp_l.size());
	topbox.GetInterpolatedPoints(circle.topSec.x.data(), circle.topSec.y.data(), temp_t.data(), temp_t.size());
	bottombox.GetInterpolatedPoints(circle.bottomSec.x.data(), circle.bottomSec.y.data(), temp_b.data(), temp_b.size());

	circle.UpdateBorder(temp_r.data(), temp_r.size(), temp_t.data(), temp_t.size(), temp_l.data(), temp_l.size(), temp_b.data(), temp_b.size());


	rightbox.SolveLaplas(1.4, EPS);
	topbox.SolveLaplas(1.4, EPS);
	leftbox.SolveLaplas(1.4, EPS);
	bottombox.SolveLaplas(1.4, EPS);

	circle.Save();
	rightbox.Save();
	topbox.Save();
	leftbox.Save();
	bottombox.Save();

	cout << "max in circle = " << circle.u.normInf() << endl;
	cout << "min in circle = " << circle.u.min() << endl;
	//cout << circle.u << endl;

	char a;
	cin >> a;
    return 0;
}

