// ShwartzSolver.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "Relax.h"
#include "LaplasSolver.h"
#include "SaveFile.h"
#include "DomainMethod.h"
using namespace std;

double Intercetion(double R, double Val)
{
	return sqrt(R*R - Val*Val);
}

int main()
{
	int SIZE = 6;
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

	CircleDomain circle("circle", SIZE, x, y, Radius, 10);
	//circle.SolveLaplas(1.1, 0.0001);
	
	double boxside = 2 * VerX;
	BoxDomain rightbox("rightbox", SIZE, VerX, -HorY, VerX + boxside, -HorY + boxside, Radius);
	
	BoxDomain topbox("topbox", SIZE, -VerX,HorY ,-VerX  + boxside,HorY  + boxside, Radius);
	BoxDomain leftbox("leftbox", SIZE, -VerX - boxside,-HorY ,-VerX,-HorY +boxside, Radius);
	BoxDomain bottombox("bottombox", SIZE, -VerX,-HorY - boxside ,-VerX + boxside, -HorY , Radius);
	
	for (size_t i = 0; i < SIZE; i++)
	{
		cout << i << "  - " << circle.IsInside(i) << endl;
	}

	double EPS = 0.00001;
	circle.SolveLaplas(0.1, 0.00001);
	cout << " L*identity = "<< circle.L*Vector<double>::Identity(circle.L.Size()) << endl;

	cout << Relax(circle.L, circle.b, circle.u, 0.78, 0.00000001) << endl;
	RelaxFastBad(&circle.L, &circle.b, &circle.u, &circle.u, 0.78, 0.00000001, SIZE);
	//cout << circle.u << endl;
	/*rightbox.SolveLaplas(1.4, EPS);
	topbox.SolveLaplas(1.4, EPS);
	leftbox.SolveLaplas(1.4, EPS);
	bottombox.SolveLaplas(1.4, EPS);
	*/
	circle.Save();
	/*rightbox.Save();
	topbox.Save();
	leftbox.Save();
	bottombox.Save();
*/
	cout << "max in circle = " << circle.u.normInf() << endl;
	//cout << circle.u << endl;

	char a;
	cin >> a;
    return 0;
}

