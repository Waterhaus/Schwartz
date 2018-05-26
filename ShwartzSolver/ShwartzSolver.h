#pragma once
#include "Vector.h";
#include "Matrix.h"
#include <vector>

using namespace std;

class Gebeat 
{
	int global_x, global_y;
	int Local_Y_size, Local_X_size;
	double Parameter;
};

class ShwartzSolver {

public:
	double Rx;
	double Ry;
	int Y_size;
	int X_size;

	double hx, hy;

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

	 ShwartzSolver()
	 {
		 Rx = Ry = 0;
		 Y_size = X_size = 0;
		 hx = hy = 0;
	 }
	 ShwartzSolver(double border_x, double border_y, int N_x, int M_y)
	 {
		 Rx = border_x;
		 Ry = border_y;
		 X_size = N_x;
		 Y_size = M_y;
		 hx = GetStep(N_x, -Rx, Rx);
		 hy = GetStep(M_y, -Ry, Ry);
	 }
};