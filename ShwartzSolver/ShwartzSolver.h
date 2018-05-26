#pragma once
#include "Vector.h";
#include "Matrix.h"
#include <vector>
#include <iostream>

using namespace std;

struct Gebeat 
{
	int global_x, global_y;
	int Local_Y_size, Local_X_size;
	Vector<double> u;
};

class ShwartzSolver {

public:
	double Rx,lx,Radius;
	double Ry,ly;
	int Y_size;
	int X_size;

	int center_x, center_y;

	double hx, hy;
	Gebeat HorizonalLine;
	Gebeat VerticalLine;
	Gebeat Circle;



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

	 int GetIndexX(double x)
	 {
		 if (x < -Rx) return 0;
		 if (x > Rx) return X_size;

		 int index = (int)((x + Rx) / hx);

		 return index;		
	 }


	 int GetIndexY(double y)
	 {
		 if (y < -Ry) return 0;
		 if (y > Ry) return Y_size;
		 int index = (int)((y + Ry) / hy);

		 return index;
	 }

	 ShwartzSolver()
	 {
		 Rx = Ry = 0;
		 Y_size = X_size = 0;
		 hx = hy = 0;
	 }
	 ShwartzSolver(double border_x, double border_y, int N_x, int M_y, double length_horiz, double length_vert, double radius)
	 {
		 Rx = border_x;
		 Ry = border_y;
		 X_size = N_x;
		 Y_size = M_y;
		 hx = GetStep(N_x, -Rx, Rx);
		 hy = GetStep(M_y, -Ry, Ry);
		 center_x = GetIndexX(0);
		 center_y = GetIndexY(0);
		 lx = length_vert;
		 ly = length_horiz;
		 Radius = radius;

		 int k_y = (int)((Ry - length_horiz) / hy);
		 HorizonalLine.Local_X_size = N_x;
		 HorizonalLine.Local_Y_size = 2*(center_y - k_y);
		 HorizonalLine.global_x = 0;
		 HorizonalLine.global_y = k_y;




		 int k_x = (int)((Rx - length_vert)/hx);
		 VerticalLine.Local_X_size = 2*(center_x - k_x);
		 VerticalLine.Local_Y_size = M_y;
		 VerticalLine.global_x = k_x;
		 VerticalLine.global_y = 0;


		 Circle.global_x = GetIndexX(-radius);
		 Circle.global_y = GetIndexY(-radius);
		 int r_x = Circle.global_x;
		 int r_y = Circle.global_y;
		 Circle.Local_X_size = 2*GetIndexX(( radius - Rx));
		 Circle.Local_Y_size = 2 * GetIndexY((radius - Ry));

	 }

	 Matrix<int> CreateMap()
	 {
		 Matrix<int> A(X_size, Y_size);

		 //horizontal
		 int bord_x = GetIndexX(- sqrt(Radius*Radius - ly*ly));
		 cout << "bx = " << bord_x << endl;
		 for (size_t x = 0; x < bord_x; x++)
		 {
			 A[HorizonalLine.global_y][x] = 1;
			 A[HorizonalLine.global_y + HorizonalLine.Local_Y_size][x] = 1;
		 }
		 bord_x = GetIndexX(sqrt(Radius*Radius - ly*ly));


		 for (size_t x = bord_x; x < X_size; x++)
		 {
			 A[HorizonalLine.global_y][x] = 1;
			 A[HorizonalLine.global_y + HorizonalLine.Local_Y_size][x] = 1;
		 }

		 //vertical
		 int bord_y = GetIndexY(-sqrt(Radius*Radius - lx*lx));
		 cout << "by = " << bord_y << endl;
		 for (size_t y = 0; y < bord_y; y++)
		 {
			 A[y][VerticalLine.global_x] = 2;
			 A[y][VerticalLine.global_x + VerticalLine.Local_X_size] = 2;
		 }
		 bord_y = GetIndexY(sqrt(Radius*Radius - lx*lx));
		 for (size_t y = bord_y; y < Y_size; y++)
		 {
			 A[y][VerticalLine.global_x] = 2;
			 A[y][VerticalLine.global_x + VerticalLine.Local_X_size] = 2;
		 }

		 //circle


		 return A;
		
	 }

};