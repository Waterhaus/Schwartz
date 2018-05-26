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
	Vector<double> *u;
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
		 Circle.Local_X_size = 2 * GetIndexX(( radius - Rx)) + 1;
		 Circle.Local_Y_size = 2 * GetIndexY((radius - Ry)) + 1;

	 }

	 Matrix<int> CreateMap()
	 {
		 Matrix<int> A(X_size, Y_size);


		 //circle
		 double EPS = 0.1;
		 double R = GetIndexX(Radius - Rx) - EPS;
		 double I = 0;
		 int index = 0;


		 for (size_t y = Circle.global_y; y < Circle.global_y + Circle.Local_Y_size; y++)
		 {
			 A[y][Circle.global_x] = 7;
			 A[y][Circle.global_x + Circle.Local_X_size] = 7;
		 }

		 for (size_t x = Circle.global_x; x < Circle.global_x + Circle.Local_X_size; x++)
		 {
			 A[Circle.global_y][x] = 7;
			 A[Circle.global_y + Circle.Local_Y_size][x] = 7;
		 }


		 
		 for (size_t y = Circle.global_y; y < Circle.global_y + Circle.Local_Y_size + 1; y++)
		 {
			 for (size_t x = Circle.global_x; x < Circle.global_x + Circle.Local_X_size + 1; x++)
			 {
				 I = sqrt((double)((y - center_y)*(y - center_y)) + (double)((x - center_x)*(x - center_x)));
				 if (I <= R)
				 {
					 A[y][x] = 1;

				 }
				 else {
					 A[y][x] = 0;

				 }

				 if ((y == Circle.global_y  || y == Circle.global_y + Circle.Local_Y_size - 1) && A[y][x] == 1)
				 {
					 A[y][x] = 3;

				 }
				 if ((x == Circle.global_x + 1|| x == Circle.global_x + Circle.Local_X_size - 1) && A[y][x] == 1)
				 {
					 A[y][x] = 3;

				 }

				 if (x > Circle.global_x)
				 {
					 if (A[y][x] == 1 && A[y][x - 1] == 0)
					 {
						 A[y][x] = 3;

					 }
					 if (A[y][x] == 0 && A[y][x - 1] == 1) {
						 A[y][x - 1] = 3;

					 }

				 }


				 index++;


			 }
		 }

		 


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

		 return A;
		
	 }

	 void CreateRigthPartH(Gebeat G, Vector<double> *b)
	 {
		 int N_x = G.Local_X_size;
		 int M_y = G.Local_Y_size;

		  
		 int lbord_x = GetIndexX(-sqrt(Radius*Radius - ly*ly));
		 int rbord_x = GetIndexX(sqrt(Radius*Radius - ly*ly));
		 int s = 0;

		 for (size_t x = 0; x < N_x; x++)
		 {
			 
			 for (size_t y = 0; y < M_y; y++)
			 {
				 s = x*N_x + y;
				 if (x == 0 || x == N_x - 1 || y == 0 || y == M_y - 1) 
				 {
					 if (x <= lbord_x || x >= rbord_x) 
					 {
						 (*b)[s] = g(x*hx - Rx, y*hy - Ry);
					 
					 }
					 else (*b)[s] = 0;
				 }

				 (*b)[s] = f(x*hx - Rx, y*hy - Ry);
			 }
		 }

		
	 }



	 Vector<double> CreateRigthPartV(Gebeat G)
	 {
		 int N_x = G.Local_X_size;
		 int M_y = G.Local_Y_size;

		 Vector<double> b(N_x*M_y - 1);

		 int lbord_y = GetIndexY(-sqrt(Radius*Radius - lx*lx));
		 int rbord_y = GetIndexY(sqrt(Radius*Radius - lx*lx));
		 int s = 0;

		 for (size_t x = 0; x < N_x; x++)
		 {

			 for (size_t y = 0; y < M_y; y++)
			 {
				 s = x*N_x + y;
				 if (x == 0 || x == N_x - 1 || y == 0 || y == M_y - 1)
				 {
					 if (y <= lbord_y || y >= rbord_y)
					 {
						 b[s] = g(x*hx - Rx, y*hy - Ry);

					 }
					 else b[s] = 0;
				 }

				 b[s] = f(x*hx - Rx, y*hy - Ry);
			 }
		 }
	 }

	 void SolveInHorizontal(Gebeat *horiz,double h_x, double h_y)
	 {
		 int N_x = horiz->Local_X_size;
		 int M_y = horiz->Local_Y_size;

		 int SIZE = N_x*M_y;
		 
		 Matrix<double> A = Laplas(N_x,M_y, h_x, h_y);
		 Vector<double> b(SIZE);
		 CreateRigthPartH((*horiz), &b);
		 Vector<double> x0(SIZE);
		 Vector<double> x(SIZE);
		 RelaxMini(&A, &b, &x0, 0.5, 0.0001, M_y);
		 
		 horiz->u = &x;
	 }


};