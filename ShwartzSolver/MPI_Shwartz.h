#pragma once
#include <mpi.h>
#include <vector>
#include "ShwartzSolver.h"

void MPI_Shwartz(int *rank, int Wproc,
	double border_x, double border_y, int N_x, int M_y, double length_horiz, 
	double length_vert, double radius)
{

	int size = 0;
	int rc = 0;
	MPI_Status Status;

	//-----------


	MPI_Comm_rank(MPI_COMM_WORLD, rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int MainProc = size - 1;


	vector<int> ranks(Wproc);
	for (int i = 0; i < Wproc; i++)
		ranks[i] = i;

	/*
	-----------------------------
	*/
	MPI_Group slaves;
	MPI_Group world;
	MPI_Comm slaves_comm;
	MPI_Comm_group(MPI_COMM_WORLD, &world);
	MPI_Group_incl(world, Wproc, ranks.data(), &slaves);
	MPI_Comm_create(MPI_COMM_WORLD, slaves, &slaves_comm);
	//-------------------
	
	ShwartzSolver solver(border_x, border_y, N_x, M_y, length_horiz, length_vert, radius);
	

	if((*rank) == 0)
	{
	
		solver.SolveInHorizontal()
	}

}