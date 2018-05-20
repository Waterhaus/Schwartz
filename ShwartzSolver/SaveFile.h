#pragma once
#include <iostream>
#include <fstream>
#include "Vector.h"
#include <string> 

using namespace std;

void SaveToFile(string file, Vector<double> *anser, int N)
{
	ofstream myfile;
	myfile.open(file + ".txt");
	myfile << (*anser).Size() <<"\n";
	myfile << N << "\n";

	for (size_t i = 0; i < anser->Size(); i++)
	{
		myfile << (*anser)[i] << "\n";
	}

	myfile.close();
	cout << "save to " + file << endl;
	


}