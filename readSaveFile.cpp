/*
 * readSaveFile.cpp
 *
 *  Created on: 7 oct. 2022
 *      Author: romainsimon
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "readSaveFile.h"

inputVar readInput(std::string path)
{
	inputVar initVar;
	std::string null;
	std::ifstream infile(path);

	if (!infile.is_open())
		std::cout << "Error opening file" ;

	infile >> initVar.temp;
	infile >> null;
	infile >> null;
	infile >> initVar.density;
	infile >> null;
	infile >> null;
	infile >> initVar.rc;
	infile >> null;
	infile >> null;
	infile >> initVar.rskin;
	infile >> null;
	infile >> null;
	infile >> initVar.rbox;
	infile >> null;
	infile >> null;
	infile >> initVar.neighUpdate;
	infile >> null;
	infile >> null;
	infile >> initVar.numberIteration;
	return initVar;
}
posRad readXYZ(std::string path)
{
	std::ifstream infile(path);

	if (!infile.is_open())
		std::cout << "Error opening file" ;

	int row{};
	infile >> row;
    int col = 4;

    std::vector<std::vector <double>> positionArray(row, std::vector<double>(3));
    std::vector<double> radiusArray(row);


	//Defining the loop for getting input from the file

	for (int r = 0; r < row; r++) //Outer loop for rows
	{
		for (int c = 0; c < col; c++) //inner loop for columns
		{
			if (c == 0)
			{
				infile >> radiusArray[r];
			}

			else
			{
				infile >> positionArray[r][c-1]; //Take input from file and put into positionArray
			}
    	}

	}

	infile.close();
	posRad initPosRad;
	initPosRad.radVector = radiusArray;
	initPosRad.posMatrix = positionArray;
	return initPosRad;
}



void saveInXYZ(std::vector<std::vector<double>>& positionArray, std::vector<double> radiusArray, std::string path)
{

	/*
	 * This function saves in a .xyz file the radius and the position of each particle.
	  */


	std::ofstream fout(path);
    fout << positionArray.size();
	fout <<  "\n";
    fout <<  "\n";
    int it {0};
	for(auto const& x : positionArray)
	{

		int j { 0 };
		for(auto const& i : x)
		{

			if (j == 2)
			{
				fout << i;
				fout << "\n";
			}
			else if (j==0)
			{
				fout << radiusArray[it];
				fout << " ";
				fout << i;
				fout << " ";
			}
			else
			{
				fout << i;
				fout << " ";
			}

			j += 1;

		}
		it += 1;
	}
	fout.close();
}

void saveEnergyTXT(double energy, std::string path)
/*
 * This function saves the system energy in a txt file
 */
{
	std::ofstream outE;
	outE.open(path, std::ios_base::app);

	outE << energy;
	outE << "\n";
	outE.close();

}

void printing(const std::vector<std::vector<double>>& matrix)
/*
 * This function prints a two dimensional array (matrix).
 */
{
	for (auto i: matrix)
	{
		for (auto j: i)
			std::cout << j << " " ;
	}
}

