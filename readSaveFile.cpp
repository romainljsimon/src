/*
 * readSaveFile.cpp
 *
 *  Created on: 7 oct. 2022
 *      Author: romainsimon
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "readSaveFile.h"

inputVar readInput(const std::string& path)
{
	inputVar initVar;
	std::string null;
	std::ifstream infile(path);

	if (!infile.is_open())
		std::cout << "Error opening file" ;

	infile >> initVar.simulationMol;
	infile >> null;
	infile >> null;
	infile >> initVar.neighborMethod;
	infile >> null;
	infile >> null;
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
	infile >> initVar.timeSteps;
	infile >> null;
	infile >> null;

	if (initVar.simulationMol == "polymer")
	{
		infile >> initVar.r0;
		infile >> null;
		infile >> null;
		infile >> initVar.feneK;
	}

	else
	{
		initVar.r0 = 0;
		initVar.feneK = 0;
	}
	return initVar;
}
posRad readXYZ(const std::string& path, const std::string& simulationMol)
{
	std::ifstream infile(path);

	if (!infile.is_open())
		std::cout << "Error opening file" ;

	int row{};
	infile >> row;
	int col { 5 };

	std::string line;
	getline(infile, line);
	getline(infile, line);


	std::vector<std::vector <double>> positionArray(row, std::vector<double>(3));
    std::vector<double> radiusArray (row);
    std::vector<int> moleculeType (row , 1);


	//Defining the loop for getting input from the file

	for (int r = 0; r < row; r++) //Outer loop for rows
	{
		for (int c = 0; c < col; c++) //inner loop for columns
		{
			if (c == 0)
			{
				infile >> moleculeType[r];
			}
			else if (c == 1)
			{
				infile >> radiusArray[r];
			}

			else
			{
				infile >> positionArray[r][c - (col - 3)]; //Take input from file and put into positionArray
			}
    	}

	}

	infile.close();
	posRad initPosRad;
	initPosRad.radVector = radiusArray;
	initPosRad.posMatrix = positionArray;
	initPosRad.moleculeType = moleculeType;
	return initPosRad;
}

std::vector<std::vector<int>> readBondsTXT(const std::string& path)
{
	std::ifstream infile(path);

	if (!infile.is_open())
		std::cout << "Error opening file";


	int row{};
	infile >> row;
	int col{2};
	std::vector<std::vector <int>> bondsMatrix(row, std::vector<int>(2));

	for (int r = 0; r < row; r++) //Outer loop for rows
	{
		for (int c = 0; c < col; c++) //inner loop for columns
		{
			infile >> bondsMatrix[r][c]; //Take input from file and put into bondsMatrix

    	}

	}

	infile.close();

	return bondsMatrix;
}

void saveInXYZ(const std::vector<std::vector<double>>& positionArray, const std::vector<double>& radiusArray,
			   const std::vector<int>& moleculeType, const double& lengthCube,  const std::string& path)
{

	/*
	 * This function saves in a .xyz file the radius and the position of each particle.
	  */


	std::ofstream fout(path);
    fout << positionArray.size();
	fout <<  "\n";
	std::string lengthStr = std::to_string(lengthCube);
    fout << "Lattice=";
    fout << '"';
    fout << lengthStr;
    fout << " 0.0 0.0 0.0 ";
    fout << lengthStr;
    fout << " 0.0 0.0 0.0 ";
    fout << lengthStr;
    fout << '"';
    fout << " Properties=type:I:1:radius:R:1:pos:R:3";
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
				fout << moleculeType[it];
				fout << " ";
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

void saveDoubleTXT(const double& number, const std::string& path)
/*
 * This function saves the system energy in a txt file
 */
{
	std::ofstream fout;
	fout.open(path, std::ios_base::app);
	fout << number;
	fout << "\n";
	fout.close();

}

void saveDisplacement(const std::vector<std::vector<double>>& dispMatrix, const std::string& path)
{
	std::ofstream fout(path);

	for(auto const& x : dispMatrix)
	{
		int j { 0 };
		for(auto const& i : x)
		{
			fout << i;
			fout << " ";
			if (j == 2)
			{
				fout << "\n";
			}

			j += 1;

		}
	}
	fout.close();
}

void printing(const std::vector<std::vector<double>>& matrix)
/*
 * This function prints a two-dimensional array (matrix).
 */
{
	for (auto i: matrix)
	{
		for (auto j: i)
			std::cout << j << " " ;
	}
}

