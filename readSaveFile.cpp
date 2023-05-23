/*
 * readSaveFile.cpp
 *
 *  Created on: 7 oct. 2022
 *      Author: Romain Simon
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "cnpy.h"
#include "readSaveFile.h"
#include <iomanip>

std::string quote( const std::string& s )
{
    std::ostringstream ss;
    ss << std::quoted( s );
    return ss.str();
}

template<typename T>
std::vector<T> flatten(std::vector<std::vector<T>> const &vec)
{
    std::vector<T> flattened;
    for (auto const &v: vec) {
        flattened.insert(flattened.end(), v.begin(), v.end());
    }
    return flattened;
}


posRad readXYZ(const std::string& path)
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

    int sizeArray{static_cast<int>(radiusArray.size())};
    long unsigned int uSizeArray{radiusArray.size()};
    cnpy::npz_save(path, "n_particles", &sizeArray, {1}, "w");
    cnpy::npz_save(path, "length_cube", &lengthCube, {1}, "a");
    cnpy::npz_save(path, "arr_radius", &radiusArray[0], {uSizeArray}, "a");
    cnpy::npz_save(path, "arr_molecule_type", &moleculeType[0], {uSizeArray}, "a");
    cnpy::npz_save(path, "arr_position", &flatten(positionArray)[0], {uSizeArray, 3}, "a");
}

void saveDoubleTXT(const double& number, const std::string& path)
/*
 * This function saves the system energy in a txt file
 */
{
	std::ofstream fOut;
	fOut.open(path, std::ios_base::app);
	fOut << number;
	fOut << "\n";
	fOut.close();
}

void saveDisplacement(const std::vector<std::vector<double>>& dispMatrix, const std::string& path)
{
    long unsigned int uSizeArray{dispMatrix.size()};
    cnpy::npy_save(path, &flatten(dispMatrix)[0], {uSizeArray, 3}, "w");
}
