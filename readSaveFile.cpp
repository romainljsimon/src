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
#include "readSaveFile.h"


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
            infile >> bondsMatrix[r][c]; //Take INPUT from file and put into bondsMatrix

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


    std::ofstream fOut(path);
    fOut << positionArray.size();
    fOut << "\n";
    std::string lengthStr = std::to_string(lengthCube);
    fOut << "Lattice=";
    fOut << '"';
    fOut << lengthStr;
    fOut << " 0.0 0.0 0.0 ";
    fOut << lengthStr;
    fOut << " 0.0 0.0 0.0 ";
    fOut << lengthStr;
    fOut << '"';
    fOut << " Properties=type:I:1:radius:R:1:pos:R:3";
    fOut << "\n";
    int it {0};
    for(auto const& x : positionArray)
    {

        int j { 0 };
        for(auto const& i : x)
        {

            if (j == 2)
            {
                fOut << i;
                fOut << "\n";
            }
            else if (j==0)
            {
                fOut << moleculeType[it];
                fOut << " ";
                fOut << radiusArray[it];
                fOut << " ";
                fOut << i;
                fOut << " ";
            }
            else
            {
                fOut << i;
                fOut << " ";
            }

            j += 1;

        }
        it += 1;
    }
    fOut.close();
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

void saveDoubleIntTXT(const double& number1, const int& number2, const std::string& path)
/*
 * This function saves the system energy in a txt file
 */
{
    std::ofstream fOut;
    fOut.open(path, std::ios_base::app);
    fOut << number1;
    fOut << " ";
    fOut << number2;
    fOut << "\n";
    fOut.close();

}

void saveDisplacement(const std::vector<std::vector<double>>& dispMatrix, const std::string& path)
{
    std::ofstream fOut(path);

    for(auto const& x : dispMatrix)
    {
        int j { 0 };
        for(auto const& i : x)
        {
            fOut << i;
            fOut << " ";
            if (j == 2)
            {
                fOut << "\n";
            }

            j += 1;

        }
    }
    fOut.close();
}

void printing(const std::vector<std::vector<double>>& matrix)
/*
 * This function prints a two-dimensional array (matrix).
 */
{
    for (const auto& i: matrix)
    {
        for (auto j: i)
            std::cout << j << " " ;
    }
}
