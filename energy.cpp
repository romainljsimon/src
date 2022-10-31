/*
 * energy.cpp
 *
 *  Created on: 3 oct. 2022
 *      Author: romainsimon
 */
#include <iostream>
#include <vector>
#include <math.h>
#include "util.h"

double ljPotential(double distance, double sigmaA, double sigmaB, double rc)
/*
 *This function is calculates the Lennard-Jones potential between two particles.
	We consider a cut off that is equal to 2.5 times the mean radius of the two particles considered
 */
{
	double sigma { (sigmaA + sigmaB) / 2 };

	if (distance > rc * sigma) // We consider a cut-off radius which is the threshold maximum distance of interaction between two particles
	{
		return 0;
	}


	else
	{
		return 4 * (pow ((sigma / distance), 12) - pow ((sigma / distance), 6));
	}
}

double squareDistancePair(std::vector<double> positionA,  std::vector<double> positionB, double lengthCube)
/*
 *This function calculates the distance between two particles considering the periodic boundary conditions.
 */
{
	double squareDistance { 0. };
	double halfLengthCube { lengthCube / 2 };

	for (int i = 0; i < static_cast<int>(positionA.size()); i++)
	{
		double diff { positionA[i] - positionB[i] };

		if (diff > halfLengthCube)
			diff -= lengthCube;

		else if (diff < - halfLengthCube)
			diff += lengthCube;

		squareDistance += pow (diff, 2);
	}
	return squareDistance;
}

double energySystem(std::vector<std::vector<double>> positionArray, std::vector<double> radiusArray, double rc, double lengthCube)
/*
 *This function calculates the total potential energy of the system ie. the Lennard-Jones potential between each pair of particles.
 */
{
	double energy { 0. };
	int positionArraySize {static_cast<int>(positionArray.size())};

	for (int i = 0; i < positionArraySize - 1; i++) //Outer loop for rows
    {

		for (int j = i + 1; j < positionArraySize; j++) //inner loop for columns
        {
			double distance { sqrt (squareDistancePair (positionArray[i], positionArray[j], lengthCube))};
        	energy += ljPotential(distance, 2 * radiusArray[j], 2 * radiusArray[i], rc);
        }
    }
	return energy;
}


double energyParticle(int wu, int indexParticle, std::vector<double> positionParticle,
		std::vector<std::vector<double>> positionArray, std::vector<int> neighborIList, std::vector<double> radiusArray,
		double rc, double rskin, double lengthCube)
{
	double energy { 0. };
	double particleRadius = radiusArray[indexParticle];
	int neighborIListSize { static_cast<int>(neighborIList.size()) };

	for (int i = 0; i < neighborIListSize; i++)
	{
		int realIndex = neighborIList[i];
		double distance { sqrt (squareDistancePair (positionParticle, positionArray[realIndex], lengthCube)) };

		if (realIndex == indexParticle)
		{
			continue;
		}

		energy += ljPotential(distance, 2 * particleRadius, 2 * radiusArray[realIndex], rc);
	}
	return energy;
}



std::vector<std::vector<int>> createNeighborList(std::vector<std::vector<double>> positionArray, double skin, double lengthCube)
{
	int positionArraySize { static_cast<int>(positionArray.size()) };
	std::vector<std::vector<int>> neighborList( positionArraySize );
	double squareSkin {pow(skin, 2)};

	for (int i = 0; i < positionArraySize - 1; i++)
	{
		std::vector<double> positionParticle = positionArray[i];

		for (int j = i + 1; j < positionArraySize; j++)
		{

			if (squareDistancePair (positionParticle, positionArray[j], lengthCube) < squareSkin)
			{
				neighborList[i].push_back(j);
				neighborList[j].push_back(i);
			}
		}

		}
return neighborList;
}



//
//bool changeList(std::vector<double> positionDouble, double skin, double lengthCube)
//{
//	int positionArraySize { static_cast<int>(positionArray.size()) };
//	int cellListSize {static_cast<int>(lengthCube / skin)};
//	std::vector<std::vector<std::vector<std::vector<int>>>> cellList (cellListSize, std::vector<std::vector<std::vector<int>>>
//            														 (cellListSize, std::vector<std::vector<int>>
//            														 (cellListSize, std::vector<int>
//                                                                     )));
//	std::vector<std::vector <int>> particleIndexCell(positionArraySize, std::vector<int>(3));
//
//
//	for (int i = 0; i < positionArraySize - 1; i++)
//	{
//		std::vector<double> positionParticle = positionArray[i];
//		std::vector<int> boxIndex static_cast<int>(divideVectorByScalar(positionParticle, skin));
//		cellList[boxIndex[0]][boxIndex[1]][boxIndex[2]].push_back(i);
//		particleIndexCell[i] = boxIndex;
//
//
//		}
//	return cellList;
//}
//








