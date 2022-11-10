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
#include <limits>


double ljPotential(double squareDistance, double sigmaA, double sigmaB, double squareRc, double shift)
/*
 *This function is calculates the Lennard-Jones potential between two particles.
	We consider a cut off that is equal to 2.5 times the mean radius of the two particles considered
 */
{
	double squareSigma { pow (((sigmaA + sigmaB) / 2 ), 2) };

	if (squareDistance > squareRc * squareSigma) // We consider a cut-off radius which is the threshold maximum distance of interaction between two particles
	{
		return 0;
	}


	else
	{
		return 4 * (pow ((squareSigma / squareDistance), 6) - pow ((squareSigma / squareDistance), 3) + shift);
	}
}

double fenePotential(double squareDistance, double sigmaA, double sigmaB, double squareR0, double feneK)
{
	double squareSigma { pow (((sigmaA + sigmaB) / 2 ), 2) };
	squareR0 = squareR0 * squareSigma;

	if (squareDistance > squareR0) // We consider a cut-off radius which is the threshold maximum distance of interaction between two particles
	{
		return std::numeric_limits<int>::max();
	}


	else
	{
		feneK = feneK / squareSigma;
		return -0.5 * feneK * squareR0 * log(1 - squareDistance / squareR0);
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

double energySystem(std::vector<std::vector<double>> positionArray, std::vector<double> radiusArray, double squareRc, double lengthCube)
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
			double squareDistance { squareDistancePair (positionArray[i], positionArray[j], lengthCube)};
        	energy += ljPotential(squareDistance, 2 * radiusArray[j], 2 * radiusArray[i], squareRc, 0.);
        }
    }
	return energy;
}


double energyParticle(int indexParticle, std::vector<double> positionParticle,
		std::vector<std::vector<double>> positionArray, std::vector<int> neighborIList, std::vector<double> radiusArray,
		double squareRc, double lengthCube)
{
	double energy { 0. };
	double particleRadius = radiusArray[indexParticle];
	int neighborIListSize { static_cast<int>(neighborIList.size()) };

	for (int i = 0; i < neighborIListSize; i++)
	{
		int realIndex = neighborIList[i];
		double squareDistance { squareDistancePair (positionParticle, positionArray[realIndex], lengthCube)};

		if (realIndex == indexParticle)
		{
			continue;
		}

		energy += ljPotential(squareDistance, 2 * particleRadius, 2 * radiusArray[realIndex], squareRc, 0.);
	}
	return energy;
}


double energySystemPolymer(std::vector<std::vector<double>> positionArray, std::vector<double> radiusArray,
						   std::vector<std::vector<int>> bondsMatrix, double squareRc, double lengthCube,
						   double squareR0, double feneK)
/*
 *This function calculates the total potential energy of the system ie. the Lennard-Jones potential between each pair of particles.
 */
{
	double energy { 0. };
	int positionArraySize {static_cast<int>(positionArray.size())};

	for (int i = 0; i < positionArraySize - 1; i++) //Outer loop for rows
    {

    	std::vector<int> bondsI ( bondsMatrix[i] );
    	int bondsISize { static_cast<int>(bondsI.size()) };

		for (int j = i + 1; j < positionArraySize; j++) //inner loop for columns
        {
			double squareDistance { squareDistancePair (positionArray[i], positionArray[j], lengthCube)};

			energy += ljPotential(squareDistance, 2 * radiusArray[j], 2 * radiusArray[i], squareRc, 127. / 4096);

        }
		for (int j = 0; j < bondsISize; j++)
		{
			int realIndex = bondsI[j];

			if ((realIndex == -1) || (realIndex < i))
			{
				continue;
			}

			double squareDistance { squareDistancePair (positionArray[i], positionArray[realIndex], lengthCube)};
			energy += fenePotential(squareDistance, 2 * radiusArray[i], 2 * radiusArray[realIndex], squareR0, feneK);

        }
    }
	return energy;
}


double energyParticlePolymer (int indexParticle, std::vector<double> positionParticle,
								std::vector<std::vector<double>> positionArray, std::vector<int> neighborIList,
								std::vector<double> radiusArray, std::vector<int> bondsI,
								double squareRc, double lengthCube, double squareR0, double feneK)
{
	double energy { 0. };
	double particleRadius = radiusArray[indexParticle];
	int neighborIListSize { static_cast<int>(neighborIList.size()) };
	int bondsISize { static_cast<int>(bondsI.size()) };

	for (int i = 0; i < neighborIListSize; i++)
	{
		int realIndex = neighborIList[i];
		double squareDistance { squareDistancePair (positionParticle, positionArray[realIndex], lengthCube)};

		if (realIndex == indexParticle)
		{
			continue;
		}

		energy += ljPotential(squareDistance, 2 * particleRadius, 2 * radiusArray[realIndex], squareRc, 127. / 4096);
	}

	for (int i = 0; i < bondsISize; i++)
	{
		int realIndex = bondsI[i];

		if ((realIndex == indexParticle) || (realIndex == -1))
		{
			continue;
		}

		double squareDistance { squareDistancePair (positionParticle, positionArray[realIndex], lengthCube)};
		energy += fenePotential(squareDistance, 2 * particleRadius, 2 * radiusArray[realIndex], squareR0, feneK);
	}

	return energy;
}

double rforce(double squareDistance, double sigmaA, double sigmaB)
{
	double squareSigma { pow (((sigmaA + sigmaB) / 2 ), 2) };
	double rforceAB { 24 * (pow ((squareSigma / squareDistance), 3) - 2 * pow ((squareSigma / squareDistance), 6))};
	return rforceAB;

}

double pressureSystem(double temp, std::vector<std::vector<double>> positionArray, std::vector<double> radiusArray, double lengthCube)
{
	double sum { 0 };
	int positionArraySize {static_cast<int>(positionArray.size())};
	for (int i = 0; i < positionArraySize - 1; i++)
	{
		for (int j = i + 1; j < positionArraySize; j++) //inner loop for columns
        {
			double squareDistance { squareDistancePair (positionArray[i], positionArray[j], lengthCube)};
			sum += rforce(squareDistance, 2 * radiusArray[j], 2 * radiusArray[i]);
        }

	}
	return 1 - sum / (3. * positionArraySize * temp);
}
//std::vector<std::vector<int>> createNeighborList(std::vector<std::vector<double>> positionArray, double skin, double lengthCube)
//{
//	int positionArraySize { static_cast<int>(positionArray.size()) };
//	std::vector<std::vector<int>> neighborList( positionArraySize );
//	double squareSkin {pow(skin, 2)};
//
//	for (int i = 0; i < positionArraySize - 1; i++)
//	{
//		std::vector<double> positionParticle = positionArray[i];
//
//		for (int j = i + 1; j < positionArraySize; j++)
//		{
//
//			if (squareDistancePair ( positionParticle, positionArray[j], lengthCube ) < squareSkin)
//			{
//				neighborList[i].push_back(j);
//				neighborList[j].push_back(i);
//			}
//		}
//
//		}
//return neighborList;
//}



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








