/*
 * energy.cpp
 *
 *  Created on: 3 oct. 2022
 *      Author: romainsimon
 */
#include <iostream>
#include <vector>
#include <math.h>

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
		return pow ((sigma / distance), 12) - pow ((sigma / distance), 6);
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

		if (diff > (halfLengthCube))
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
	std::vector<double> positionA(3);
	std::vector<double> positionB(3);
	int positionArraySize {static_cast<int>(positionArray.size())};
	for (int i = 0; i < positionArraySize - 1; i++) //Outer loop for rows
    {
		positionA =  {positionArray[i][0], positionArray[i][1], positionArray[i][2]} ;

		for (int j = i + 1; j < positionArraySize; j++) //inner loop for columns
        {
			positionB = {positionArray[j][0], positionArray[j][1], positionArray[j][2]};
			double distance { sqrt (squareDistancePair (positionA, positionB, lengthCube))};
        	energy += ljPotential(distance, 2 * radiusArray[j], 2 *  radiusArray[i], rc);
        }
    }
	return energy;
}

