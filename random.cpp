/*
 * random.cpp
 *
 *  Created on: 3 oct. 2022
 *      Author: Romain Simon
 */
#include <random>
#include <iostream>
#include "util.h"


double randomDoubleGenerator(const double& min, const double& max)
{
	// This function returns a random double in the interval [min, max).
	// For Fixed seed comment these two next lines and un comment the 3d !!
    static std::random_device rd;									// Used to obtain the seed for the following random number generator
	static std::mt19937 mt(rd());
    //static std::mt19937 mt;
    std::uniform_real_distribution<double> dist(min, max);          // The period of this generator is 2^9337-1. This there shouldn
    return dist(mt);												// This allows to transform the unsigned int mt into a random double between [min, max).
}

int randomIntGenerator(const int& min, const int& max)
{
	// This function returns a random int from the {min, min + 1, ... , max  - 1, max} set.
	return randomDoubleGenerator(static_cast<int>(min), static_cast<int>(max) + 1);
}

std::vector<double> randomVectorDoubleGenerator(const int& vectorSize, const double& min, const double& max)
{
	std::vector<double> randomVector(vectorSize);
	for (int i = 0; i < vectorSize; i++)
	{
		randomVector[i] = randomDoubleGenerator(min, max);
	}
	return randomVector;
}

double randomNormalDoubleGenerator()
{
    static std::random_device rd;									// Used to obtain the seed for the following random number generator
    static std::mt19937 mt(rd());
    //static std::mt19937 mt;
    std::normal_distribution<double> dist(0, 1);          // The period of this generator is 2^9337-1. This there shouldn
    return dist(mt);
}

std::vector<double> randomUnitSphereVectorGenerator(const int& vectorSize)
{
    std::vector<double> randomVector(vectorSize);

    for (int i = 0; i < vectorSize; i++)
    {
        randomVector[i] = randomNormalDoubleGenerator();
    }
    vectorNormalization(randomVector.begin(), randomVector.end());
    return randomVector;
}

std::vector<double> randomRotationMatrixGenerator(const int& vectorSize)
{
    std::vector<double> randomRotationVector{randomUnitSphereVectorGenerator(vectorSize)};
    double randomAngle {randomDoubleGenerator(-M_PI,M_PI)};
    const double c = cos(randomAngle);
    const double s = sin(randomAngle);
    const double& rx {randomRotationVector[0]};
    const double& ry {randomRotationVector[1]};
    const double& rz {randomRotationVector[2]};
    std::vector<double> rotationMatrix{rx * rx * (1. - c) + c, rx * ry * (1. - c) - rz * s, rx * rz * (1. - c) + ry * s,
                                       rx * ry * (1. - c) + rz * s, ry * ry * (1. - c) + c, ry * rz * (1. - c) - rx * s,
                                       rx * rz * (1. - c) - ry * s, ry * rz * (1. - c) + rx * s, rz * rz * (1. - c) + c};
    return rotationMatrix;
}

