/*
 * random.cpp
 *
 *  Created on: 3 oct. 2022
 *      Author: romainsimon
 */
#include <stdlib.h>
#include <random>
#include <iostream>



double randomDoubleGenerator(const double& min, const double& max)
{
	// This function returns a random double in the interval [min, max).
	// For Fixed seed comment these two next lines and un comment the 3d !!
    // static std::random_device rd;									// Used to obtain the seed for the following random number generator
	// static std::mt19937 mt(rd());
    static std::mt19937 mt;
    std::uniform_real_distribution<double> dist(min, max);          // The period of this generator is 2^9337-1. This there shouldn
    return dist(mt);												// This allows to transform the unisgned int mt into a random double between [min, max).
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
