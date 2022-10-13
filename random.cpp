/*
 * random.cpp
 *
 *  Created on: 3 oct. 2022
 *      Author: romainsimon
 */

#include <stdlib.h>
#include <random>

double randomDoubleGenerator(double min, double max)
{
	// This function returns a random double in the interval [min, max).
	std::random_device rd;									// Used to obtain the seed for the following random number generator
    std::mt19937 mt(rd());  								// Mersenne Twister 19937 generator: seeded by rd(). mt is a unsigned int.
    std::uniform_real_distribution<double> dist(min, max);
	return dist(mt);										// This allows to transform the unisgned int mt into a random double between [min, max).
}

int randomIntGenerator(int min, int max)
{
	// This function returns a random int from the {min, min + 1, ... , max  - 1, max} set.
	return randomDoubleGenerator(static_cast<int>(min), static_cast<int>(max) + 1);
}
